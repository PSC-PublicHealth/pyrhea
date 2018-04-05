#! /usr/bin/env python

###################################################################################
# Copyright   2015, Pittsburgh Supercomputing Center (PSC).  All Rights Reserved. #
# =============================================================================== #
#                                                                                 #
# Permission to use, copy, and modify this software and its documentation without #
# fee for personal use within your organization is hereby granted, provided that  #
# the above copyright notice is preserved in all copies and that the copyright    #
# and this permission notice appear in supporting documentation.  All other       #
# restrictions and obligations are defined in the GNU Affero General Public       #
# License v3 (AGPL-3.0) located at http://www.gnu.org/licenses/agpl-3.0.html  A   #
# copy of the license is also provided in the top level of the source directory,  #
# in the file LICENSE.txt.                                                        #
#                                                                                 #
###################################################################################

import os.path
from random import random
import math
from scipy.stats import lognorm, expon
import logging

import pyrheabase
import pyrheautils
import schemautils
from stats import CachedCDFGenerator, lognormplusexp, BayesTree
from facilitybase import CareTier, DiagClassA, PthStatus
from facilitybase import NURSINGQueue, VENTQueue, SKILNRSQueue
from facilitybase import PatientOverallHealth, Facility, PatientAgent, ForcedStateWard
from facilitybase import PatientStatusSetter, FacilityManager, tierToQueueMap
from hospital import estimateWork as hospitalEstimateWork
from hospital import ClassASetter, OverallHealthSetter
from hospital import buildChangeTree, biasTransfers

logger = logging.getLogger(__name__)

category = 'VSNF'
_schema = 'vsnffacts_ChicagoLand_schema.yaml'
_constants_values = '$(MODELDIR)/constants/vsnf_constants.yaml'
_constants_schema = 'vsnf_ChicagoLand_constants_schema.yaml'
_constants = None
_validator = None


class NursingWard(ForcedStateWard):
    def __init__(self, name, patch, nBeds):
        super(NursingWard, self).__init__(name, patch, CareTier.NURSING, nBeds)
        

class SkilNrsWard(ForcedStateWard):
    def __init__(self, name, patch, nBeds):
        super(SkilNrsWard, self).__init__(name, patch, CareTier.SKILNRS, nBeds)


class VentWard(ForcedStateWard):
    def __init__(self, name, patch, nBeds):
        super(VentWard, self).__init__(name, patch, CareTier.VENT, nBeds)


class VentSNF(Facility):
    def __init__(self, descr, patch, policyClasses=None, categoryNameMapper=None):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr,
                          descr, patch,
                          reqQueueClasses=[SKILNRSQueue, NURSINGQueue, VENTQueue],
                          policyClasses=policyClasses,
                          categoryNameMapper=categoryNameMapper)
        descr = self.mapDescrFields(descr)
        _c = _constants

        losModel = descr['losModel']
        assert losModel['pdf'] == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)', \
            "Unexpected losModel form %s for %s!" % (losModel['pdf'], descr['abbrev'])

        self.initialResidentFrac = (1.0 - losModel['parms'][0])
        self.initialUnhealthyFracByTier = {}
        self.initialNonResidentFrailFracByTier = {}
        nameToTierD = {val:key for key,val in CareTier.names.items()}
        for ent in _c['initialUnhealthyFracByTier']:
            self.initialUnhealthyFracByTier[nameToTierD[ent['tier']]] = ent['frac']['value']
        for ent in _c['initialNonResidentFrailFracByTier']:
            self.initialNonResidentFrailFracByTier[nameToTierD[ent['tier']]] = ent['frac']['value']

        totDsch = float(descr['totalDischarges']['value'])
        totTO = sum([elt['count']['value'] for elt in descr['totalTransfersOut']])
        ttoD = {elt['category']: elt['count']['value'] for elt in descr['totalTransfersOut']}
        fracNotTransferred = (float(totDsch) - float(totTO)) / float(totDsch)
        assert fracNotTransferred >= 0.0, '%s has more transfers than discharges' % self.name

        lclRates = {}
        lclRates['death'] = _c['deathRate']['value']
        assert fracNotTransferred >= lclRates['death'], '%s has more deaths than non-transfers'
        lclRates['home'] = fracNotTransferred - lclRates['death']

        expectedKeys = ['HOSPITAL', 'LTAC', 'NURSINGHOME', 'VSNF']
        if totTO:
            for key in ttoD.keys():
                assert key in expectedKeys, \
                    '%s records outgoing transfers to unexpected category %s' % (self.name, key)
            for key in expectedKeys:
                nTransfers = ttoD[key] if key in ttoD else 0
                lclRates[key.lower()] = float(nTransfers) / float(totDsch)
        else:
            for key in expectedKeys:
                lclRates[key.lower()] = 0.0
            logger.warning('%s has no transfers out', self.name)

        for subCat, cat, key in [('icu', 'hospital', 'hospTransferToICURate'),
                                 ('vent', 'vsnf', 'vsnfTransferToVentRate'),
                                 ('skilnrs', 'vsnf', 'vsnfTransferToSkilNrsRate')]:
            lclRates[subCat] = _c[key]['value'] * lclRates[cat]
            lclRates[cat] -= lclRates[subCat]

        # VSNF-specific care tiers have been separated out, so any remaining vsnf fraction
        # is actually just nursing.
        lclRates['nursinghome'] += lclRates['vsnf']
        lclRates['vsnf'] = 0.0

        # The second element of the tuples built below is the pthRates for each tier.
        # pthRates is the equivalent of lclRates but for pathogen carriers.  It must
        # be initialized lazily because the pathogen has not yet been defined.
        self.rateD = {}
        for tier in [CareTier.NURSING, CareTier.SKILNRS, CareTier.VENT]:
            self.rateD[tier] = (lclRates.copy(), None)

        nBeds = int(descr['nBeds']['value'])
        nBedsVent = int(nBeds * descr['fracVentBeds']['value'])
        nBedsSkil = int(nBeds * descr['fracSkilledBeds']['value'])
        nBedsOther = nBeds - (nBedsVent + nBedsSkil)
        self.addWard(NursingWard('%s_%s_%s' % (category, patch.name, descr['abbrev']),
                                 patch, nBedsOther))
        self.addWard(SkilNrsWard('%s_%s_%s' % (category, patch.name, descr['abbrev']),
                                 patch, nBedsSkil))
        self.addWard(VentWard('%s_%s_%s' % (category, patch.name, descr['abbrev']),
                              patch, nBedsVent))

        k, mu, sigma, lmda = losModel['parms']
        self.cachedCDF = CachedCDFGenerator(lognormplusexp(s=sigma, mu=mu, k=k, lmda=lmda))
        self.treeCache = {}

    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        self.treeCache = {}
        for tier, tpl in self.rateD.items():
            lclRates, pthRates = tpl
            self.rateD[tier] = (lclRates, None)  # force recalculation of biases

    def getOrderedCandidateFacList(self, patientAgent, oldTier, newTier, timeNow):
        """Specialized to prioritize transfers to our own wards if possible"""
        facAddrList = super(VentSNF, self).getOrderedCandidateFacList(patientAgent, oldTier, newTier,
                                                                      timeNow)
        if (newTier in [CareTier.NURSING, CareTier.SKILNRS, CareTier.VENT]
            and oldTier != newTier):
            queueClass = tierToQueueMap[newTier]
            lclList = [rQ for rQ in self.reqQueues if isinstance(rQ, queueClass)]
            lclList = [q.getGblAddr() for q in lclList]
            facAddrList = lclList + facAddrList
        return facAddrList

    def getStatusChangeTree(self, patientAgent, startTime, timeNow):
        patientStatus = patientAgent.getStatus()
        ward = patientAgent.ward
        treatment = patientAgent.getTreatmentProtocol()
        careTier = ward.tier
        dCA = patientAgent.getDiagnosis().diagClassA
        assert careTier in [CareTier.NURSING, CareTier.SKILNRS, CareTier.VENT],\
            "VSNFs do not provide this CareTier: found %s" % careTier
        if not ((careTier == CareTier.SKILNRS and dCA == DiagClassA.NEEDSSKILNRS)
                or (careTier == CareTier.VENT and dCA == DiagClassA.NEEDSVENT)
                or (careTier == CareTier.NURSING and (dCA == DiagClassA.NEEDSREHAB
                                                      or patientAgent.getDiagnosis().overall
                                                      == PatientOverallHealth.FRAIL))):
            # This patient doesn't belong in this ward
            logger.warning('fac %s patient: %s careTier %s with status %s startTime: %s: '
                           'this patient should be gone by now'
                           % (self.name, patientAgent.name, CareTier.names[careTier],
                              DiagClassA.names[patientStatus.diagClassA], startTime))
            return BayesTree(PatientStatusSetter())

        lclRates, pthRates = self.rateD[careTier]
        if not pthRates:
            # Lazy initialization
            infectiousAgent = self.getWards(careTier)[0].iA
            lclRates, pthRates = biasTransfers(lclRates,
                                               infectiousAgent,
                                               CareTier.ICU,
                                               self.abbrev, self.category,
                                               careTier)
            self.rateD[careTier] = (lclRates, pthRates)

        biasFlag = (patientStatus.pthStatus not in [PthStatus.CLEAR, PthStatus.RECOVERED])
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA,
               careTier, biasFlag)
        #
        # This implementation applies the same LOS distributions regardless of CareTier,
        # so we can ignore the tier in most of what follows
        #
        if key in self.treeCache:
            return self.treeCache[key]
        else:
            if biasFlag:
                changeTree = buildChangeTree(pthRates)
            else:
                changeTree = buildChangeTree(lclRates)
            tree = BayesTree(changeTree,
                             PatientStatusSetter(),
                             self.cachedCDF.intervalProb, tag='LOS')

            self.treeCache[key] = tree
            return tree

    def getInitialOverallHealth(self, ward, timeNow):  # @UnusedVariable
        tier = ward.tier
        if random() <= self.initialResidentFrac:
            return PatientOverallHealth.FRAIL
        else:
            rn = random()
            if rn <= self.initialUnhealthyFracByTier[tier]:
                return PatientOverallHealth.UNHEALTHY
            elif rn <= (self.initialUnhealthyFracByTier[tier]
                        + self.initialNonResidentFrailFracByTier[tier]):
                return PatientOverallHealth.FRAIL
            else:
                return PatientOverallHealth.HEALTHY


def _populate(fac, descr, patch):
    assert 'meanPop' in descr, \
        "Nursing home description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = descr['meanPop']['value']
    if 'nBeds' in descr and meanPop > descr['nBeds']['value']:
        logger.warning('Nursing Home %s meanPop %s > nBeds %s'
                       % (descr['abbrev'], meanPop, descr['nBeds']['value']))
        meanPop = descr['nBeds']['value']
    losModel = descr['losModel']
    assert losModel['pdf'] == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)', \
        "Unexpected losModel form %s for %s!" % (losModel['pdf'], descr['abbrev'])
    # The following is approximate, but adequate...
    residentFrac = (1.0 - losModel['parms'][0])
    agentList = []
    nVent = int(round(descr['fracVentAdmissions']['value'] * meanPop))
    nSkil = int(round(descr['fracSkilledAdmissions']['value'] * meanPop))
    nOther = int(round(meanPop - (nVent + nSkil)))
    for count, tier in [(nVent, CareTier.VENT),
                        (nSkil, CareTier.SKILNRS),
                        (nOther, CareTier.NURSING)]:
        for i in xrange(count):
            ward = fac.manager.allocateAvailableBed(tier)
            if ward is None:
                logger.warn('Ran out of beds populating %s tier %s after %s of %s!',
                            descr['abbrev'], CareTier.names[tier], i, count)
                break
            a = PatientAgent('PatientAgent_%s_%s_%d' % (CareTier.names[tier], ward._name, i),
                             patch, ward)
            ward.lock(a)
            fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                                  fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                                  None)
            agentList.append(a)
    return agentList


def generateFull(facilityDescr, patch, policyClasses=None, categoryNameMapper=None):
    fac = VentSNF(facilityDescr, patch, policyClasses=policyClasses,
                  categoryNameMapper=categoryNameMapper)
    return [fac], fac.getWards(), _populate(fac, facilityDescr, patch)


def estimateWork(descr):
    return hospitalEstimateWork(descr)


def checkSchema(facilityDescr):
    global _validator
    if _validator is None:
        _validator = schemautils.getValidator(_schema)
    nErrors = sum([1 for e in _validator.iter_errors(facilityDescr)])  # @UnusedVariable
    return nErrors


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(_constants_values,
                                         _constants_schema)
