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
from facilitybase import DiagClassA, CareTier, DiagClassA
from facilitybase import NURSINGQueue, VENTQueue, SKILNRSQueue
from facilitybase import PatientOverallHealth, Facility, PatientAgent, ForcedStateWard
from facilitybase import PatientStatusSetter, FacilityManager, tierToQueueMap
from hospital import estimateWork as hospitalEstimateWork
from hospital import ClassASetter, OverallHealthSetter

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

        totDsch = float(descr['totalDischarges']['value'])
        totTO = sum([elt['count']['value'] for elt in descr['totalTransfersOut']])
        ttoD = {elt['category']: elt['count']['value'] for elt in descr['totalTransfersOut']}
        fracNotTransferred = (float(totDsch) - float(totTO)) / float(totDsch)
        assert fracNotTransferred >= 0.0, '%s has more transfers than discharges' % self.name

        self.lclRates = {}
        self.lclRates['death'] = _c['deathRate']['value']
        assert fracNotTransferred >= self.lclRates['death'], '%s has more deaths than non-transfers'
        self.lclRates['home'] = fracNotTransferred - self.lclRates['death']

        expectedKeys = ['HOSPITAL', 'LTAC', 'NURSINGHOME', 'VSNF']
        if totTO:
            for key in ttoD.keys():
                assert key in expectedKeys, \
                    '%s records outgoing transfers to unexpected category %s' % (self.name, key)
            for key in expectedKeys:
                nTransfers = ttoD[key] if key in ttoD else 0
                self.lclRates[key.lower()] = float(nTransfers) / float(totDsch)
        else:
            for key in expectedKeys:
                self.lclRates[key.lower()] = 0.0
            logger.warning('%s has no transfers out', self.name)

        for subCat, cat, key in [('icu', 'hospital', 'hospTransferToICURate'),
                                 ('vent', 'vsnf', 'vsnfTransferToVentRate'),
                                 ('skilnrs', 'vsnf', 'vsnfTransferToSkilNrsRate')]:
            self.lclRates[subCat] = _c[key]['value'] * self.lclRates[cat]
            self.lclRates[cat] -= self.lclRates[subCat]

        # VSNF-specific care tiers have been separated out, so any remaining vsnf fraction
        # is actually just nursing.
        self.lclRates['nursinghome'] += self.lclRates['vsnf']
        self.lclRates['vsnf'] = 0.0

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

    def getOrderedCandidateFacList(self, oldTier, newTier, timeNow):
        """Specialized to prioritize transfers to our own wards if possible"""
        facAddrList = super(VentSNF, self).getOrderedCandidateFacList(oldTier, newTier,
                                                                      timeNow)
        if (newTier in [CareTier.NURSING, CareTier.SKILNRS, CareTier.VENT]
            and oldTier != newTier):
            queueClass = tierToQueueMap[newTier]
            lclList = [rQ for rQ in self.reqQueues if isinstance(rQ, queueClass)]
            lclList = [q.getGblAddr() for q in lclList]
            facAddrList = lclList + facAddrList
        return facAddrList

    def getStatusChangeTree(self, patientStatus, ward, treatment, startTime, timeNow):
        careTier = ward.tier
        assert careTier in [CareTier.NURSING, CareTier.SKILNRS, CareTier.VENT],\
            "VSNFs do not provide this CareTier: found %s" % careTier
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        _c = _constants
        #
        # This implementation applies the same LOS distributions regardless of CareTier,
        # so we can ignore the tier in most of what follows
        #
        if key in self.treeCache:
            return self.treeCache[key]
        else:
            adverseProb = (self.lclRates['death']
                           + self.lclRates['hospital']
                           + self.lclRates['icu']
                           + self.lclRates['ltac']
                           + self.lclRates['nursinghome']
                           + self.lclRates['vsnf']
                           + self.lclRates['vent']
                           + self.lclRates['skilnrs'])
            if adverseProb > 0.0:
                adverseTree = BayesTree.fromLinearCDF([(self.lclRates['death']
                                                        / adverseProb,
                                                        ClassASetter(DiagClassA.DEATH)),
                                                       (self.lclRates['nursinghome'] / adverseProb,
                                                        ClassASetter(DiagClassA.NEEDSREHAB)),
                                                       (self.lclRates['skilnrs'] / adverseProb,
                                                        ClassASetter(DiagClassA.NEEDSSKILNRS)),
                                                       (self.lclRates['vent'] / adverseProb,
                                                        ClassASetter(DiagClassA.NEEDSVENT)),
                                                       (self.lclRates['ltac'] / adverseProb,
                                                        ClassASetter(DiagClassA.NEEDSLTAC)),
                                                       (self.lclRates['hospital'] / adverseProb,
                                                        ClassASetter(DiagClassA.SICK)),
                                                       (self.lclRates['icu'] / adverseProb,
                                                        ClassASetter(DiagClassA.VERYSICK))],
                                                      tag='FATE')
                tree = BayesTree(BayesTree(adverseTree,
                                           ClassASetter(DiagClassA.HEALTHY),
                                           adverseProb),
                                 PatientStatusSetter(),
                                 self.cachedCDF.intervalProb, tag='LOS')
            else:
                tree = BayesTree(ClassASetter(DiagClassA.HEALTHY),
                                 PatientStatusSetter(),
                                 self.cachedCDF.intervalProb, tag='LOS')
            self.treeCache[key] = tree
            return tree

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
            if random() <= residentFrac:
                a._status = a._status._replace(overall=PatientOverallHealth.FRAIL)
            ward.lock(a)
            fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                                  fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                                  0)
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
