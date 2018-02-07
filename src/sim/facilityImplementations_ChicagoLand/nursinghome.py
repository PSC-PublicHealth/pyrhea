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

_rhea_svn_id_ = "$Id$"

import os.path
from random import random
import math
from scipy.stats import lognorm, expon
import logging

import pyrheabase
import pyrheautils
import schemautils
from stats import CachedCDFGenerator, lognormplusexp, BayesTree
from facilitybase import DiagClassA, CareTier, TreatmentProtocol, NURSINGQueue
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent
from facilitybase import PatientStatusSetter, buildTimeTupleList
from hospital import estimateWork as hospitalEstimateWork
from hospital import ClassASetter, OverallHealthSetter

logger = logging.getLogger(__name__)

category = 'NURSINGHOME'
_schema = 'nursinghomefacts_ChicagoLand_schema.yaml'
_constants_values = '$(MODELDIR)/constants/nursinghome_constants.yaml'
_constants_schema = 'nursinghome_ChicagoLand_constants_schema.yaml'
_constants = None
_validator = None


class NursingWard(Ward):
    pass
#     def handlePatientArrival(self, patientAgent, timeNow):
#         super(NursingWard, self).handlePatientArrival(patientAgent, timeNow)
#         # Sometimes the transfer tables send a patient in good health- they must
#         # be here for rehab.
#         if (patientAgent.getStatus().diagClassA == DiagClassA.WELL
#                 and not patientAgent.getTreatment('rehab')):
#             print '######## %s fixed %s %s at %s' % (self._name, patientAgent.name,
#                                                buildTimeTupleList(patientAgent, timeNow), timeNow)
#             patientAgent.setTreatment(rehab=True)

class NursingHome(Facility):
    def __init__(self, descr, patch, policyClasses=None, categoryNameMapper=None):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr,
                          descr, patch,
                          reqQueueClasses=[NURSINGQueue],
                          policyClasses=policyClasses,
                          categoryNameMapper=categoryNameMapper)
        descr = self.mapDescrFields(descr)
        _c = _constants
        nBeds = int(descr['nBeds']['value'])
        losModel = descr['losModel']
        assert losModel['pdf'] == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)', \
            "Unexpected losModel form %s for %s!" % (losModel['pdf'], descr['abbrev'])

        self.initialResidentFrac = (1.0 - losModel['parms'][0])
        self.initialUnhealthyFrac = _c['initialUnhealthyFrac']['value'] # frac of non-residents

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


        lMP = losModel['parms']
        self.rehabCachedCDF = CachedCDFGenerator(lognorm(lMP[2],
                                                         scale=math.exp(lMP[1])))
        self.rehabTreeCache = {}
        self.frailCachedCDF = CachedCDFGenerator(expon(scale=1.0/lMP[3]))
        self.frailTreeCache = {}
        self.frailRehabTreeCache = {}
        self.addWard(NursingWard('%s_%s_%s' % (category, patch.name, descr['abbrev']),
                                 patch, CareTier.NURSING, nBeds))

    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        self.rehabTreeCache = {}        
        self.frailRehabTreeCache = {}        
        self.frailTreeCache = {}        

    def getStatusChangeTree(self, patientStatus, ward, treatment, startTime, timeNow):
        careTier = ward.tier
        assert careTier == CareTier.NURSING, \
            "Nursing homes only offer CareTier 'NURSING'; found %s" % careTier
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        _c = _constants
        if patientStatus.overall == PatientOverallHealth.FRAIL:
            if treatment.rehab:
                if key in self.frailRehabTreeCache:
                    return self.frailRehabTreeCache[key]
                else:
                    # If no major changes happen, allow for the patient's completing rehab
                    innerTree = BayesTree(ClassASetter(DiagClassA.WELL),
                                          PatientStatusSetter(),
                                          self.rehabCachedCDF.intervalProb(*key))
            else:
                if key in self.frailTreeCache:
                    return self.frailTreeCache[key]
                else:
                    innerTree = PatientStatusSetter()  # No change of state
            unhealthySetter = OverallHealthSetter(PatientOverallHealth.UNHEALTHY)
            changeTree = BayesTree.fromLinearCDF([(self.lclRates['death'],
                                                   ClassASetter(DiagClassA.DEATH)),
                                                  (self.lclRates['nursinghome'],
                                                   ClassASetter(DiagClassA.NEEDSREHAB)),
                                                  (self.lclRates['vent'],
                                                   ClassASetter(DiagClassA.NEEDSVENT)),
                                                  (self.lclRates['skilnrs'],
                                                   ClassASetter(DiagClassA.NEEDSSKILNRS)),
                                                  (self.lclRates['ltac'],
                                                   ClassASetter(DiagClassA.NEEDSLTAC)),
                                                  (self.lclRates['hospital'],
                                                   ClassASetter(DiagClassA.SICK)),
                                                  (self.lclRates['icu'],
                                                   ClassASetter(DiagClassA.VERYSICK)),
                                                  (self.lclRates['home'],
                                                   unhealthySetter)], tag='FATE')

            tree = BayesTree(changeTree, innerTree,
                             self.frailCachedCDF.intervalProb, tag='LOS')
            if treatment.rehab:
                self.frailRehabTreeCache[key] = tree
            else:
                self.frailTreeCache[key] = tree
            return tree
        else:  # overall non-FRAIL
            if treatment.rehab:
                if key in self.rehabTreeCache:
                    return self.rehabTreeCache[key]
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
                                                               (self.lclRates['nursinghome']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.NEEDSREHAB)),
                                                               (self.lclRates['vent']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.NEEDSVENT)),
                                                               (self.lclRates['skilnrs']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.NEEDSSKILNRS)),
                                                               (self.lclRates['ltac']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.NEEDSLTAC)),
                                                               (self.lclRates['hospital']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.SICK)),
                                                               (self.lclRates['icu']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.VERYSICK))],
                                                              tag='FATE')
                        tree = BayesTree(BayesTree(adverseTree,
                                                   ClassASetter(DiagClassA.WELL),
                                                   adverseProb),
                                         PatientStatusSetter(),
                                         self.rehabCachedCDF.intervalProb, tag='LOS')
                    else:
                        tree = BayesTree(ClassASetter(DiagClassA.WELL),
                                         PatientStatusSetter(),
                                         self.rehabCachedCDF.intervalProb, tag='LOS')
                    self.rehabTreeCache[key] = tree
                    return tree
            else:  # healthy and not rehab
                if patientStatus.diagClassA in ([DiagClassA.NEEDSLTAC, DiagClassA.SICK,
                                                 DiagClassA.VERYSICK, DiagClassA.NEEDSVENT,
                                                 DiagClassA.NEEDSSKILNRS, DiagClassA.WELL]):
                    logger.warning('fac %s status: %s careTier %s startTime: %s: '
                                   'this patient should be gone by now'
                                   % (self.name, str(patientStatus), CareTier.names[careTier],
                                      startTime))
                    return BayesTree(PatientStatusSetter())
                else:
                    raise RuntimeError('Patients with non-FRAIL overall health should only be'
                                       ' in NURSING care for rehab')

    def getInitialOverallHealth(self, ward, timeNow):  # @UnusedVariable
        if random() <= self.initialResidentFrac:
            return PatientOverallHealth.FRAIL
        else:
            if random() <= self.initialUnhealthyFrac:
                return PatientOverallHealth.UNHEALTHY
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
    for i in xrange(int(round(meanPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.NURSING)
        assert ward is not None, 'Ran out of beds populating %(abbrev)s!' % descr
        a = PatientAgent('PatientAgent_NURSING_%s_%d' % (ward._name, i), patch, ward)
        if a.getStatus().overall == PatientOverallHealth.FRAIL:
            a.setStatus(homeAddr=ward.getGblAddr())  # They live here
        else:
            a.setTreatment(rehab=True)  # They must be here for rehab
        ward.lock(a)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              None)
        agentList.append(a)
    return agentList


def generateFull(facilityDescr, patch, policyClasses=None, categoryNameMapper=None):
    fac = NursingHome(facilityDescr, patch, policyClasses=policyClasses,
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
