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
import math
from scipy.stats import lognorm
import logging
from random import random

import pyrheabase
import pyrheautils
from facilitybase import DiagClassA, CareTier, TreatmentProtocol
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent
from facilitybase import PatientStatusSetter, HOSPQueue, ICUQueue, tierToQueueMap
from facilitybase import FacilityManager
from stats import CachedCDFGenerator, BayesTree
import schemautils

category = 'HOSPITAL'
_schema = 'hospitalfacts_schema.yaml'
_constants_values = '$(MODELDIR)/constants/hospital_constants.yaml'
_constants_schema = 'hospital_ChicagoLand_constants_schema.yaml'
_validator = None
_constants = None

logger = logging.getLogger(__name__)


class ClassASetter(PatientStatusSetter):
    def __init__(self, newClassA):
        self.newClassA = newClassA

    def set(self, patientStatus, timeNow):
        return (patientStatus._replace(diagClassA=self.newClassA, startDateA=timeNow)
                ._replace(relocateFlag=True))
    
    def __str__(self):
        return 'PatientStatusSetter(classA <- %s)' % DiagClassA.names[self.newClassA]

    def __repr__(self):
        return 'PatientStatusSetter(classA <- %s)' % DiagClassA.names[self.newClassA]


class OverallHealthSetter(PatientStatusSetter):
    def __init__(self, newOverallHealth):
        self.newOverallHealth = newOverallHealth

    def set(self, patientStatus, timeNow):
        return patientStatus._replace(overall=self.newOverallHealth, startDateA=timeNow)

    def __str__(self):
        return ('PatientStatusSetter(overallHealth <- %s)'
                % PatientOverallHealth.names[self.newOverallHealth])


class HospitalManager(FacilityManager):
    def allocateAvailableBed(self, requestedTier, strict=False):
        """
        If a HOSP tier bed is requested, maybe upgrade to ICU to simulate patients who
        are triaged to the ICU on arrival.
        """
        if strict:
            tier = requestedTier
        else:
            if requestedTier == CareTier.HOSP and random() < _constants['fracTriageHOSPToICU']['value']:
                tier = CareTier.ICU
            else:
                tier = requestedTier
        return super(HospitalManager, self).allocateAvailableBed(tier)


class Hospital(Facility):
    def __init__(self, descr, patch, policyClasses=None, categoryNameMapper=None):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr,
                          descr, patch,
                          reqQueueClasses=[HOSPQueue, ICUQueue],
                          managerClass=HospitalManager,
                          policyClasses=policyClasses,
                          categoryNameMapper=categoryNameMapper)
        descr = self.mapDescrFields(descr)
        bedsPerWard = _constants['bedsPerWard']['value']
        bedsPerICUWard = _constants['bedsPerWard']['value']
        
        _c = _constants
        totDsch = float(descr['totalDischarges']['value'])
        totTO = sum([elt['count']['value'] for elt in descr['totalTransfersOut']])
        ttoD = {elt['category']: elt['count']['value'] for elt in descr['totalTransfersOut']}
        fracNotTransferred = (float(totDsch) - float(totTO)) / float(totDsch)
        assert fracNotTransferred >= 0.0, '%s has more transfers than discharges' % self.name
        
        self.lclRates = {}
        self.lclRates['death'] = _c['hospDischargeViaDeathFrac']['value']
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

        # All patients getting rehab must be in the transfer streams to rehab-providing facilities
#         fracWhoMayHaveRehab = self.lclRates['nursinghome'] + self.lclRates['vsnf']
#         assert fracWhoMayHaveRehab >= _c['fracOfDischargesRequiringRehab']['value'], \
#             ("%s has only %s of discharges going to rehab facilities but should rehab %s of patients" %
#              (self.name, fracWhoMayHaveRehab, _c['fracOfDischargesRequiringRehab']['value']))

        self.icuDischargeViaDeathFrac = _constants['icuDischargeViaDeathFrac']['value']
        if 'nBeds' in descr:
            allBeds = descr['nBeds']['value']
            if 'nBedsICU' in descr:
                icuBeds = descr['nBedsICU']['value']
            else:
                icuBeds = descr['fracAdultPatientDaysICU']['value'] * allBeds
            icuWards = int(icuBeds/bedsPerICUWard) + 1
            nonICUBeds = max(allBeds - icuBeds, 0)
            nonICUWards = int(float(nonICUBeds)/bedsPerWard) + 1
        else:
            meanPop = descr['meanPop']['value']
            meanICUPop = meanPop * descr['fracAdultPatientDaysICU']['value']
            meanNonICUPop = meanPop - meanICUPop
            icuWards = int(math.ceil(meanICUPop / bedsPerICUWard))
            nonICUWards = int(math.ceil(meanNonICUPop / bedsPerWard))

        icuBeds = icuWards * bedsPerICUWard
        nonICUBeds = nonICUWards * bedsPerWard

        assert descr['losModel']['pdf'] == 'lognorm(mu=$0,sigma=$1)', \
            'Hospital %(abbrev)s LOS PDF is not lognorm(mu=$0,sigma=$1)' % descr
        self.hospCachedCDF = CachedCDFGenerator(lognorm(descr['losModel']['parms'][1],
                                                        scale=math.exp(descr['losModel']
                                                                       ['parms'][0])))
        self.hospTreeCache = {}
        if descr['meanLOSICU']['value'] == 0.0:
            self.icuCachedCDF = None
        else:
            mean = descr['meanLOSICU']['value']
            sigma = _constants['icuLOSLogNormSigma']['value']
            mu = math.log(mean) - (0.5 * sigma * sigma)
            self.icuCachedCDF = CachedCDFGenerator(lognorm(sigma, scale=math.exp(mu)))
        self.icuTreeCache = {}

        for i in xrange(icuWards):
            self.addWard(Ward(('%s_%s_%s_%s_%d' %
                               (category, patch.name, descr['abbrev'], 'ICU', i)),
                              patch, CareTier.ICU, bedsPerICUWard))
        for i in xrange(nonICUWards):
            self.addWard(Ward(('%s_%s_%s_%s_%d' %
                               (category, patch.name, descr['abbrev'], 'HOSP', i)),
                              patch, CareTier.HOSP, bedsPerWard))

    def getOrderedCandidateFacList(self, oldTier, newTier, timeNow):
        """Specialized to restrict transfers to being between our own HOSP and ICU"""
        queueClass = tierToQueueMap[newTier]
        if ((oldTier == CareTier.ICU and newTier == CareTier.HOSP)
                or (oldTier == CareTier.HOSP and newTier == CareTier.ICU)):
            qList = [rQ for rQ in self.reqQueues if isinstance(rQ, queueClass)]
#             print '%s: clause 1' % self.abbrev
            return [q.getGblAddr() for q in qList]
        else:
            facAddrList = super(Hospital, self).getOrderedCandidateFacList(oldTier, newTier,
                                                                           timeNow)
#             print '%s: clause 2: %s %s' % (self.abbrev, CareTier.names[newTier], facAddrList[:3])
        return facAddrList

    def getStatusChangeTree(self, patientStatus, ward, treatment, startTime, timeNow):
        assert treatment == TreatmentProtocol.NORMAL, \
            ("Hospitals only offer 'NORMAL' treatment; found %s" %
             TreatmentProtocol.names[treatment])
        careTier = ward.tier
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        _c = _constants
        if careTier == CareTier.HOSP:
            if key in self.hospTreeCache:
                return self.hospTreeCache[key]
            else:
                changeProb = self.hospCachedCDF.intervalProb(*key)
                changeTree = BayesTree.fromLinearCDF([(self.lclRates['death'],
                                                       ClassASetter(DiagClassA.DEATH)),
                                                      (self.lclRates['hospital'],
                                                       ClassASetter(DiagClassA.SICK)),
                                                      (self.lclRates['ltac'],
                                                       ClassASetter(DiagClassA.NEEDSLTAC)),
                                                      (self.lclRates['nursinghome'],
                                                       ClassASetter(DiagClassA.NEEDSREHAB)),
                                                      (self.lclRates['vsnf'],
                                                       ClassASetter(DiagClassA.NEEDSSKILNRS)),
                                                      (self.lclRates['home'],
                                                       ClassASetter(DiagClassA.HEALTHY)),
                                                      ])

                tree = BayesTree(changeTree,
                                 PatientStatusSetter(),
                                 changeProb)
                self.hospTreeCache[key] = tree
                return tree
        elif careTier == CareTier.ICU:
            if key in self.icuTreeCache:
                return self.icuTreeCache[key]
            else:
                changeProb = self.icuCachedCDF.intervalProb(*key)
                tree = BayesTree(BayesTree.fromLinearCDF([(self.icuDischargeViaDeathFrac,
                                                           ClassASetter(DiagClassA.DEATH)),
                                                          ((1.0 - self.icuDischargeViaDeathFrac),
                                                           ClassASetter(DiagClassA.SICK))]),
                                 PatientStatusSetter(),
                                 changeProb)
                self.icuTreeCache[key] = tree
                return tree
        else:
            raise RuntimeError('Hospitals do not provide care tier %s' % careTier)

def _populate(fac, descr, patch):
    assert 'meanPop' in descr, \
        "Hospital description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = float(descr['meanPop']['value'])
    if 'nBeds' in descr and meanPop > descr['nBeds']['value']:
        logger.warning('Hospital %s meanPop %s > nBeds %s'
                       % (descr['abbrev'], meanPop, descr['nBeds']['value']))
        meanPop = descr['nBeds']['value']
    meanICUPop = meanPop * descr['fracAdultPatientDaysICU']['value']
    meanHospPop = meanPop - meanICUPop
    agentList = []
    for i in xrange(int(round(meanICUPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.ICU, strict=True)
        assert ward is not None, 'Ran out of ICU beds populating %(abbrev)s!' % descr
        a = PatientAgent('PatientAgent_ICU_%s_%d' % (ward._name, i), patch, ward)
        ward.lock(a)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              0)
        agentList.append(a)
    for i in xrange(int(round(meanHospPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.HOSP, strict=True)
        assert ward is not None, 'Ran out of HOSP beds populating %(abbrev)s!' % descr
        a = PatientAgent('PatientAgent_HOSP_%s_%d' % (ward._name, i), patch, ward)
        ward.lock(a)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              0)
        agentList.append(a)
    return agentList


def generateFull(facilityDescr, patch, policyClasses=None, categoryNameMapper=None):
    fac = Hospital(facilityDescr, patch, policyClasses=policyClasses,
                   categoryNameMapper=categoryNameMapper)
    return [fac], fac.getWards(), _populate(fac, facilityDescr, patch)


def estimateWork(facRec):
    if 'meanPop' in facRec:
        return facRec['meanPop']['value']
    elif 'nBeds' in facRec:
        return facRec['nBeds']['value']
    else:
        logger.warn('Cannot estimate work for %(abbrev)s' % facRec)
        return 0


def checkSchema(facilityDescr):
    global _validator
    if _validator is None:
        _validator = schemautils.getValidator(_schema)
    nErrors = sum([1 for e in _validator.iter_errors(facilityDescr)])  # @UnusedVariable
    return nErrors


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(_constants_values, _constants_schema)
