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
import math
from scipy.stats import lognorm
import logging
from random import shuffle

import pyrheabase
import pyrheautils
from facilitybase import DiagClassA, CareTier, TreatmentProtocol
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent
from facilitybase import PatientStatusSetter, HOSPQueue, ICUQueue, tierToQueueMap
from stats import CachedCDFGenerator, BayesTree
import schemautils

category = 'HOSPITAL'
_schema = 'hospitalfacts_schema.yaml'
_constants_values = '$(MODELDIR)/constants/hospital_constants.yaml'
_constants_schema = 'hospital_constants_schema.yaml'
_validator = None
_constants = None

logger = logging.getLogger(__name__)


class ClassASetter(PatientStatusSetter):
    def __init__(self, newClassA):
        self.newClassA = newClassA

    def set(self, patientStatus, timeNow):
        return patientStatus._replace(diagClassA=self.newClassA, startDateA=timeNow)

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


class Hospital(Facility):
    def __init__(self, descr, patch, policyClasses=None, categoryNameMapper=None):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr,
                          descr, patch,
                          reqQueueClasses=[HOSPQueue, ICUQueue],
                          policyClasses=policyClasses,
                          categoryNameMapper=categoryNameMapper)
        descr = self.mapDescrFields(descr)
        bedsPerWard = _constants['bedsPerWard']['value']
        bedsPerICUWard = _constants['bedsPerWard']['value']
        self.hospDischargeViaDeathFrac = _constants['hospDischargeViaDeathFrac']['value']
        nTotalTransfersOut = sum([v['count']['value'] for v in descr['totalTransfersOut']])
        self.hospTotalTransferFrac = (float(nTotalTransfersOut)
                                      / float(descr['totalDischarges']['value']))
        self.fracDischargeHealthy = 1.0 - (self.hospTotalTransferFrac
                                           + self.hospDischargeViaDeathFrac)
        transferFrac = {}
        for v in descr['totalTransfersOut']:
            implCat = v['category']
            if implCat not in transferFrac:
                transferFrac[implCat] = 0.0
            transferFrac[implCat] += (self.hospTotalTransferFrac * float(v['count']['value'])
                                      / float(nTotalTransfersOut))
        self.fracTransferHosp = (transferFrac['HOSPITAL'] if 'HOSPITAL' in transferFrac else 0.0)
        self.fracTransferLTAC = (transferFrac['LTAC'] if 'LTAC' in transferFrac else 0.0)
        self.fracTransferNH = (transferFrac['NURSINGHOME'] if 'NURSINGHOME' in transferFrac
                               else 0.0)
        self.icuDischargeViaDeathFrac = _constants['icuDischargeViaDeathFrac']['value']
        if 'nBeds' in descr:
            icuBeds = descr['fracAdultPatientDaysICU']['value'] * descr['nBeds']['value']
            icuWards = int(icuBeds/bedsPerICUWard) + 1
            nonICUBeds = max(descr['nBeds']['value'] - icuBeds, 0)
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

    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        self.hospTreeCache = {}
        self.icuTreeCache = {}

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
        assert not treatment.rehab, \
            ("Hospitals do not provide rehab treatment; found %s" % treatment)
        careTier = ward.tier
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        _c = _constants
        if careTier == CareTier.HOSP:
            if key in self.hospTreeCache:
                return self.hospTreeCache[key]
            else:
                changeProb = self.hospCachedCDF.intervalProb(*key)
                rehabFrac = _c['fracOfDischargesRequiringRehab']['value']
                changeTree = BayesTree.fromLinearCDF([(self.hospDischargeViaDeathFrac,
                                                       ClassASetter(DiagClassA.DEATH)),
                                                      (self.fracTransferHosp,
                                                       ClassASetter(DiagClassA.SICK)),
                                                      (self.fracTransferLTAC,
                                                       ClassASetter(DiagClassA.NEEDSLTAC)),
                                                      (((self.fracTransferNH
                                                         + self.fracDischargeHealthy)
                                                        * rehabFrac),
                                                       ClassASetter(DiagClassA.NEEDSREHAB)),
                                                      (((self.fracTransferNH
                                                         + self.fracDischargeHealthy)
                                                        * (1.0 - rehabFrac)),
                                                       ClassASetter(DiagClassA.HEALTHY))])
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

    def prescribe(self, ward, patientId, patientDiagnosis, patientTreatment):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
                return (CareTier.HOME, patientTreatment._replace(rehab=False))
            else:
                return (CareTier.NURSING, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
            return (CareTier.NURSING, patientTreatment._replace(rehab=True))
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSLTAC:
            return (CareTier.LTAC, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.SICK:
            return (CareTier.HOSP, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.VERYSICK:
            return (CareTier.ICU, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.DEATH:
            return (None, patientTreatment._replace(rehab=False))
        else:
            raise RuntimeError('Unknown DiagClassA %s' % str(patientDiagnosis.diagClassA))


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
        ward = fac.manager.allocateAvailableBed(CareTier.ICU)
        assert ward is not None, 'Ran out of ICU beds populating %(abbrev)s!' % descr
        a = PatientAgent('PatientAgent_ICU_%s_%d' % (ward._name, i), patch, ward)
        ward.lock(a)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              0)
        agentList.append(a)
    for i in xrange(int(round(meanHospPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.HOSP)
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
_constants = pyrheautils.importConstants(_constants_values,
                                         _constants_schema)
