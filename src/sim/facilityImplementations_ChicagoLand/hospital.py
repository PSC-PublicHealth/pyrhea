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
import numpy as np
from scipy.stats import lognorm
import logging
from random import random

import pyrheabase
import pyrheautils
from facilitybase import DiagClassA, CareTier
from facilitybase import PatientOverallHealth, Facility, PatientAgent, ForcedStateWard
from facilitybase import PatientStatusSetter, ClassASetter, HOSPQueue, ICUQueue, tierToQueueMap
from facilitybase import FacilityManager, PthStatus
from facilitybase import Ward as BaseWard
from stats import CachedCDFGenerator, BayesTree
import schemautils

category = 'HOSPITAL'
_schema = 'hospitalfacts_schema.yaml'
_constants_values = '$(MODELDIR)/constants/hospital_constants.yaml'
_constants_schema = 'hospital_ChicagoLand_constants_schema.yaml'
_validator = None
_constants = None

logger = logging.getLogger(__name__)

tierKeyMap = {CareTier.HOME: 'home',
              CareTier.NURSING: 'nursinghome',
              CareTier.LTAC: 'ltac',
              CareTier.HOSP: 'hospital',
              CareTier.ICU: 'icu',
              CareTier.VENT: 'vent',
              CareTier.SKILNRS: 'skilnrs'}

knownTierKeys = tierKeyMap.values() + ['death']


def pickWardSizes(nBeds, bedsPerWard):
    """
    Given a number of beds and a nominal number of beds per ward, generate an
    iterable giving a series of ward sizes which 'comes close to' the requested
    total bed count and ward size.

    Both input parameters are assumed to be integers.
    """
    nWards = int((nBeds//bedsPerWard) + 1)
    if nWards > 1:
        delta = (bedsPerWard * nWards) - nBeds
        if delta > 0.5 * bedsPerWard:
            # Better to reduce the ward count by 1 and spread out the extra
            # beds evenly
            nWards -= 1
        else:
            # Better to have a few under-full wards
            pass
        samps = np.full(nWards, bedsPerWard)
        conv = nBeds - np.sum(samps)
        while True:
            if conv > 0:
                samps += np.random.poisson(float(conv)/ float(nWards), nWards)
            elif conv < 0:
                samps -= np.random.poisson(float(-conv)/ float(nWards), nWards)
            else:
                break
            conv = nBeds - np.sum(samps)
        return samps
    else:
        return np.asarray([nBeds])


def buildChangeTree(lclRates, forceRelocateDiag=None):
    """
    Given a dict with probabilities for the standard fates (CareTiers plus death),
    return a Bayes tree implementing it.  Note that the probabilities must sum to 1.0 .
    """
    for key in lclRates:
        if not (key in knownTierKeys or lclRates[key] == 0.0):
            raise RuntimeError('Unknown transfer rate key %s' % key)

    ratePairL = []
    for key, tier in [('death', DiagClassA.DEATH), ('hospital', DiagClassA.SICK),
                      ('icu', DiagClassA.VERYSICK), ('ltac', DiagClassA.NEEDSLTAC),
                      ('nursinghome', DiagClassA.NEEDSREHAB), ('vent', DiagClassA.NEEDSVENT),
                      ('skilnrs', DiagClassA.NEEDSSKILNRS), ('home', DiagClassA.WELL)]:
        if key in lclRates:
            if tier == forceRelocateDiag:
                ratePairL.append((lclRates[key], ClassASetter(tier, forceRelocate=True)))
            else:
                ratePairL.append((lclRates[key], ClassASetter(tier)))
    changeTree = BayesTree.fromLinearCDF(ratePairL, tag='FATE')

    return changeTree


def biasTransfers(lclRates, infectiousAgent, scaleThisTier,
                  fromAbbrev, fromCategory, fromTier):
    """
    Split the lclRates table into two copies, one for healthy and the other for carriers.
    """
    enhancementScale = infectiousAgent.getRelativeProb(PthStatus.COLONIZED,
                                                       fromTier, scaleThisTier)
    estPthFrac = sum([infectiousAgent.getEstimatedPrevalence(pthStatus, fromAbbrev,
                                                             fromCategory,
                                                             fromTier)
                      for pthStatus in PthStatus.names
                      if pthStatus not in [PthStatus.CLEAR, PthStatus.RECOVERED]])
    scaleMe = tierKeyMap[scaleThisTier]
    beta = estPthFrac
    alpha = lclRates[scaleMe]
    k = enhancementScale * alpha
    x = (alpha - k * beta)/(1.0 - beta)
    if x < 0.0:
        logger.error('%s: Requested rescaling of rates for COLONIZED is too large; correcting!',
                     fromAbbrev)
        logger.error('beta = %s, alpha = %s, enhancementScale = %s, k = %s, x = %s',
                     beta, alpha, enhancementScale, k, x)
        x = 0.0
        k = alpha/beta

    if alpha == 0.0 or alpha == 1.0:
        newLclRates = lclRates.copy()
        pthRates = lclRates.copy()
    else:
        newLclRates = {}
        pthRates = {}
        for key, val in lclRates.items():
            if key == 'death':
                newLclRates[key] = pthRates[key] = val
            elif key == scaleMe:
                newLclRates[key] = (x/alpha) * val
                pthRates[key] = (k/alpha) * val
            elif key in knownTierKeys:  # but not ignored or scaled above
                pthRates[key] = ((1.0 - k)/(1.0-alpha)) * val
                newLclRates[key] = ((1.0 - x)/(1.0-alpha)) * val
            elif val == 0.0:
                pass
            else:
                raise RuntimeError('unexpected lclRates key %s' % key)
    return newLclRates, pthRates


class OverallHealthSetter(PatientStatusSetter):
    def __init__(self, newOverallHealth):
        self.newOverallHealth = newOverallHealth

    def set(self, patientStatus, timeNow):
        return patientStatus._replace(overall=self.newOverallHealth, startDateA=timeNow)

    def __str__(self):
        return ('PatientStatusSetter(overallHealth <- %s)'
                % PatientOverallHealth.names[self.newOverallHealth])


class Ward(BaseWard):
    def handlePatientDeparture(self, patientAgent, timeNow):
        """An opportunity for derived classes to customize the departure processing of patients"""
        super(Ward, self).handlePatientDeparture(patientAgent, timeNow)
#        patientAgent._status._replace(mostRecentHosp=self.fac.abbrev,
#                                      mostRecentHospReleaseTime=timeNow)
        for tP in self.fac.treatmentPolicies:
            tP.handlePatientDeparture(self, patientAgent, timeNow)
        self.miscCounters['departures'] += 1

class ICUWard(ForcedStateWard):
    def __init__(self, name, patch, nBeds):
        super(ICUWard, self).__init__(name, patch, CareTier.ICU, nBeds)
        
    def handlePatientArrival(self, patientAgent, timeNow):
        """
        Patient is arriving by triage, so force the patient's internal status and diagnosis
        to correspond to those of the ward.
        """
        self.forceState(patientAgent, CareTier.ICU, DiagClassA.VERYSICK, timeNow=timeNow)
        super(ICUWard, self).handlePatientArrival(patientAgent, timeNow)


class HospitalManager(FacilityManager):
    def allocateAvailableBed(self, requestedTier, strict=False):
        """
        If a HOSP tier bed is requested, maybe upgrade to ICU to simulate patients who
        are triaged to the ICU on arrival.
        """
        if strict:
            tier = requestedTier
        else:
            if (requestedTier == CareTier.HOSP
                and random() < _constants['fracTriageHOSPToICU']['value']):
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
        bedsPerICUWard = _constants['bedsPerICUWard']['value']

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

        for subCat, cat, key in [('icu', 'hospital', 'hospTransferToICURate'),
                                 ('vent', 'vsnf', 'vsnfTransferToVentRate'),
                                 ('skilnrs', 'vsnf', 'vsnfTransferToSkilNrsRate')]:
            self.lclRates[subCat] = _c[key]['value'] * self.lclRates[cat]
            self.lclRates[cat] -= self.lclRates[subCat]

        # VSNF-specific care tiers have been separated out, so any remaining vsnf fraction
        # is actually just nursing.
        self.lclRates['nursinghome'] += self.lclRates['vsnf']
        self.lclRates['vsnf'] = 0.0

        # pthRates is the equivalent of lclRates but for pathogen carriers.  It must
        # be initialized lazily because the pathogen has not yet been defined.
        self.pthRates = None

        self.icuDischargeViaDeathFrac = _constants['icuDischargeViaDeathFrac']['value']

        bedCountMultiplier = (_c['bedCountMultiplier']['value'] if 'bedCountMultiplier' in _c
                              else 1.0)
        if 'nBeds' in descr:
            allBeds = int(bedCountMultiplier * descr['nBeds']['value'])
        else:
            meanPop = descr['meanPop']['value']
            nBeds = int(bedCountMultiplier * 1.1 * meanPop)
        if 'nBedsICU' in descr:
            icuBeds = int(round(bedCountMultiplier * descr['nBedsICU']['value']))
        else:
            icuBeds = int(round(descr['fracAdultPatientDaysICU']['value'] * allBeds))
        nonICUBeds = max(allBeds - icuBeds, 0)

        assert descr['losModel']['pdf'] == 'lognorm(mu=$0,sigma=$1)', \
            'Hospital %(abbrev)s LOS PDF is not lognorm(mu=$0,sigma=$1)' % descr
        scaledLOSParms = self.scaleLOS(descr['losModel'],
                                       (descr['scaleLengthOfStay']['value'] if
                                        'scaleLengthOfStay' in descr else None))
        self.hospCachedCDF = CachedCDFGenerator(lognorm(scaledLOSParms[1],
                                                        scale=math.exp(scaledLOSParms[0])))
        self.hospTreeCache = {}
        if descr['meanLOSICU']['value'] == 0.0:
            self.icuCachedCDF = None
        else:
            # Construct a losModel struct from the mean and sigma
            mean = descr['meanLOSICU']['value']
            sigma = _constants['icuLOSLogNormSigma']['value']
            mu = math.log(mean) - (0.5 * sigma * sigma)
            icuLM = {'parms': [mu, sigma], 'pdf': 'lognorm(mu=$0,sigma=$1)',
                     'prov': 'from mean and sigma'}
            scaledICUParms = self.scaleLOS(icuLM,
                                           (descr['scaleLengthOfStay']['value'] if
                                            'scaleLengthOfStay' in descr else None))
            self.icuCachedCDF = CachedCDFGenerator(lognorm(scaledICUParms[1],
                                                           scale=math.exp(scaledICUParms[0])))
        self.icuTreeCache = {}

        for i, bedCt in enumerate(pickWardSizes(icuBeds, bedsPerICUWard)):
            self.addWard(ICUWard(('%s_%s_%s_%s_%d' %
                                  (category, patch.name, descr['abbrev'], 'ICU', i)),
                                 patch, bedCt))
        for i, bedCt in enumerate(pickWardSizes(nonICUBeds, bedsPerWard)):
            self.addWard(Ward(('%s_%s_%s_%s_%d' %
                               (category, patch.name, descr['abbrev'], 'HOSP', i)),
                              patch, CareTier.HOSP, bedCt))

    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        self.hospTreeCache = {}
        self.icuTreeCache = {}

    def getOrderedCandidateFacList(self, patientAgent, oldTier, newTier, modifierDct, timeNow):
        """Specialized to restrict transfers to being between our own HOSP and ICU"""
        queueClass = tierToQueueMap[newTier]
        if ((oldTier == CareTier.ICU and newTier == CareTier.HOSP)
                or (oldTier == CareTier.HOSP and newTier == CareTier.ICU)):
            qList = [rQ for rQ in self.reqQueues if isinstance(rQ, queueClass)]
#             print '%s: clause 1' % self.abbrev
            return [q.getGblAddr() for q in qList]
        else:
            facAddrList = super(Hospital, self).getOrderedCandidateFacList(patientAgent,
                                                                           oldTier, newTier,
                                                                           modifierDct, timeNow)
#             print '%s: clause 2: %s %s' % (self.abbrev, CareTier.names[newTier], facAddrList[:3])
        return facAddrList

    def getStatusChangeTree(self, patientAgent, modifierDct, startTime, timeNow):
        patientStatus = patientAgent.getStatus()
        ward = patientAgent.ward
        treatment = patientAgent.getTreatmentProtocol()
        if not self.pthRates:
            # Lazy initialization
            infectiousAgent = self.getWards(CareTier.HOSP)[0].iA
            self.lclRates, self.pthRates = biasTransfers(self.lclRates,
                                                         infectiousAgent,
                                                         CareTier.LTAC,
                                                         self.abbrev, self.category,
                                                         CareTier.HOSP)

        careTier = ward.tier
        if careTier == CareTier.HOSP:
            if patientAgent.getDiagnosis().diagClassA == DiagClassA.SICK:
                biasFlag = (patientStatus.pthStatus not in [PthStatus.CLEAR, PthStatus.RECOVERED])
                key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA,
                       biasFlag)
                if key in self.hospTreeCache:
                    return self.hospTreeCache[key]
                else:
                    if biasFlag:
                        changeTree = buildChangeTree(self.pthRates,
                                                     forceRelocateDiag=DiagClassA.SICK)
                    else:
                        changeTree = buildChangeTree(self.lclRates,
                                                     forceRelocateDiag=DiagClassA.SICK)

                    tree = BayesTree(changeTree,
                                     PatientStatusSetter(),
                                     self.hospCachedCDF.intervalProb, tag='LOS')
                    self.hospTreeCache[key] = tree
                    return tree
            else:
                # This patient doesn't belong in this ward
                logger.warning('fac %s patient: %s careTier %s overall %s with status %s '
                               'startTime: %s: this patient should be gone by now',
                               self.name, patientAgent.name, CareTier.names[careTier],
                               PatientOverallHealth.names[patientStatus.overall],
                               DiagClassA.names[patientStatus.diagClassA], startTime)
                return BayesTree(PatientStatusSetter())

        elif careTier == CareTier.ICU:
            if patientAgent.getDiagnosis().diagClassA == DiagClassA.VERYSICK:
                key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
                if key in self.icuTreeCache:
                    return self.icuTreeCache[key]
                else:
                    tree = BayesTree(BayesTree.fromLinearCDF([(self.icuDischargeViaDeathFrac,
                                                               ClassASetter(DiagClassA.DEATH)),
                                                              ((1.0 - self.icuDischargeViaDeathFrac),
                                                               ClassASetter(DiagClassA.SICK))],
                                                             tag='FATE'),
                                     PatientStatusSetter(),
                                     self.icuCachedCDF.intervalProb, tag='LOS')
                    self.icuTreeCache[key] = tree
                    return tree
            else:
                # This patient doesn't belong in this ward
                logger.warning('fac %s patient: %s careTier %s with status %s startTime: %s: '
                               'this patient should be gone by now'
                               % (self.name, patientAgent.name, CareTier.names[careTier],
                                  DiagClassA.names[patientStatus.diagClassA], startTime))
                return BayesTree(PatientStatusSetter())
        else:
            raise RuntimeError('Hospitals do not provide care tier %s' % careTier)

    def getInitialOverallHealth(self, ward, timeNow):  # @UnusedVariable
        # Hospital patients have by definition been in a hospital in the last year.
        return PatientOverallHealth.UNHEALTHY


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
        ward.handlePatientArrival(a, None)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              None)
        agentList.append(a)
    for i in xrange(int(round(meanHospPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.HOSP, strict=True)
        assert ward is not None, 'Ran out of HOSP beds populating %(abbrev)s!' % descr
        a = PatientAgent('PatientAgent_HOSP_%s_%d' % (ward._name, i), patch, ward)
        ward.lock(a)
        ward.handlePatientArrival(a, None)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              None)
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
