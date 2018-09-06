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
from scipy.stats import lognorm, expon, weibull_min
import logging

import pyrheabase
import pyrheautils
import schemautils
from quilt.peopleplaces import FutureMsg
from stats import CachedCDFGenerator, lognormplusexp, BayesTree
from facilitybase import DiagClassA, CareTier, TreatmentProtocol, NURSINGQueue
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent
from facilitybase import PatientStatusSetter, buildTimeTupleList, FacilityManager
from facilitybase import findQueueForTier
from hospital import estimateWork as hospitalEstimateWork
from hospital import ClassASetter, OverallHealthSetter
from hospital import buildChangeTree

logger = logging.getLogger(__name__)

category = 'NURSINGHOME'
_schema = 'nursinghomefacts_ChicagoLand_schema.yaml'
_constants_values = '$(MODELDIR)/constants/nursinghome_constants.yaml'
_constants_schema = 'nursinghome_ChicagoLand_constants_schema.yaml'
_constants = None
_validator = None


class NursingWard(Ward):
    def handlePatientArrival(self, patientAgent, timeNow):
        super(NursingWard, self).handlePatientArrival(patientAgent, timeNow)
        # Set this to be the patient's "home" so held rooms can be found later
        patientAgent.setStatus(homeAddr=findQueueForTier(CareTier.NURSING,
                                                         self.fac.reqQueues).getGblAddr())


class CancelHoldMsg(FutureMsg):
    idCounter = 0

    @classmethod
    def nextId(cls):
        rslt = cls.idCounter
        cls.idCounter += 1
        return rslt

    def __init__(self, baseName, patch, payload, destAddr, arrivalTime, debug=False):
        super(CancelHoldMsg, self).__init__(baseName + '_cancel_hold_%d' % self.nextId(),
                                            patch, payload, destAddr, arrivalTime,
                                            debug=debug)


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
        if 'bedCountMultiplier' in _c:
            nBeds = int(nBeds * _c['bedCountMultiplier']['value'])
        losModel = descr['losModel']
        lMP = losModel['parms']
        if losModel['pdf'] == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)':
            self.initialResidentFrac = (1.0 - losModel['parms'][0])
            self.rehabCachedCDF = CachedCDFGenerator(lognorm(lMP[2],
                                                             scale=math.exp(lMP[1])))
            self.frailCachedCDF = CachedCDFGenerator(expon(scale=1.0/lMP[3]))
        elif losModel['pdf'] == '$0*weibull(k=$1, lmda=$2)+(1-$0)*weibull(k=$3, lmda=$4)':
            self.initialResidentFrac = (1.0 - losModel['parms'][0])
            self.rehabCachedCDF = CachedCDFGenerator(weibull_min(lMP[1], scale=lMP[2]))
            self.frailCachedCDF = CachedCDFGenerator(weibull_min(lMP[3], scale=lMP[4]))
        else:
            raise RuntimeError("Unexpected losModel form %s for %s!" % (losModel['pdf'],
                                                                        descr['abbrev']))

        self.initialUnhealthyFrac = _c['initialUnhealthyFrac']['value'] # frac of non-residents
        self.initialNonResidentFrailFrac = _c['initialNonResidentFrailFrac']['value']

        initialFrac = (self.initialResidentFrac
                       + (1.0-self.initialResidentFrac)*self.initialNonResidentFrailFrac)
        frailBeds = int(math.ceil(initialFrac * nBeds))
        self.bedAllocDict = {'frail_beds': frailBeds, 'non_frail_beds': (nBeds - frailBeds),
                             'frail': 0, 'non_frail': 0,
                             'frail_held': 0, 'non_frail_held': 0}
        self.bedHoldTime = _c['heldBedDurationDays']['value']

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

        self.rehabTreeCache = {}
        self.frailTreeCache = {}
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
        self.frailTreeCache = {}

    def checkBedAllocDict(self, healthKey, heldKey, bedsKey):
        """Check that the bed allocation table has been initialized and is sensible."""
        bAD = self.bedAllocDict
        assert bAD[heldKey] >= 0, '%s %s held bed count is less than 0' % (self.name, healthKey)
        assert bAD[healthKey] >= 0, '%s %s occupancy is less than 0' % (self.name, healthKey)
        assert bAD[bedsKey] >= 0, '%s %s bed count is less than 0' % (self.name, healthKey)
        assert bAD[heldKey] + bAD[healthKey] <= bAD[bedsKey], ('%s %s too many beds committed'
                                                               % (self.name, healthKey))

    def handleBedRequestResponse(self, ward, payload, timeNow):  # @UnusedVariable
        """
        This routine is called in the time slice of the Facility Manager when the manager
        responds to a request for a bed.  If the request was denied, 'ward' will be None.
        The return value is the tuple (newPayload, ward).  Returning None for the ward
        value of the tuple will cause the request to be denied.

        If ward is not None, this facility is in the process of accepting the bed request
        and the associated patient will be arriving later in the day.  Thus we cache
        the transfer info which may be associated with the patient.
        """
        outgoingPayload, ward = super(NursingHome, self).handleBedRequestResponse(ward,
                                                                                  payload,
                                                                                  timeNow)
        if ward:
            # We may yet decline if it would violate our overall health ratio
            transferInfoDict = payload[3]
            pId = transferInfoDict['patientId']
            pOH = transferInfoDict['overallHealth']
            if self.patientRecordExists(pId):
                pRec = self.getPatientRecord(pId)
                if 'bedHeld' in pRec.noteD and pRec.noteD['bedHeld']:
                    bedHeld = True
                else:
                    bedHeld = False
            else:
                bedHeld = False
            if pOH == PatientOverallHealth.FRAIL:
                healthKey, heldKey, bedsKey = 'frail', 'frail_held', 'frail_beds'
            else:
                healthKey, heldKey, bedsKey = 'non_frail', 'non_frail_held', 'non_frail_beds'
            self.checkBedAllocDict(healthKey, heldKey, bedsKey)
            if bedHeld:
                # print '%s: reentered held bed! %s %s' % (self.name, healthKey, self.bedAllocDict)
                pass
            else:
                if (self.bedAllocDict[healthKey] < (self.bedAllocDict[bedsKey]
                                                    - self.bedAllocDict[heldKey])):
                    self.bedAllocDict[heldKey] += 1  # commit the bed
                else:
                    ward = None  # decline the request
                    # print '%s: request declined %s %s' % (self.name, healthKey, self.bedAllocDict)
            if ward is None:
                del self.arrivingPatientTransferInfoDict[pId]  # entry added by superclass
        else:
            pass
            # healthKey = ('frail' if payload[3]['overallHealth'] == PatientOverallHealth.FRAIL
            #              else 'non_frail')
            # print '%s %s ward is None %s' % (self.name, healthKey, self.bedAllocDict)
        return outgoingPayload, ward

    def getMsgPayload(self, msgType, patientAgent):
        payload = super(NursingHome, self).getMsgPayload(msgType, patientAgent)
        if issubclass(msgType, pyrheabase.DepartureMsg):
            # stash the overall health and patient Id
            payload = (patientAgent.getStatus().overall, patientAgent.id,
                       patientAgent.tier, payload)  # tier is that of destination at this point
        elif issubclass(msgType, pyrheabase.ArrivalMsg):
            # stash overall health
            payload = (patientAgent.getStatus().overall, patientAgent.id, payload)
        return payload

    def handleIncomingMsg(self, msgType, payload, timeNow):
        if issubclass(msgType, pyrheabase.DepartureMsg):
            pOH, pId, destTier, payload = payload  # strip out overall health and patient Id.
            if pOH == PatientOverallHealth.FRAIL:
                healthKey, heldKey, bedsKey = 'frail', 'frail_held', 'frail_beds'
            else:
                healthKey, heldKey, bedsKey = 'non_frail', 'non_frail_held', 'non_frail_beds'
            self.checkBedAllocDict(healthKey, heldKey, bedsKey)
            self.bedAllocDict[healthKey] -= 1
            if destTier not in [CareTier.HOME, CareTier.NURSING] and destTier is not None:
                self.bedAllocDict[heldKey] += 1
                with self.getPatientRecord(pId, timeNow) as pRec:
                    pRec.noteD['bedHeld'] = True
                cancelHoldMsg = CancelHoldMsg(self.name, self.manager.patch,
                                              (pOH, pId, timeNow),
                                              self.reqQueues[0].getGblAddr(),
                                              timeNow + self.bedHoldTime)
                self.manager.patch.launch(cancelHoldMsg, timeNow)
            return super(NursingHome, self).handleIncomingMsg(msgType, payload, timeNow)
        elif issubclass(msgType, pyrheabase.ArrivalMsg):
            pOH, pId, payload = payload # strip out overall health and patient id
            if pOH == PatientOverallHealth.FRAIL:
                (healthKey, heldKey, bedsKey,
                 otherHealthKey, otherBedsKey) = ('frail', 'frail_held', 'frail_beds',
                                                  'non_frail', 'non_frail_beds')
            else:
                (healthKey, heldKey, bedsKey,
                 otherHealthKey, otherBedsKey) = ('non_frail', 'non_frail_held', 'non_frail_beds',
                                                  'frail', 'frail_beds')
            bAD = self.bedAllocDict
            bAD[healthKey] += 1
            if bAD[heldKey] > 0:
                bAD[heldKey] -= 1
            else:
                if timeNow > 0:
                    logger.error('%s unexpected arrival of %s %s', self.name, healthKey, pId)
                    self.checkBedAllocDict(healthKey, heldKey, bedsKey)
            if bAD[healthKey] > bAD[bedsKey]:
                delta = bAD[healthKey] - bAD[bedsKey]
                assert bAD[otherHealthKey] <= bAD[otherBedsKey] + delta, ('%s overflow arrival'
                                                                          % self.name)
                bAD[bedsKey] += delta
                bAD[otherBedsKey] -= delta
            with self.getPatientRecord(pId, timeNow) as pRec:
                pRec.noteD['bedHeld'] = False
            return super(NursingHome, self).handleIncomingMsg(msgType, payload, timeNow)
        elif issubclass(msgType, CancelHoldMsg):
            pOH, pId, launchTime = payload
            if pOH == PatientOverallHealth.FRAIL:
                (healthKey, heldKey, bedsKey,
                 otherHealthKey, otherBedsKey) = ('frail', 'frail_held', 'frail_beds',
                                                  'non_frail', 'non_frail_beds')
            else:
                (healthKey, heldKey, bedsKey,
                 otherHealthKey, otherBedsKey) = ('non_frail', 'non_frail_held', 'non_frail_beds',
                                                  'frail', 'frail_beds')
            with self.getPatientRecord(pId, timeNow) as pRec:
                if pRec.departureDate == launchTime and pRec.noteD['bedHeld']:
                    pRec.noteD['bedHeld'] = False
                    bAD = self.bedAllocDict
                    bAD[heldKey] -= 1
                    assert bAD[heldKey] >= 0, '%s cannot find bed hold to cancel'
            return timeNow

    def getStatusChangeTree(self, patientAgent, modifierDct, startTime, timeNow):
        patientStatus = patientAgent.getStatus()
        patientDiagnosis = patientAgent.getDiagnosis()
        ward = patientAgent.ward
        treatment = patientAgent.getTreatmentProtocol()
        careTier = ward.tier
        assert careTier == CareTier.NURSING, \
            "Nursing homes only offer CareTier 'NURSING'; found %s" % careTier
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        _c = _constants
        if patientDiagnosis.overall == PatientOverallHealth.FRAIL:
            if key in self.frailTreeCache:
                return self.frailTreeCache[key]
            else:
                innerTree = ClassASetter(DiagClassA.WELL)  # We declare them done with any rehab
                changeTree = buildChangeTree(self.lclRates, forceRelocateDiag=DiagClassA.NEEDSREHAB)
                tree = BayesTree(changeTree, innerTree,
                                 self.frailCachedCDF.intervalProb, tag='LOS')
                self.frailTreeCache[key] = tree
                return tree

        else:  # overall health is not FRAIL, hence HEALTHY or UNHEALTHY
            if treatment.rehab:
                if key in self.rehabTreeCache:
                    return self.rehabTreeCache[key]
                else:
                    changeTree = buildChangeTree(self.lclRates,
                                                 forceRelocateDiag=DiagClassA.NEEDSREHAB)
                    tree = BayesTree(changeTree,
                                     PatientStatusSetter(),
                                     self.rehabCachedCDF.intervalProb, tag='LOS')
                    self.rehabTreeCache[key] = tree
                    return tree
            else:  # not frail and not rehab
                if patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
                    raise RuntimeError('%s has a REHAB diagnosis but no treatment?')
            logger.warning('fac %s patient: %s careTier %s overall %s with status %s startTime: %s: '
                           'this patient should be gone by now'
                           % (self.name, patientAgent.name, CareTier.names[careTier],
                              PatientOverallHealth.names[patientStatus.overall],
                              DiagClassA.names[patientStatus.diagClassA], startTime))
            return BayesTree(PatientStatusSetter())

    def getInitialOverallHealth(self, ward, timeNow):  # @UnusedVariable
        if random() <= self.initialResidentFrac:
            return PatientOverallHealth.FRAIL
        else:
            rn = random()
            if rn <= self.initialUnhealthyFrac:
                return PatientOverallHealth.UNHEALTHY
            elif rn <= (self.initialUnhealthyFrac + self.initialNonResidentFrailFrac):
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
    # The following is approximate, but adequate...
    agentList = []
    logger.debug('%s before _populate %s: %s', fac.name, int(round(meanPop)), fac.bedAllocDict)
    for i in xrange(int(round(meanPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.NURSING)
        assert ward is not None, 'Ran out of beds populating %(abbrev)s!' % descr
        a = PatientAgent('PatientAgent_NURSING_%s_%d' % (ward._name, i), patch, ward)
        a.setStatus(homeAddr=findQueueForTier(CareTier.NURSING, fac.reqQueues).getGblAddr())  # They live here
        if a.getStatus().overall != PatientOverallHealth.FRAIL:
            a.setTreatment(rehab=True)  # They must be here for rehab
        ward.lock(a)
        ward.handlePatientArrival(a, None)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              None)
        agentList.append(a)
    logger.debug('%s after _populate %s: %s', fac.name, int(round(meanPop)), fac.bedAllocDict)
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
