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

from random import randint, choice
import logging
import math
from collections import defaultdict
import cPickle as pickle
import numpy as np

from phacsl.utils.notes.statval import HistoVal
from phacsl.utils.collections.phacollections import SingletonMetaClass
import pyrheabase
from typebase import CareTier, PatientOverallHealth, DiagClassA
from typebase import TreatmentProtocol, TREATMENT_DEFAULT  # @UnusedImport
from typebase import PatientStatus, PatientDiagnosis  # @UnusedImport
from stats import BayesTree
from pathogenbase import PthStatus
from registry import Registry  # @UnusedImport
from labwork import LabWork, LabWorkMsg

from policybase import TransferDestinationPolicy, TreatmentPolicy, DiagnosticPolicy

LOGGER = logging.getLogger(__name__)
infectionLogger = logging.getLogger("infectionTracking")

HackBedMultiplier = 1

class PatientStatusSetter(object):
    def __init__(self):
        pass

    def set(self, patientStatus, timeNow):  # @UnusedVariable
        """This base class just returns a copy"""
        return PatientStatus._make(patientStatus)

    def __str__(self):
        return 'PatientStatusSetter()'


class ClassASetter(PatientStatusSetter):
    def __init__(self, newClassA, forceRelocate=False):
        super(ClassASetter, self).__init__()
        self.newClassA = newClassA
        self.forceRelocate = forceRelocate

    def set(self, patientStatus, timeNow):
        if self.forceRelocate:
            if patientStatus.diagClassA == self.newClassA:
                return (patientStatus._replace(relocateFlag=True))
            else:
                return (patientStatus._replace(diagClassA=self.newClassA, startDateA=timeNow)
                        ._replace(relocateFlag=True))
        else:
            if patientStatus.diagClassA == self.newClassA:
                return (patientStatus)
            else:
                return patientStatus._replace(diagClassA=self.newClassA, startDateA=timeNow)

    def __str__(self):
        if self.forceRelocate:
            return ('PatientStatusSetter(classA <- %s, forceRelocate=True)'
                    % DiagClassA.names[self.newClassA])
        else:
            return 'PatientStatusSetter(classA <- %s)' % DiagClassA.names[self.newClassA]

    def __repr__(self):
        return ('PatientStatusSetter(classA <- %s, forceRelocate=%s)'
                % (DiagClassA.names[self.newClassA], self.forceRelocate))


class PthStatusSetter(PatientStatusSetter):
    def __init__(self, newPthStatus):
        self.newPthStatus = newPthStatus

    def set(self, patientStatus, timeNow):
        return patientStatus._replace(pthStatus=self.newPthStatus, startDatePth=timeNow)

    def __str__(self):
        return 'PatientStatusSetter(pthStatus <- %s)' % PthStatus.names[self.newPthStatus]


class PatientCumStats(object):
    def __init__(self):
        nEntries = len(PatientOverallHealth.names)
        self.valV = np.zeros([2 * nEntries], dtype=np.int32)  # upper half is for total ever

    def incrPatient(self, patientAgent):
        nEntries = len(PatientOverallHealth.names)
        pOH = patientAgent.getStatus().overall
        self.valV[pOH] += 1
        self.valV[pOH + nEntries] += 1

    def decrPatient(self, patientAgent):
        pOH = patientAgent.getStatus().overall
        self.valV[pOH] -= 1

    def popTotal(self):
        nEntries = len(PatientOverallHealth.names)
        return np.sum(self.valV[:nEntries])

    def popEver(self):
        nEntries = len(PatientOverallHealth.names)
        return np.sum(self.valV[nEntries:])

    def popByOH(self, patientOverallHealth):
        return self.valV[patientOverallHealth]


class Ward(pyrheabase.Ward):
    def __init__(self, name, patch, tier, nBeds):
        pyrheabase.Ward.__init__(self, name, patch, tier, int(HackBedMultiplier*nBeds))
        self.checkInterval = 1  # check health daily
        self.iA = None  # infectious agent
        self.miscCounters = defaultdict(lambda: 0)
        self.cumStats = PatientCumStats()

        try:
            self.wardNum = int(name.split('_')[-1])
        except ValueError:
            self.wardNum = 0

    def getPatientList(self):
        return self.getLiveLockedAgents()

    def setInfectiousAgent(self, iA):
        self.iA = iA

    def initializePatientPthState(self):
        for p in self.getPatientList():
            self.iA.initializePatientState(p)

    def initializePatientTreatment(self):
        for p in self.getPatientList():
            for tP in self.fac.treatmentPolicies:
                tP.initializePatientTreatment(self, p)

    def handlePatientArrival(self, patientAgent, timeNow):
        """An opportunity for derived classes to customize the arrival processing of patients"""
        patientAgent.setStatus(justArrived=True)
        self.cumStats.incrPatient(patientAgent)
        self.miscCounters['arrivals'] += 1
        if patientAgent.getStatus().pthStatus == PthStatus.COLONIZED:
            self.miscCounters['creArrivals'] += 1
        if patientAgent.id in self.fac.arrivingPatientTransferInfoDict.keys():
            transferInfo = self.fac.arrivingPatientTransferInfoDict[patientAgent.id]
            del self.fac.arrivingPatientTransferInfoDict[patientAgent.id]
        else:
            transferInfo = self.fac.getBedRequestPayload(patientAgent, self.tier)[-1]

        self.fac.diagnosticPolicy.handlePatientArrival(self, patientAgent, transferInfo,
                                                       timeNow)
        patientAgent.addHistoryEntry(self, timeNow)
        for tP in self.fac.treatmentPolicies:
            tP.handlePatientArrival(self, patientAgent, transferInfo, timeNow)

    def handlePatientDeparture(self, patientAgent, timeNow):
        """An opportunity for derived classes to customize the departure processing of patients"""
        for tP in self.fac.treatmentPolicies:
            tP.handlePatientDeparture(self, patientAgent, timeNow)
        self.fac.diagnosticPolicy.handlePatientDeparture(self, patientAgent, timeNow)
        self.miscCounters['departures'] += 1
        self.cumStats.decrPatient(patientAgent)


class ForcedStateWard(Ward):
    def forceState(self, patientAgent, careTier, diagClassA, timeNow=None):

        patientAgent.setStatus(diagClassA=diagClassA)
        patientAgent._diagnosis = self.fac.diagnose(patientAgent.ward,
                                                    patientAgent.id,
                                                    patientAgent.getStatus(),
                                                    patientAgent.getDiagnosis(),
                                                    timeNow=timeNow)
        newTier, patientAgent._treatment = self.fac.prescribe(self, patientAgent.id,
                                                              patientAgent.getDiagnosis(),
                                                              patientAgent.getTreatmentProtocol(),
                                                              {},
                                                              timeNow=timeNow)[0:2]
        assert newTier == self.tier, ('%s %s %s tried to force state of %s to match but failed'
                                      % (self._name, CareTier.names[careTier],
                                         DiagClassA.names[diagClassA],
                                         patientAgent.name))


class FacilityManager(pyrheabase.FacilityManager):
    pass


class FacRequestQueue(pyrheabase.FacRequestQueue):
    def getInfo(self):
        return (super(FacRequestQueue, self).getInfo(),
                self._lockingAgent.fac.abbrev, self._lockingAgent.fac.coords)


class BirthQueue(FacRequestQueue):
    pass


class BirthMsg(pyrheabase.SimpleMsg):
    pass


class ICUQueue(FacRequestQueue):
    pass


class HOSPQueue(FacRequestQueue):
    pass


class LTACQueue(FacRequestQueue):
    pass


class NURSINGQueue(FacRequestQueue):
    pass


class HOMEQueue(FacRequestQueue):
    pass

class VENTQueue(FacRequestQueue):
    pass

class SKILNRSQueue(FacRequestQueue):
    pass


tierToQueueMap = {CareTier.HOME: HOMEQueue,
                  CareTier.NURSING: NURSINGQueue,
                  CareTier.LTAC: LTACQueue,
                  CareTier.HOSP: HOSPQueue,
                  CareTier.ICU: ICUQueue,
                  CareTier.VENT: VENTQueue,
                  CareTier.SKILNRS: SKILNRSQueue}


def findQueueForTier(careTier, queueList):
    """
    Given a list of queues of various types (for example Facility.reqQueues), find
    and return the first which matches the given CareTier.
    """
    for queue in queueList:
        if isinstance(queue, tierToQueueMap[careTier]):
            return queue
    return None


class PatientRecord(object):
    _boolProps = ['isFrail',
                  'carriesPth',  # Carries the pathogen being simulated
                  'carriesOther', # Carries some other contagious pathogen
                  ]

    def __init__(self, patientId, arrivalDate, #owningFacility,
                 isFrail, carriesPth=False, carriesOther=False):
        self.patientId = patientId
        self.arrivalDate = arrivalDate
        #self.fac = owningFacility
        self.departureDate = None
        self.prevVisits = 0
        for propN in PatientRecord._boolProps:
            setattr(self, propN, locals()[propN])
        self.noteD = {}

    def __enter__(self):
        assert hasattr(self, '_owningFac'), ('PatientRecord can only be a context if'
                                             ' created via Facility.getPatientRecord')
        return self

    def __exit__(self, tp, val, tb):  # @UnusedVariable
        self._owningFac.mergePatientRecord(self.patientId, self, None)
    
    def forgetPathogenInfo(self):
        """
        Scenarios sometimes require that patient information be 'lost', but we don't want to
        lose it completely for bookkeeping reasons.  This method loses only the pathogen-
        related parts of the record.
        """
        self.carriesPth = False
        self.carriesOther = False

    @property
    def isContagious(self):
        return (self.carriesPth or self.carriesOther)

    def merge(self, otherRec):
        # print 'merge %s and %s %s' % (str(self), otherRec, dir(otherRec))
        assert self.patientId == otherRec.patientId, 'Cannot merge records for different patients'
        if otherRec.arrivalDate > self.arrivalDate:     # newer visit
            for propN in PatientRecord._boolProps + ['arrivalDate', 'departureDate',
                                                     'prevVisits']:
                setattr(self, propN, getattr(otherRec, propN)) # new rec overrides
            self.noteD.update(otherRec.noteD)  # new rec overrides
        elif otherRec.arrivalDate < self.arrivalDate: # earlier visit
            assert (otherRec.departureDate is not None
                    and (self.departureDate is None
                         or otherRec.departureDate <= self.arrivalDate)), 'record-keeping inconsistency'
            # Bool fields are more recent and thus stay the same
            # merge notes, keeping more recent
            newD = otherRec.noteD
            newD.update(self.noteD)
            self.noteD = newD

        else:  # same visit
            if self.departureDate is None:
                self.departureDate = otherRec.departureDate
            else:
                if (otherRec.departureDate is not None
                        and otherRec.departureDate > self.departureDate):
                    self.departureDate = otherRec.departureDate
            self.prevVisits = max(self.prevVisits, otherRec.prevVisits)
            for propN in PatientRecord._boolProps:
                setattr(self, propN, (getattr(self, propN) or getattr(otherRec, propN)))
            self.noteD.update(otherRec.noteD)  # new rec overrides

        for key, val in self.noteD.items():
            if val is None:
                del self.noteD[key]

    def __str__(self):
        return '<patient %s, %s -> %s, %s, %s>' % (self.patientId,
                                                   self.arrivalDate,
                                                   self.departureDate,
                                                   'Frail' if self.isFrail else 'not frail',
                                                   ('Contagious' if self.isContagious
                                                    else 'NonContagious'))

class PatientStats(object):
    def __init__(self, fac):
        self.fac = fac

    @property
    def currentOccupancy(self):
        return sum([ward.cumStats.popTotal() for ward in self.fac.getWards()])

    @property
    def totalOccupancy(self):
        return sum([ward.cumStats.popEver() for ward in self.fac.getWards()])

    def __str__(self):
        dct = defaultdict(int)
        for ward in self.fac.getWards():
            dct[ward.tier] += ward.cumStats.popTotal()
        bedStr = ', '.join(['%s: %d' % (CareTier.names[tier], ct)
                            for tier, ct in dct.items()])
        return '<occupied beds %s>' % bedStr


class MissingPatientRecordError(RuntimeError):
    pass


class Facility(pyrheabase.Facility):
    def __init__(self, name, descr, patch, reqQueueClasses=None, policyClasses=None,
                 managerClass=FacilityManager, categoryNameMapper=None):
        """
        If provided, categoryNameMapper should be a function with the signature:

           implCat = categoryNameMapper(descrCat)

        where implCat is a facility implementation category name (e.g. NURSINGHOME) and
        descrCat is a facility description category name (e.g. SNF).
        """
        pyrheabase.Facility.__init__(self, name, patch, managerClass=managerClass,
                                     reqQueueClasses=reqQueueClasses)
        if categoryNameMapper is None:
            self.categoryNameMapper = (lambda(descrCat): descrCat)
        else:
            self.categoryNameMapper = categoryNameMapper
#         self.category = self.categoryNameMapper(descr['category'])
        self.category = descr['category']
        if categoryNameMapper:
            self.implCategory = categoryNameMapper(self.category)
        else:
            self.implCategory = self.category
        self.abbrev = descr['abbrev']
        if 'longitude' in descr and 'latitude' in descr:
            self.coords = (descr['longitude'], descr['latitude'])
        else:
            self.coords = (None, None)
        self.noteHolder = None
        self.idCounter = 0
        self.patientDataDict = {}
        self.patientStats = PatientStats(self)
        transferDestinationPolicyClass = TransferDestinationPolicy
        treatmentPolicyClasses = []
        diagnosticPolicyClass = DiagnosticPolicy
        #print "PolicyClasses for {0} is {1}".format(name, policyClasses)
        if policyClasses is not None:
            for pC in policyClasses:
                if issubclass(pC, TransferDestinationPolicy):
                    transferDestinationPolicyClass = pC
                if issubclass(pC, TreatmentPolicy):
                    treatmentPolicyClasses.append(pC)
                if issubclass(pC, DiagnosticPolicy):
                    diagnosticPolicyClass = pC
        if not treatmentPolicyClasses:
            treatmentPolicyClasses = [TreatmentPolicy]
        self.transferDestinationPolicy = transferDestinationPolicyClass(patch,
                                                                        self.categoryNameMapper)
        self.treatmentPolicies = [treatmentPolicyClass(patch, self.categoryNameMapper)
                                  for treatmentPolicyClass in treatmentPolicyClasses]
        self.diagnosticPolicy = diagnosticPolicyClass(self, patch, self.categoryNameMapper)
        self.arrivingPatientTransferInfoDict = {}

    def finalizeBuild(self, facDescription):
        """
        This provides an opportunity for the facility to do any bookkeeping necessary after all
        wards and agents have been fully initialized.
        """
        pass

    def __str__(self):
        return '<%s>' % self.name

    def setNoteHolder(self, noteHolder):
        self.noteHolder = noteHolder
        self.noteHolder.addNote({'name': self.name})

    def getNoteHolder(self):
        return self.noteHolder

    def getPatientRecord(self, patientId, timeNow=None):
        if patientId in self.patientDataDict:
            pR = pickle.loads(self.patientDataDict[patientId])
            pR._owningFac = self
            return pR
        else:
            # Create a new blank record
            if timeNow is None:
                raise MissingPatientRecordError('Record not found and cannot create a new one')
            pR = PatientRecord(patientId, timeNow, isFrail=False)
            self.patientDataDict[patientId] = pickle.dumps(pR, 2)  # keep a copy
            pR._owningFac = self
            return pR

    def forgetPatientRecord(self, patientId):
        """
        The facility forgets it ever saw this patient.  Used to implement record-keeping errors.
        Completely removing the record causes too many bookkeeping problems, so we just wipe the
        pathogen-related parts.
        """
        if patientId in self.patientDataDict:
            with self.getPatientRecord(patientId) as pRec:
                pRec.forgetPathogenInfo()

    def getPatientRecords(self):
        """In case someone wants to exhaustively search patient records"""
        for pStr in self.patientDataDict.values():
            yield pickle.loads(pStr)

    def mergePatientRecord(self, patientId, newPatientRec, timeNow):
        patientRec = self.getPatientRecord(patientId, timeNow)
        patientRec.merge(newPatientRec)
        delattr(patientRec, '_owningFac')
        self.patientDataDict[patientId] = pickle.dumps(patientRec, 2)

    def patientRecordExists(self, patientId):
        return patientId in self.patientDataDict

    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        pass

    def getOrderedCandidateFacList(self, patientAgent, oldTier, newTier, modifierDct, timeNow):
        oCFL = self.transferDestinationPolicy.getOrderedCandidateFacList(self,
                                                                         patientAgent,
                                                                         oldTier, newTier,
                                                                         modifierDct,
                                                                         timeNow)
        return oCFL

    def getMsgPayload(self, msgType, patientAgent):
        if issubclass(msgType, (pyrheabase.ArrivalMsg, pyrheabase.DepartureMsg)):
            innerPayload = super(Facility, self).getMsgPayload(msgType, patientAgent)
            return ((patientAgent.id, patientAgent.ward.tier,
                     patientAgent.getStatus().overall == PatientOverallHealth.FRAIL,
                     patientAgent.name),
                    innerPayload)
        elif issubclass(msgType, BirthMsg):
            return patientAgent.getStatus().overall
        else:
            raise RuntimeError('%s: payload request for unknown message type %s'
                               % (self.name, msgType.__name__))

    def handleIncomingMsg(self, msgType, payload, timeNow):
        if issubclass(msgType, pyrheabase.ArrivalMsg):
            myPayload, innerPayload = payload
            timeNow = super(Facility, self).handleIncomingMsg(msgType, innerPayload, timeNow)
            patientId, tier, isFrail, patientName = myPayload  # @UnusedVariable
            patientRec = self.getPatientRecord(patientId,
                                               (timeNow if timeNow is not None else 0))
            patientRec.isFrail = isFrail
            if timeNow is not None:
                if patientRec.arrivalDate != timeNow:
                    LOGGER.debug('Patient %s has returned to %s', patientId, self.name)
                    patientRec.arrivalDate = timeNow
                    patientRec.prevVisits += 1
                if timeNow != 0:  # Exclude initial populations
                    nh = self.getNoteHolder()
                    if nh:
                        nh.addNote({(CareTier.names[tier] + '_arrivals'): 1})
            patientRec.departureDate = None
            self.mergePatientRecord(patientId, patientRec,
                                    (timeNow if timeNow is not None else 0))
        elif issubclass(msgType, pyrheabase.DepartureMsg):
            if timeNow is None:
                raise RuntimeError('Only arrival messages should happen before execution starts')
            myPayload, innerPayload = payload
            timeNow = super(Facility, self).handleIncomingMsg(msgType, innerPayload, timeNow)
            patientId, tier, isFrail, patientName = myPayload # @UnusedVariable
            if self.patientRecordExists(patientId):
                with self.getPatientRecord(patientId, timeNow=timeNow) as patientRec:
                    patientRec.departureDate = timeNow
                    lengthOfStay = timeNow - patientRec.arrivalDate
                    losKey = (CareTier.names[tier] + '_LOS')
                    nh = self.getNoteHolder()
                    if losKey not in nh:
                        nh.addNote({losKey: HistoVal([])})
                    nh.addNote({(CareTier.names[tier] + '_departures'): 1, losKey: lengthOfStay})
            else:
                # this can happen if the scenario calls for the patientRec to be 'lost'
                LOGGER.error('%s has no record of patient %s', self.name, patientId)
                nh = self.getNoteHolder()
                nh.addNote({(CareTier.names[tier] + '_departures'): 1})
        elif issubclass(msgType, BirthMsg):
            if timeNow is None:
                raise RuntimeError('Only arrival messages should happen before execution starts')
            ward = self.manager.allocateAvailableBed(CareTier.HOME)
            assert ward is not None, 'Ran out of beds with birth in %s!' % self.name
            a = PatientAgent('PatientAgent_%s_%d_birth' % (ward._name, self.idCounter),
                             self.manager.patch, ward)
            self.idCounter += 1
            a.setStatus(homeAddr=findQueueForTier(ward.tier, self.reqQueues).getGblAddr())
            a.setStatus(overall=payload)
            ward.handlePatientArrival(a, timeNow)
            self.handleIncomingMsg(pyrheabase.ArrivalMsg,
                                   self.getMsgPayload(pyrheabase.ArrivalMsg, a),
                                   timeNow)
            self.manager.patch.launch(a, timeNow)
            if timeNow != 0:  # exclude initial populations
                nh = self.getNoteHolder()
                nh.addNote({'births': 1})
        elif issubclass(msgType, LabWorkMsg):
            LabWork.handleLabMsg(self, msgType, payload, timeNow)
        else:
            raise RuntimeError('%s: got unknown message type %s' % (self.name,
                                                                    msgType .__name__))
        return timeNow

    def getBedRequestPayload(self, patientAgent, desiredTier):
        """
        The return value defines the contents of the BedRequest payload and must
        be parsed in the other BedRequest methods of the facility.

        The format used here is (number of bounces, tier, originating facility, and
        transferInfo dictionary).
        """
        transferInfoDict = {'patientId': patientAgent.id,
                            'overallHealth': patientAgent.getStatus().overall}
        try:
            transferInfoDict = self.diagnosticPolicy.sendPatientTransferInfo(self,
                                                                             patientAgent,
                                                                             transferInfoDict)
        except MissingPatientRecordError:
            pass
        for tP in self.treatmentPolicies:
            try:
                transferInfoDict = tP.sendPatientTransferInfo(self, patientAgent,
                                                              transferInfoDict)
            except MissingPatientRecordError:
                pass
        return (0, desiredTier, self.abbrev, transferInfoDict)  # num bounces, tier, originating fac

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
        nBounces, tier, oldSenderAbbrev, transferInfoDict = payload  # @UnusedVariable
        if ward is not None:
            assert 'patientId' in transferInfoDict, 'Transfer info lacks patientId'
            pId = transferInfoDict['patientId']
            self.arrivingPatientTransferInfoDict[pId] = transferInfoDict.copy()
        else:
            pass
        # updated number of bounces and sender
        return (nBounces + 1, tier, self.abbrev, transferInfoDict), ward

    def handleBedRequestFate(self, ward, payload, timeNow):  # @UnusedVariable
        """
        This routine is called in the time slice of the Facility Manager when the manager
        receives the final response to a bed request it initiated. If the search for a bed
        failed, 'ward' will be None.
        """
        nBounces, tier, senderAbbrev, transferInfoDict = payload  # @UnusedVariable
        if ward is None:
            nFail = 1
            nSuccess = 0
        else:
            nFail = 0
            nSuccess = 1
        nh = self.getNoteHolder()
        bounceKey = CareTier.names[tier] + '_bounce_histo'
        transferKey = '%s_transfer' % senderAbbrev
        if bounceKey not in nh:
            nh.addNote({bounceKey: HistoVal([])})
        fromToKey = 'fromTo:%s:%s' % (senderAbbrev, tier)
        self.getNoteHolder().addNote({(CareTier.names[tier] + '_found'): nSuccess,
                                      (CareTier.names[tier] + '_notfound'): nFail,
                                      bounceKey: nBounces,
                                      transferKey: nSuccess,
                                      fromToKey: nSuccess})

    def diagnose(self, ward, patientId, patientStatus, oldDiagnosis, timeNow=None):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.
        """
        return self.diagnosticPolicy.diagnose(ward, patientId, patientStatus,
                                              oldDiagnosis, timeNow=timeNow)

    def prescribe(self, ward, patientId, patientDiagnosis, patientTreatment, 
                  modifierDct, timeNow=None):
        """
        This returns a tuple of either form (careTier, patientTreatment).
        modifierDct is a back-channel for down-stream communication and may be modified
        by the routines to which it is a parameter.
        """
        if patientDiagnosis.relocateFlag:
            modifierDct[pyrheabase.TierUpdateModKey.FORCE_MOVE] = True
        else:
            modifierDct = {}
        careTier = None
        for tP in self.treatmentPolicies:
            newTier, patientTreatment = tP.prescribe(ward, patientId,
                                                     patientDiagnosis, patientTreatment,
                                                     modifierDct, timeNow=timeNow)
            if careTier and (newTier != careTier):
                raise RuntimeError(('Treatment policies at %s prescribe different careTiers'
                                    % self.abbrev)
                                   (' for %s %s' % (patientDiagnosis, patientTreatment)))
            careTier = newTier
        return careTier, patientTreatment

    def getInitialOverallHealth(self, ward, timeNow):  # @UnusedVariable
        """
        To support multiple patient categories in the general population, this allows the
        facility to randomly assign an overall health category (HEALTHY, UNHEALTHY, ...)
        in accordance with a model for the facility population at simulation start time.
        Returns a PatientOverallHealth
        """
        return PatientOverallHealth.HEALTHY

    def diagnosisFromCareTier(self, careTier, overallHealth, timeNow):
        """
        If I look at a patient under a given care tier, what would I expect their diagnostic class
        to be?  This is used for patient initialization purposes.
        """
        return self.diagnosticPolicy.initializePatientDiagnosis(careTier, overallHealth, timeNow)

    def statusFromCareTier(self, careTier, overallHealth, patientDiagnosis, timeNow):  # @UnusedVariable
        return PatientStatus(overallHealth,
                             patientDiagnosis.diagClassA, 0,
                             patientDiagnosis.pthStatus, 0,
                             False, True, True, None,
#                            '', -1
                             )


    def getStatusChangeTree(self, patientAgent, modifierDct, startTime, timeNow):  # @UnusedVariable
        """
        Return a Bayes tree the traversal of which yields a patientStatus update.

        Elements of the tree are either a terminal patientStatus value or a 3-tuple
        of the form (k opt1 opt2), where opt1 and opt2 are likewise elements of the tree
        and k is a float in the range 0 to 1. The tree is evaluated by generating a random
        number each time a 3-tuple is encountered and choosing the first or second opt
        if the random value is <= or > k respectively.  When a value is encountered which
        is not a 3-tuple, that value is returned as the result of the traversal.
        """
        return BayesTree(PatientStatusSetter())

    def mapDescrFields(self, descr):
        """
        Descr is a record describing the facility, some fields of which contain
        category names in terms of 'descr' categories rather than 'impl' categories.
        Rebuild the record with 'impl' category names.
        """
        newDescr = descr.copy()
        for field in ['category']:
            newDescr[field] = self.categoryNameMapper(descr[field])
        for field in ['totalTransfersOut']:
            if field in descr:
                newD = {}
                for elt in descr[field]:
                    newCat = self.categoryNameMapper(elt['category'])
                    if newCat in newD:
                        val, prov = newD[newCat]
                        val += elt['count']['value']
                        prov += elt['count']['prov']
                        newD[newCat] = (val, prov)
                    else:
                        newD[newCat] = (elt['count']['value'], elt['count']['prov'])
                newL = []
                for cat, (val, prov) in newD.items():
                    newL.append({'category': cat, 'count': {'value': val, 'prov': prov}})
                newDescr[field] = newL
        return newDescr

    def scaleLOS(self, losModel, scaleLOS=1.0):
        if scaleLOS and scaleLOS != 1.0:
            assert losModel['pdf'] == 'lognorm(mu=$0,sigma=$1)', ('scaling PDF %s is not implemented'
                                                                  % losModel['pdf'])
            muPrime = losModel['parms'][0] + math.log(scaleLOS)
            return [muPrime, losModel['parms'][1]]
        else:
            return losModel['parms']


def decodeHistoryEntry(histEntry):
    return {"time": histEntry[0],
            "abbrev": histEntry[1],
            "category": histEntry[2],
            "careTier": histEntry[3]}

def buildTimeTupleList(agent, timeNow):
    """
    This routine is supposed to access the patient's history to produce a list of the form
    [(lengthOfStay, facAbbrev, facCategory, careTier)...] in most-recent-first order,
    *including* the current facility stay as the first entry.

    return [(2, 'foo', 'COMMUNITY', 'HOME'), (21, 'bar', 'HOSPITAL', 'HOSP')]
    """

    ret = []
    lastTime = timeNow if timeNow is not None else 0  # We can get entries from before the sim starts
    for histEntry in reversed(agent.agentHistory):
        t,abbrev,cat,tier = histEntry
        if t is None:
            t = 0  # We can get entries from before the sim starts
        ret.append((lastTime - t, abbrev, cat, tier))
        lastTime = t

    return ret


class PatientAgent(pyrheabase.PatientAgent):
    idCounters = defaultdict(int) # to provide a reliable identifier for each patient.

    def __init__(self, name, patch, ward, timeNow=0, debug=False):
        pyrheabase.PatientAgent.__init__(self, name, patch, ward, timeNow=timeNow, debug=debug)

        pOH = self.ward.fac.getInitialOverallHealth(ward, timeNow)
        self._diagnosis = self.ward.fac.diagnosisFromCareTier(self.ward.tier, pOH, timeNow)
        self._status = self.ward.fac.statusFromCareTier(self.ward.tier, pOH, self.getDiagnosis(),
                                                        timeNow)
        abbrev = self.ward.fac.abbrev
        self.id = (abbrev, PatientAgent.idCounters[abbrev])
        PatientAgent.idCounters[abbrev] += 1
        ward.fac.getPatientRecord(self.id, timeNow=timeNow)  # force creation of a blank record
        newTier, self._treatment = self.ward.fac.prescribe(self.ward,  # @UnusedVariable
                                                           self.id,
                                                           self.getDiagnosis(),
                                                           TREATMENT_DEFAULT,
                                                           {},
                                                           timeNow=timeNow)[0:2]

        self.lastUpdateTime = timeNow
        self.logger = logging.getLogger(__name__ + '.PatientAgent')
        self.agentHistory = []
        self.addHistoryEntry(self.ward, timeNow)

    @classmethod
    def allocateIds(cls, fac, count):
        cls.idCounters[fac.abbrev] += count

    def addHistoryEntry(self, ward, timeNow):
        self.agentHistory.append((timeNow, ward.fac.abbrev, ward.fac.category, ward.tier))

    def printSummary(self):
        print '%s as of %s' % (self.name, self.lastUpdateTime)
        print '  status: %s ' % str(self._status)
        print '  diagnosis: %s' % str(self._diagnosis)
        print '  treatment: %s' % self._treatment

    def getStatus(self):
        """Accessor for private status"""
        return self._status

    def setStatus(self, **kwargs):
        """
        keyword arguments are elements of PatientStatus, for example 'overall'.  The
        associated values must match the keyword type.
        """
        self._status = self._status._replace(**kwargs)

    def getPthStatus(self):
        """Accessor for private status pathogen info"""
        return self._status.pthStatus

    def getDiagnosis(self):
        """Accessor for private diagnosis"""
        return self._diagnosis

    def setDiagnosis(self, **kwargs):
        """
        keyword arguments are elements of PatientDiagnosis, for example 'overall'.  The
        associated values must match the keyword type.
        """
        self._diagnosis = self._diagnosis._replace(**kwargs)

    def getPthDiagnosis(self):
        """Accessor for private diagnosis pathogen info"""
        return self._diagnosis.pthStatus

    def setTreatment(self, **kwargs):
        """
        keyword arguments are elements of PatientTreatment, for example 'rehab'.
        The associated values must be boolean.
        """
        self._treatment = self._treatment._replace(**kwargs)

    def getTreatment(self, key):
        """
        key is one of the elements of PatientTreatment, for example 'rehab'.
        Returns a boolean for the state of that treatment element for this patient.
        """
        return self._treatment._asdict()[key]

    def getTreatmentProtocol(self):
        """
        Return the whole TreatmentProtocol instance for this patient.
        """
        return self._treatment

    def updateDiseaseState(self, treatment, facility, modifierDct, timeNow):  # @UnusedVariable
        """This should embody healing, community-acquired infection, etc."""
        dT = timeNow - self.lastUpdateTime
        previousStatus = self.getStatus()
        if dT > 0:  # moving from one ward to another can trigger two updates the same day
            treeL = []
            try:
                for src in [self.ward.fac, self.ward.iA]:
                    tree = src.getStatusChangeTree(self, modifierDct, self.lastUpdateTime, timeNow)
                    treeL.append(tree)
            except Exception, e:
                self.logger.critical('Got exception %s on patient %s (id %s) for from %s'
                                     % (str(e), self.name, self.id, str(src)))
                raise
            try:
                treeL = self.ward.iA.filterStatusChangeTrees(treeL, self,
                                                             self.lastUpdateTime, timeNow)
            except Exception, e:
                self.logger.critical('Got exception %s on patient %s (id %s) filtering trees'
                                     % (str(e), self.name, self.id))
                for tree in treeL:
                    print tree.tree
                    print tree.tagTree
                raise
            for tree in treeL:
                setter = tree.traverse()
                self._status = setter.set(self._status, timeNow)
            #print "Patient status at {0} is {1}".format(facility.abbrev, self._status.pthStatus == PthStatus.CLEAR)
            if (previousStatus.pthStatus != PthStatus.COLONIZED
                and self.getStatus().pthStatus == PthStatus.COLONIZED):
                #print "New Infection at {0}".format(self.ward.fac.abbrev)
                self.ward.miscCounters['newColonizationsSinceLastChecked'] += 1
                infectionLogger.info("%s newly colonized at %d in fac %s, tier %s, ward %d"%(
                    self.name, timeNow, self.ward.fac.abbrev, CareTier.names[self.ward.tier], self.ward.wardNum))
            if (previousStatus.pthStatus == PthStatus.COLONIZED
                and self.getStatus().pthStatus != PthStatus.COLONIZED):
                infectionLogger.info("%s decolonized at %d in fac %s, tier %s, ward %d"%(
                    self.name, timeNow, self.ward.fac.abbrev, CareTier.names[self.ward.tier], self.ward.wardNum))
            if self.getTreatment('creBundle'):
                self.ward.miscCounters['creBundlesHandedOut'] += 1

    def updateEverything(self, modifierDct, timeNow):
        """
        modifierDict is literally a back channel for downstream communication.  It may be updated
        in the routines to which it is a parameter.
        """
#         if self.lastUpdateTime != timeNow:
        if True:
            self.updateDiseaseState(self.getTreatmentProtocol(), self.ward.fac, modifierDct, timeNow)
            if self.getStatus().diagClassA == DiagClassA.DEATH:
                LOGGER.debug('%s died at %s at time %s', '%s_%s'%self.id, self.ward.fac.name, timeNow)
                infectionLogger.info("%s died at %d in fac %s, tier %s, ward %s with status %s"%(
                    self.name, timeNow, self.ward.fac.abbrev, CareTier.names[self.ward.tier],
                    self.ward.wardNum, PthStatus.names[self.getStatus().pthStatus]))
                return None
            self._diagnosis = self.ward.fac.diagnose(self.ward, self.id, self.getStatus(),
                                                     self.getDiagnosis(), timeNow=timeNow)
            newTier, self._treatment = self.ward.fac.prescribe(self.ward, self.id, self.getDiagnosis(),
                                                               self.getTreatmentProtocol(), modifierDct,
                                                               timeNow=timeNow)
            self.setStatus(justArrived=False)
            if self.debug:
                self.logger.debug('%s: status -> %s, diagnosis -> %s, treatment -> %s',
                                  self.name, self.getStatus(), self.getDiagnosis(),
                                  self.getTreatmentProtocol())
                self.logger.debug('%s: record at %s is %s',
                                  self.name, self.ward._name,
                                  self.ward.fac.getPatientRecord(self.id))
                self.logger.debug('%s: modifiers are %s',
                                  self.name,
                                  {pyrheabase.TierUpdateModKey.names[flg]: val
                                   for flg, val in modifierDct.items()})
            self.lastUpdateTime = timeNow
            return newTier
        else:
            return self.ward.tier

    def getPostArrivalPauseTime(self, timeNow):  # @UnusedVariable
        if self.ward.checkInterval > 1:
            return randint(0, self.ward.checkInterval-1)
        else:
            return 0

    def handleTierUpdate(self, modifierDict, timeNow):
        newTier = self.updateEverything(modifierDict, timeNow)
        return newTier

    def handleDeath(self, timeNow):
        assert self.getStatus().diagClassA == DiagClassA.DEATH, '%s is not dead?' % self.name
        self.ward.fac.noteHolder.addNote({"death": 1})
        if self.debug:
            self.logger.debug('Alas poor %s! %s', self.name, timeNow)
        facAddr = choice([tpl[1] for tpl in self.patch.serviceLookup('BirthQueue')])
        self.patch.launch(BirthMsg(self.name + '_birthMsg',
                                   self.patch,
                                   self.ward.fac.getMsgPayload(BirthMsg, self),
                                   facAddr),
                          timeNow)

    def getNewLocAddr(self, timeNow):
        newAddr = super(PatientAgent, self).getNewLocAddr(timeNow)
        if newAddr != self.locAddr:
            # It's about to move, so clear relocateFlag
            self.setDiagnosis(relocateFlag=False)
            self.setStatus(relocateFlag=False)
        return newAddr

    def getCandidateFacilityList(self, timeNow, modifierDct, newTier):
        cFL = self.ward.fac.getOrderedCandidateFacList(self, self.ward.tier, newTier, modifierDct, timeNow)
        return cFL

    def __getstate__(self):
        d = pyrheabase.PatientAgent.__getstate__(self)
        d['status'] = self._status
        d['diagnosis'] = self._diagnosis
        d['treatment'] = self._treatment
        d['lastUpdateTime'] = self.lastUpdateTime
        d['id'] = self.id
        d['agentHistory'] = self.agentHistory
        return d

    def __setstate__(self, d):
        pyrheabase.PatientAgent.__setstate__(self, d)
        self._status = d['status']
        self._diagnosis = d['diagnosis']
        self._treatment = d['treatment']
        self.lastUpdateTime = d['lastUpdateTime']
        self.id = d['id']
        self.agentHistory = d['agentHistory']
