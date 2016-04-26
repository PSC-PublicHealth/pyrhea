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

from random import randint, shuffle, choice
import pyrheabase
from phacsl.utils.collections.phacollections import enum, namedtuple
import logging
from phacsl.utils.notes.statval import HistoVal
from stats import BayesTree
from pathogenbase import Status as PthStatus, defaultStatus as defaultPthStatus
from pathogenbase import Pathogen  # @UnusedImport

logger = logging.getLogger(__name__)

CareTier = enum('HOME', 'NURSING', 'LTAC', 'HOSP', 'ICU')

PatientOverallHealth = enum('HEALTHY', 'FRAIL')

DiagClassA = enum('HEALTHY', 'NEEDSREHAB', 'NEEDSLTAC', 'SICK', 'VERYSICK', 'DEATH')

TreatmentProtocol = enum('NORMAL', 'REHAB')  # for things like patient isolation

PatientStatus = namedtuple('PatientStatus',
                           ['overall',              # one of patientOverallHealth
                            'diagClassA',           # one of DiagClassA
                            'startDateA',           # date diagClassA status was entered
                            'pthStatus',            # one of PthStatus
                            'startDatePth'          # date PthStatus status was entered
                            ],
                           field_types=[PatientOverallHealth, DiagClassA, None, PthStatus, None])

PatientDiagnosis = namedtuple('PatientDiagnosis',
                              ['overall',              # one of PatientOverallHealth
                               'diagClassA',           # one of DiagClassA
                               'pthStatus',            # one of PthStatus
                               ],
                              field_types=[PatientOverallHealth, DiagClassA, PthStatus])


class PatientStatusSetter(object):
    def __init__(self):
        pass

    def set(self, patientStatus, timeNow):
        """This base class just returns a copy"""
        return PatientStatus._make(patientStatus)

    def __str__(self):
        return 'PatientStatusSetter()'


class PthStatusSetter(PatientStatusSetter):
    def __init__(self, newPthStatus):
        self.newPthStatus = newPthStatus

    def set(self, patientStatus, timeNow):
        return patientStatus._replace(pthStatus=self.newPthStatus, startDateA=timeNow)

    def __str__(self):
        return 'PatientStatusSetter(pthStatus <- %s)' % PthStatus.names[self.newPthStatus]


class Ward(pyrheabase.Ward):
    def __init__(self, name, patch, tier, nBeds):
        pyrheabase.Ward.__init__(self, name, patch, tier, nBeds)
        self.checkInterval = 1  # check health daily
        self.iA = None  # infectious agent

    def getPatientList(self):
        return self._lockingAgentList[1:self._nLocks]

    def setInfectiousAgent(self, iA):
        self.iA = iA

    def initializePatientPthState(self):
        for p in self.getPatientList():
            self.iA.initializePatientState(p)


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


tierToQueueMap = {CareTier.HOME: HOMEQueue,
                  CareTier.NURSING: NURSINGQueue,
                  CareTier.LTAC: LTACQueue,
                  CareTier.HOSP: HOSPQueue,
                  CareTier.ICU: ICUQueue}


class TransferDestinationPolicy(object):
    def __init__(self, patch):
        self.patch = patch

    def getOrderedCandidateFacList(self, facility, oldTier, newTier, timeNow):
        queueClass = tierToQueueMap[newTier]
        facAddrList = [tpl[1] for tpl in self.patch.serviceLookup(queueClass.__name__)]
        shuffle(facAddrList)
        return facAddrList


class Facility(pyrheabase.Facility):
    class PatientRecord(object):
        def __init__(self, patientID, arrivalDate, isFrail):
            self.patientID = patientID
            self.arrivalDate = arrivalDate
            self.departureDate = None
            self.prevVisits = 0
            self.isFrail = isFrail

        def __str__(self):
            return '<patient %s, %s -> %s, %s>' % (self.patientID,
                                                   self.arrivalDate,
                                                   self.departureDate,
                                                   'Frail' if self.isFrail else 'Healthy')

    def __init__(self, name, descr, patch, reqQueueClasses=None, policyClasses=None):
        pyrheabase.Facility.__init__(self, name, patch, managerClass=FacilityManager,
                                     reqQueueClasses=reqQueueClasses)
        self.category = descr['category']
        self.abbrev = descr['abbrev']
        if 'longitude' in descr and 'latitude' in descr:
            self.coords = (descr['longitude'], descr['latitude'])
        else:
            self.coords = (None, None)
        self.noteHolder = None
        self.idCounter = 0
        self.patientDataDict = {}
        transferDestinationPolicyClass = TransferDestinationPolicy
        if policyClasses is not None:
            for pC in policyClasses:
                if issubclass(pC, TransferDestinationPolicy):
                    transferDestinationPolicyClass = pC
        self.transferDestinationPolicy = transferDestinationPolicyClass(patch)

    def setNoteHolder(self, noteHolder):
        self.noteHolder = noteHolder
        self.noteHolder.addNote({'name': self.name})

    def getNoteHolder(self):
        return self.noteHolder

    def getOrderedCandidateFacList(self, oldTier, newTier, timeNow):
        return self.transferDestinationPolicy.getOrderedCandidateFacList(self,
                                                                         oldTier, newTier,
                                                                         timeNow)

    def getMsgPayload(self, msgType, patientAgent):
        if issubclass(msgType, (pyrheabase.ArrivalMsg, pyrheabase.DepartureMsg)):
            innerPayload = super(Facility, self).getMsgPayload(msgType, patientAgent)
            return ((patientAgent.id, patientAgent.ward.tier,
                     patientAgent._status.overall == PatientOverallHealth.FRAIL),
                    innerPayload)
        elif issubclass(msgType, BirthMsg):
            return patientAgent._status.overall
        else:
            raise RuntimeError('%s: payload request for unknown message type %s'
                               % (self.name, msgType.__name__))

    def handleIncomingMsg(self, msgType, payload, timeNow):
        if issubclass(msgType, pyrheabase.ArrivalMsg):
            myPayload, innerPayload = payload
            timeNow = super(Facility, self).handleIncomingMsg(msgType, innerPayload, timeNow)
            patientID, tier, isFrail = myPayload
            if patientID in self.patientDataDict:
                logger.info('Patient %s has returned to %s' % (patientID, self.name))
                patientRec = self.patientDataDict[patientID]
                patientRec.prevVisits += 1
                patientRec.arrivalDate = timeNow
                patientRec.departureDate = None
            else:
                patientRec = Facility.PatientRecord(patientID, timeNow, isFrail)
                self.patientDataDict[patientID] = patientRec
            if timeNow != 0:  # Exclude initial populations
                nh = self.getNoteHolder()
                if nh:
                    nh.addNote({(CareTier.names[tier] + '_arrivals'): 1})
        elif issubclass(msgType, pyrheabase.DepartureMsg):
            myPayload, innerPayload = payload
            timeNow = super(Facility, self).handleIncomingMsg(msgType, innerPayload, timeNow)
            patientID, tier, isFrail = myPayload
            if patientID not in self.patientDataDict:
                logger.error('%s has no record of patient %s' % (self.name, patientID))
            patientRec = self.patientDataDict[patientID]
            patientRec.departureDate = timeNow
            if patientRec.arrivalDate != 0:  # exclude initial populations
                lengthOfStay = timeNow - patientRec.arrivalDate
                losKey = (CareTier.names[tier] + '_LOS')
                nh = self.getNoteHolder()
                if losKey not in nh:
                    nh.addNote({losKey: HistoVal([])})
                nh.addNote({(CareTier.names[tier] + '_departures'): 1, losKey: lengthOfStay})
        elif issubclass(msgType, BirthMsg):
            ward = self.manager.allocateAvailableBed(CareTier.HOME)
            assert ward is not None, 'Ran out of beds with birth in %s!' % self.name
            a = PatientAgent('PatientAgent_%s_birth' % ward._name, self.manager.patch, ward)
            a._status = a._status._replace(overall=payload)
            ward.lock(a)
            self.handleIncomingMsg(pyrheabase.ArrivalMsg,
                                   self.getMsgPayload(pyrheabase.ArrivalMsg, a),
                                   timeNow)
            self.manager.patch.launch(a, timeNow)
            if timeNow != 0:  # exclude initial populations
                nh = self.getNoteHolder()
                nh.addNote({'births': 1})
        else:
            raise RuntimeError('%s: got unknown message type %s' % (self.name,
                                                                    msgType .__name__))
        return timeNow

    def getBedRequestPayload(self, patientAgent, desiredTier):
        return (0, desiredTier, self.abbrev)  # number of bounces, tier, originating fac

    def handleBedRequestResponse(self, ward, payload, timeNow):
        """
        This routine is called in the time slice of the Facility Manager when the manager
        responds to a request for a bed.  If the request was denied, 'ward' will be None.
        The return value of this message becomes the new payload.
        """
        nBounces, tier, oldSenderAbbrev = payload  # @UnusedVariable
        return (nBounces + 1, tier, self.abbrev)  # updated number of bounces and sender

    def handleBedRequestFate(self, ward, payload, timeNow):
        """
        This routine is called in the time slice of the Facility Manager when the manager
        responds to a request for a bed.  If the search for a bed failed, 'ward' will be None.
        """
        nBounces, tier, senderAbbrev = payload
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
        self.getNoteHolder().addNote({(CareTier.names[tier] + '_found'): nSuccess,
                                      (CareTier.names[tier] + '_notfound'): nFail,
                                      bounceKey: nBounces,
                                      transferKey: nSuccess})

    def diagnose(self, patientStatus, oldDiagnosis):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.
        """
        return PatientDiagnosis(patientStatus.overall,
                                patientStatus.diagClassA,
                                patientStatus.pthStatus)

    def prescribe(self, patientDiagnosis, patientTreatment):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
                return (CareTier.HOME, TreatmentProtocol.NORMAL)
            elif patientDiagnosis.overall == PatientOverallHealth.FRAIL:
                return (CareTier.NURSING, TreatmentProtocol.NORMAL)
            else:
                raise RuntimeError('Unknown overall health %s' % str(patientDiagnosis.overall))
        if patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
            return (CareTier.NURSING, TreatmentProtocol.REHAB)
        elif patientDiagnosis.diagClassA == DiagClassA.SICK:
            return (CareTier.HOSP, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.VERYSICK:
            return (CareTier.ICU, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.DEATH:
            return (None, TreatmentProtocol.NORMAL)
        else:
            raise RuntimeError('Unknown DiagClassA %s' % str(patientDiagnosis.diagClassA))

    def diagnosisFromCareTier(self, careTier):
        """
        If I look at a patient under a given care tier, what would I expect their diagnostic class
        to be?  This is used for patient initialization purposes.
        """
        if careTier == CareTier.HOME:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.HEALTHY, defaultPthStatus)
        elif careTier == CareTier.NURSING:
            return PatientDiagnosis(PatientOverallHealth.FRAIL,
                                    DiagClassA.HEALTHY, defaultPthStatus)
        elif careTier == CareTier.LTAC:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.NEEDSLTAC, defaultPthStatus)
        elif careTier == CareTier.HOSP:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.SICK, defaultPthStatus)
        elif careTier == CareTier.ICU:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.VERYSICK, defaultPthStatus)
        else:
            raise RuntimeError('Unknown care tier %s' % careTier)

    def getStatusChangeTree(self, patientStatus, ward, treatment, startTime, timeNow):
        """
        Return a Bayes tree the traversal of which yields a patientStatus.

        Elements of the tree are either a terminal patientStatus value or a 3-tuple
        of the form (k opt1 opt2), where opt1 and opt2 are likewise elements of the tree
        and k is a float in the range 0 to 1. The tree is evaluated by generating a random
        number each time a 3-tuple is encountered and choosing the first or second opt
        if the random value is <= or > k respectively.  When a value is encountered which
        is not a 3-tuple, that value is returned as the result of the traversal.
        """
        return BayesTree(PatientStatusSetter())


class PatientAgent(pyrheabase.PatientAgent):
    idCounter = 0  # to provide a reliable identifier for each patient

    def __init__(self, name, patch, ward, timeNow=0, debug=False):
        pyrheabase.PatientAgent.__init__(self, name, patch, ward, timeNow=timeNow, debug=debug)
        self._diagnosis = self.ward.fac.diagnosisFromCareTier(self.ward.tier)
        self._status = PatientStatus(PatientOverallHealth.HEALTHY,
                                     self._diagnosis.diagClassA, 0,
                                     self._diagnosis.pthStatus, 0)
        newTier, self._treatment = self.ward.fac.prescribe(self._diagnosis,  # @UnusedVariable
                                                           TreatmentProtocol.NORMAL)
        self.lastUpdateTime = timeNow
        self.id = (patch.patchId, PatientAgent.idCounter)
        PatientAgent.idCounter += 1
        self.logger = logging.getLogger(__name__ + '.PatientAgent')

    def printSummary(self):
        print '%s as of %s' % (self.name, self.lastUpdateTime)
        print '  status: %s ' % str(self._status)
        print '  diagnosis: %s' % str(self._diagnosis)
        print '  treatment: %s' % TreatmentProtocol.names[self._treatment]

    def updateDiseaseState(self, treatment, facility, timeNow):
        """This should embody healing, community-acquired infection, etc."""
        dT = timeNow - self.lastUpdateTime
        if dT > 0:  # moving from one ward to another can trigger two updates the same day
            for src in [self.ward.fac, self.ward.iA]:
                try:
                    tree = src.getStatusChangeTree(self._status, self.ward, self._treatment,
                                                   self.lastUpdateTime, timeNow)
                    setter = tree.traverse()
                    self._status = setter.set(self._status, timeNow)
                except Exception, e:
                    print 'Got exception %s on patient %s' % (str(e), self.name)
                    self.logger.critical('Got exception %s on patient %s (id %s)'
                                         % (str(e), self.name, self.id))
                    raise

    def updateEverything(self, timeNow):
        self.updateDiseaseState(self._treatment, self.ward.fac, timeNow)
        if self._status.diagClassA == DiagClassA.DEATH:
            return None
        self._diagnosis = self.ward.fac.diagnose(self._status, self._diagnosis)
        newTier, self._treatment = self.ward.fac.prescribe(self._diagnosis, self._treatment)
        if self.debug:
            self.logger.debug('%s: status -> %s, diagnosis -> %s, treatment -> %s'
                              % (self.name, self._status, self._diagnosis, self._treatment))
        self.lastUpdateTime = timeNow
        return newTier

    def getPostArrivalPauseTime(self, timeNow):
        if self.ward.checkInterval > 1:
            return randint(0, self.ward.checkInterval-1)
        else:
            return 0

    def handleTierUpdate(self, timeNow):
        return self.updateEverything(timeNow)

    def handlePatientDeath(self, timeNow):
        assert self._status.diagClassA == DiagClassA.DEATH, '%s is not dead?' % self.name
        self.ward.fac.noteHolder.addNote({"death": 1})
        if self.debug:
            self.logger.debug('Alas poor %s! %s' % (self.name, timeNow))
        facAddr = choice([tpl[1] for tpl in self.patch.serviceLookup('BirthQueue')])
        self.patch.launch(BirthMsg(self.name + '_birthMsg',
                                   self.patch,
                                   self.ward.fac.getMsgPayload(BirthMsg, self),
                                   facAddr),
                          timeNow)

    def getCandidateFacilityList(self, timeNow, newTier):
        return self.ward.fac.getOrderedCandidateFacList(self.ward.tier, newTier, timeNow)

    def __getstate__(self):
        d = pyrheabase.PatientAgent.__getstate__(self)
        d['status'] = self._status
        d['diagnosis'] = self._diagnosis
        d['treatment'] = self._treatment
        d['lastUpdateTime'] = self.lastUpdateTime
        d['id'] = self.id
        return d

    def __setstate__(self, d):
        pyrheabase.PatientAgent.__setstate__(self, d)
        self._status = d['status']
        self._diagnosis = d['diagnosis']
        self._treatment = d['treatment']
        self.lastUpdateTime = d['lastUpdateTime']
        self.id = d['id']
