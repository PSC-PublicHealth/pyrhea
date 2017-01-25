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

from random import randint, shuffle, choice
import pyrheabase
from phacsl.utils.collections.phacollections import enum, namedtuple
import logging
from phacsl.utils.notes.statval import HistoVal
from quilt.netinterface import GblAddr
from stats import BayesTree
from pathogenbase import PthStatus, defaultPthStatus, Pathogen
from collections import defaultdict

logger = logging.getLogger(__name__)

CareTier = enum('HOME', 'NURSING', 'LTAC', 'HOSP', 'ICU', 'VENT', 'SKILNRS')

PatientOverallHealth = enum('HEALTHY', 'FRAIL')

DiagClassA = enum('HEALTHY', 'NEEDSREHAB', 'NEEDSLTAC', 'SICK', 'VERYSICK', 'DEATH',
                  'NEEDSVENT', 'NEEDSSKILNRS')

TreatmentProtocol = enum('NORMAL', 'REHAB')  # for things like patient isolation

PatientStatus = namedtuple('PatientStatus',
                           ['overall',              # one of PatientOverallHealth
                            'diagClassA',           # one of DiagClassA
                            'startDateA',           # date diagClassA status was entered
                            'pthStatus',            # one of PthStatus
                            'startDatePth',         # date PthStatus status was entered
                            'relocateFlag',         # true if patient needs relocation
                            'canClear',             # true if patient can spontaneously clear infection
                            'homeAddr'              # GblAddr of patient's home tract or NH
                            ],
                           field_types=[PatientOverallHealth, DiagClassA, None, PthStatus, None,
                                        bool, bool, GblAddr])

PatientDiagnosis = namedtuple('PatientDiagnosis',
                              ['overall',              # one of PatientOverallHealth
                               'diagClassA',           # one of DiagClassA
                               'pthStatus',            # one of PthStatus
                               'relocateFlag'          # true if patient needs relocation
                               ],
                              field_types=[PatientOverallHealth, DiagClassA, PthStatus, bool])


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


class ForcedStateWard(Ward):
    def forceState(self, patientAgent, careTier, diagClassA):
        
        patientAgent._status =  patientAgent._status._replace(diagClassA=diagClassA)
        patientAgent._diagnosis = self.fac.diagnose(patientAgent._status, patientAgent._diagnosis)
        newTier, patientAgent._treatment = self.fac.prescribe(patientAgent._diagnosis,
                                                              patientAgent._treatment)[0:2]
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


class TransferDestinationPolicy(object):
    def __init__(self, patch, categoryNameMapper):
        self.patch = patch
        self.categoryNameMapper = categoryNameMapper

    def getOrderedCandidateFacList(self, facility, oldTier, newTier, timeNow):
        raise RuntimeError('Base TransferDestinationPolicy was called for %s' % facility.name)


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
        self.transferDestinationPolicy = transferDestinationPolicyClass(patch, self.categoryNameMapper)
        
    def __str__(self):
        return '<%s>' % self.name

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
                     patientAgent._status.overall == PatientOverallHealth.FRAIL,
                     patientAgent.name),
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
            patientID, tier, isFrail, patientName = myPayload
            if patientID in self.patientDataDict:
                logger.debug('Patient %s has returned to %s' % (patientID, self.name))
                patientRec = self.patientDataDict[patientID]
                patientRec.prevVisits += 1
                patientRec.arrivalDate = timeNow
                patientRec.departureDate = None
            else:
#                patientRec = Facility.PatientRecord(patientID, timeNow, isFrail)
                patientRec = PatientRecord(patientID, timeNow, isFrail)
                self.patientDataDict[patientID] = patientRec
            if timeNow != 0:  # Exclude initial populations
                nh = self.getNoteHolder()
                if nh:
                    nh.addNote({(CareTier.names[tier] + '_arrivals'): 1})
        elif issubclass(msgType, pyrheabase.DepartureMsg):
            myPayload, innerPayload = payload
            timeNow = super(Facility, self).handleIncomingMsg(msgType, innerPayload, timeNow)
            patientID, tier, isFrail, patientName = myPayload
            if patientID not in self.patientDataDict:
                logger.error('%s has no record of patient %s' % (self.name, patientID))
            patientRec = self.patientDataDict[patientID]
            patientRec.departureDate = timeNow
            #if patientRec.arrivalDate != 0:  # exclude initial populations
            if True:  # include initial populations
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
        receives the final response to a bed request it initiated. If the search for a bed
        failed, 'ward' will be None.
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
        fromToKey = 'fromTo:%s:%s' % (senderAbbrev, tier)
        self.getNoteHolder().addNote({(CareTier.names[tier] + '_found'): nSuccess,
                                      (CareTier.names[tier] + '_notfound'): nFail,
                                      bounceKey: nBounces,
                                      transferKey: nSuccess,
                                      fromToKey: nSuccess})

    def diagnose(self, patientStatus, oldDiagnosis):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.
        """
        return PatientDiagnosis(patientStatus.overall,
                                patientStatus.diagClassA,
                                patientStatus.pthStatus,
                                patientStatus.relocateFlag)

    def prescribe(self, patientDiagnosis, patientTreatment):
        """
        This returns a tuple of either form (careTier, patientTreatment) or
        (careTier, patientTreatment, modifierList)
        """
        if patientDiagnosis.relocateFlag:
            modifierList = [pyrheabase.TierUpdateModFlag.FORCE_MOVE]
        else:
            modifierList = []
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
                return (CareTier.HOME, TreatmentProtocol.NORMAL, modifierList)
            elif patientDiagnosis.overall == PatientOverallHealth.FRAIL:
                return (CareTier.NURSING, TreatmentProtocol.NORMAL, modifierList)
            else:
                raise RuntimeError('Unknown overall health %s' % str(patientDiagnosis.overall))
        if patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
            return (CareTier.NURSING, TreatmentProtocol.REHAB, modifierList)
        elif patientDiagnosis.diagClassA == DiagClassA.SICK:
            return (CareTier.HOSP, TreatmentProtocol.NORMAL, modifierList)
        elif patientDiagnosis.diagClassA == DiagClassA.VERYSICK:
            return (CareTier.ICU, TreatmentProtocol.NORMAL, modifierList)
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSLTAC:
            return (CareTier.LTAC, TreatmentProtocol.NORMAL, modifierList)
        elif patientDiagnosis.diagClassA == DiagClassA.DEATH:
            return (None, TreatmentProtocol.NORMAL, modifierList)
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSVENT:
            return (CareTier.VENT, TreatmentProtocol.NORMAL, modifierList)
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSSKILNRS:
            return (CareTier.SKILNRS, TreatmentProtocol.NORMAL, modifierList)
        else:
            raise RuntimeError('Unknown DiagClassA %s' % str(patientDiagnosis.diagClassA))

    def diagnosisFromCareTier(self, careTier):
        """
        If I look at a patient under a given care tier, what would I expect their diagnostic class
        to be?  This is used for patient initialization purposes.
        """
        if careTier == CareTier.HOME:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.HEALTHY, defaultPthStatus, False)
        elif careTier == CareTier.NURSING:
            return PatientDiagnosis(PatientOverallHealth.FRAIL,
                                    DiagClassA.HEALTHY, defaultPthStatus, False)
        elif careTier == CareTier.LTAC:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.NEEDSLTAC, defaultPthStatus, False)
        elif careTier == CareTier.HOSP:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.SICK, defaultPthStatus, False)
        elif careTier == CareTier.ICU:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.VERYSICK, defaultPthStatus, False)
        elif careTier == CareTier.VENT:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.NEEDSVENT, defaultPthStatus, False)
        elif careTier == CareTier.SKILNRS:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.NEEDSSKILNRS, defaultPthStatus, False)
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
    

class PatientAgent(pyrheabase.PatientAgent):
    idCounters = defaultdict(int) # to provide a reliable identifier for each patient.

    def __init__(self, name, patch, ward, timeNow=0, debug=False):
        pyrheabase.PatientAgent.__init__(self, name, patch, ward, timeNow=timeNow, debug=debug)
        self._diagnosis = self.ward.fac.diagnosisFromCareTier(self.ward.tier)
        self._status = PatientStatus(PatientOverallHealth.HEALTHY,
                                     self._diagnosis.diagClassA, 0,
                                     self._diagnosis.pthStatus, 0,
                                     False, True, None)
        newTier, self._treatment = self.ward.fac.prescribe(self._diagnosis,  # @UnusedVariable
                                                           TreatmentProtocol.NORMAL)[0:2]
        self.lastUpdateTime = timeNow
        abbrev = self.ward.fac.abbrev
        self.id = (abbrev, PatientAgent.idCounters[abbrev])
        PatientAgent.idCounters[abbrev] += 1
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
            treeL = []
            try:
                for src in [self.ward.fac, self.ward.iA]:
                    tree = src.getStatusChangeTree(self._status, self.ward, self._treatment,
                                                   self.lastUpdateTime, timeNow)
                    treeL.append(tree)
            except Exception, e:
                print ('Got exception %s on patient %s tree for %s' %
                       (str(e), self.name, str(src)))
                self.logger.critical('Got exception %s on patient %s (id %s) for tree %s'
                                     % (str(e), self.name, self.id, str(src)))
                raise
            try:
                treeL = self.ward.iA.filterStatusChangeTrees(treeL, self._status, self.ward, self._treatment,
                                                             self.lastUpdateTime, timeNow)
            except Exception, e:
                print ('Got exception %s on patient %s filtering trees' %
                       (str(e), self.name))
                self.logger.critical('Got exception %s on patient %s (id %s) filtering trees'
                                     % (str(e), self.name, self.id))
                for tree in treeL:
                    print tree.tree
                    print tree.tagTree
                raise
            for tree in treeL:
                setter = tree.traverse()
                self._status = setter.set(self._status, timeNow)

    def updateEverything(self, timeNow):
        self.updateDiseaseState(self._treatment, self.ward.fac, timeNow)
        if self._status.diagClassA == DiagClassA.DEATH:
            return None, []
        self._diagnosis = self.ward.fac.diagnose(self._status, self._diagnosis)
        tpl = self.ward.fac.prescribe(self._diagnosis, self._treatment)
        newTier, self._treatment = tpl[0:2]
        if len(tpl) == 3:
            modifierList = tpl[2]
        else:
            modifierList = []
        if self.debug:
            self.logger.debug('%s: status -> %s, diagnosis -> %s, treatment -> %s'
                              % (self.name, self._status, self._diagnosis,
                                 TreatmentProtocol.names[self._treatment]))
            self.logger.debug('%s: modifiers are %s',
                              self.name,
                              [pyrheabase.TierUpdateModFlag.names[flg] for flg in modifierList])
        self.lastUpdateTime = timeNow
        return newTier, modifierList

    def getPostArrivalPauseTime(self, timeNow):
        if self.ward.checkInterval > 1:
            return randint(0, self.ward.checkInterval-1)
        else:
            return 0

    def handleTierUpdate(self, timeNow):
        newTier, modifierList = self.updateEverything(timeNow)
        return newTier, modifierList

    def handleDeath(self, timeNow):
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

    def getNewLocAddr(self, timeNow):
        newAddr = super(PatientAgent, self).getNewLocAddr(timeNow)
        if newAddr != self.locAddr:
            # It's about to move, so clear relocateFlag
            self._diagnosis = self._diagnosis._replace(relocateFlag=False)
            self._status = self._status._replace(relocateFlag=False)
        return newAddr
    
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
