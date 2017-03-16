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
import cPickle as pickle
from random import random

logger = logging.getLogger(__name__)

CareTier = enum('HOME', 'NURSING', 'LTAC', 'HOSP', 'ICU', 'VENT', 'SKILNRS')

PatientOverallHealth = enum('HEALTHY', 'FRAIL')

DiagClassA = enum('HEALTHY', 'NEEDSREHAB', 'NEEDSLTAC', 'SICK', 'VERYSICK', 'DEATH',
                  'NEEDSVENT', 'NEEDSSKILNRS')

TreatmentProtocol = namedtuple('TreatmentProtocol',
                               ['rehab',
                                'contactPrecautions',
                                'creBundle'
                                ],
                               field_types=[bool, bool, bool])

TREATMENT_DEFAULT = TreatmentProtocol(rehab=False, contactPrecautions=False, creBundle=False)

PatientStatus = namedtuple('PatientStatus',
                           ['overall',              # one of PatientOverallHealth
                            'diagClassA',           # one of DiagClassA
                            'startDateA',           # date diagClassA status was entered
                            'pthStatus',            # one of PthStatus
                            'startDatePth',         # date PthStatus status was entered
                            'relocateFlag',         # true if patient needs relocation
                            'justArrived',          # true on patient's first day in new location
                            'canClear',             # true if patient can spontaneously clear infection
                            'homeAddr'              # GblAddr of patient's home tract or NH
                            ],
                           field_types=[PatientOverallHealth, DiagClassA, None, PthStatus, None,
                                        bool, bool, bool, GblAddr])

PatientDiagnosis = namedtuple('PatientDiagnosis',
                              ['overall',              # one of PatientOverallHealth
                               'diagClassA',           # one of DiagClassA
                               'startDateA',           # Date diagClassA was entered
                               'pthStatus',            # one of PthStatus
                               'relocateFlag'          # true if patient needs relocation
                               ],
                              field_types=[PatientOverallHealth, DiagClassA, None, PthStatus, bool])


class PatientStatusSetter(object):
    def __init__(self):
        pass

    def set(self, patientStatus, timeNow):
        """This base class just returns a copy"""
        return PatientStatus._make(patientStatus)

    def __str__(self):
        return 'PatientStatusSetter()'


class ClassASetter(PatientStatusSetter):
    def __init__(self, newClassA):
        super(ClassASetter, self).__init__()
        self.newClassA = newClassA

    def set(self, patientStatus, timeNow):
        return (patientStatus._replace(diagClassA=self.newClassA, startDateA=timeNow)
                ._replace(relocateFlag=True))

    def __str__(self):
        return 'PatientStatusSetter(classA <- %s)' % DiagClassA.names[self.newClassA]

    def __repr__(self):
        return 'PatientStatusSetter(classA <- %s)' % DiagClassA.names[self.newClassA]


class PthStatusSetter(PatientStatusSetter):
    def __init__(self, newPthStatus):
        self.newPthStatus = newPthStatus

    def set(self, patientStatus, timeNow):
        return patientStatus._replace(pthStatus=self.newPthStatus, startDatePth=timeNow)

    def __str__(self):
        return 'PatientStatusSetter(pthStatus <- %s)' % PthStatus.names[self.newPthStatus]


class Ward(pyrheabase.Ward):
    def __init__(self, name, patch, tier, nBeds):
        pyrheabase.Ward.__init__(self, name, patch, tier, nBeds)
        self.checkInterval = 1  # check health daily
        self.iA = None  # infectious agent
        self.miscCounters = defaultdict(lambda: 0)
        ### Shawn Sucks and is paranoid
        self.miscCounters['passiveDaysOnCP'] = 0
        self.miscCounters['swabDaysOnCP'] = 0
        self.miscCounters['otherDaysOnCP'] = 0
        self.miscCounters['xdroDaysOnCP'] = 0

    def getPatientList(self):
        return self._lockingAgentList[:self._nLocks]

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
        patientAgent._status = patientAgent._status._replace(justArrived=True)
        self.miscCounters['arrivals'] += 1
        if patientAgent._status.pthStatus == PthStatus.COLONIZED:
            self.miscCounters['creArrivals'] += 1
        for tP in self.fac.treatmentPolicies:
            tP.handlePatientArrival(self, patientAgent, timeNow)


    def handlePatientDeparture(self, patientAgent, timeNow):
        """An opportunity for derived classes to customize the departure processing of patients"""
        for tP in self.fac.treatmentPolicies:
            tP.handlePatientDeparture(self, patientAgent, timeNow)
        self.miscCounters['departures'] += 1


class ForcedStateWard(Ward):
    def forceState(self, patientAgent, careTier, diagClassA):
        
        patientAgent._status =  patientAgent._status._replace(diagClassA=diagClassA)
        patientAgent._diagnosis = self.fac.diagnose(patientAgent, patientAgent._diagnosis)
        newTier, patientAgent._treatment = self.fac.prescribe(self, patientAgent, patientAgent._diagnosis,
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


class Policy(object):
    def __init__(self, patch, categoryNameMapper):
        self.patch = patch
        self.categoryNameMapper = categoryNameMapper


class DiagnosticPolicy(Policy):
    def diagnose(self, patient, oldDiagnosis, fac=None):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.
        """
        patientStatus = patient._status
        return PatientDiagnosis(patientStatus.overall,
                                patientStatus.diagClassA,
                                patientStatus.startDateA,
                                PthStatus.CLEAR,
                                patientStatus.relocateFlag)
    
    def initializePatientDiagnosis(self, ward, timeNow):
        careTier = ward.tier
        if careTier == CareTier.HOME:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.HEALTHY, timeNow, defaultPthStatus, False)
        elif careTier == CareTier.NURSING:
            return PatientDiagnosis(PatientOverallHealth.FRAIL,
                                    DiagClassA.HEALTHY, timeNow, defaultPthStatus, False)
        elif careTier == CareTier.LTAC:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.NEEDSLTAC, timeNow, defaultPthStatus, False)
        elif careTier == CareTier.HOSP:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.SICK, timeNow, defaultPthStatus, False)
        elif careTier == CareTier.ICU:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.VERYSICK, timeNow, defaultPthStatus, False)
        elif careTier == CareTier.VENT:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.NEEDSVENT, timeNow, defaultPthStatus, False)
        elif careTier == CareTier.SKILNRS:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.NEEDSSKILNRS, timeNow, defaultPthStatus, False)
        else:
            raise RuntimeError('Unknown care tier %s' % careTier)

    def setValue(self, key, val):
        """
        Setting values may be useful for changing phases in a scenario, for example. The
        values that can be set are treatment-specific; attempting to set an incorrect value
        is an error.
        
        The base class doesn't know how to set any values.
        """
        raise RuntimeError('Class %s does not know how to set the value %s'
                           % (type(self).__name__, key))


class TreatmentPolicy(Policy):
    def initializePatientTreatment(self, ward, patient):
        """
        This is called on patients at time zero, when they are first assigned to the
        ward in which they start the simulation.
        """
        raise RuntimeError('Base TreatmentPolicy was called for %s' % ward._name)

    def handlePatientArrival(self, ward, patient, timeNow):
        """
        This is called on patients when they arrive at a ward.
        """
        raise RuntimeError('Base TreatmentPolicy was called for %s' % ward._name)

    def handlePatientDeparture(self, ward, patient, timeNow):
        """
        This is called on patients when they depart from a ward.
        """
        raise RuntimeError('Base TreatmentPolicy was called for %s' % ward._name)

    def prescribe(self, ward, patient, patientDiagnosis, patientTreatment, modifierList):
        """
        This returns a tuple of form (careTier, patientTreatment).
        modifierList is for functional modifiers, like pyrheabase.TierUpdateModFlag.FORCE_MOVE,
        and is not generally relevant to the decisions made by this method.
        """
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
                return (CareTier.HOME, patientTreatment._replace(rehab=False))
            elif patientDiagnosis.overall == PatientOverallHealth.FRAIL:
                return (CareTier.NURSING, patientTreatment._replace(rehab=False))
            else:
                raise RuntimeError('Unknown overall health %s' % str(patientDiagnosis.overall))
        if patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
            newTreatment = patientTreatment._replace(rehab=True)
            return (CareTier.NURSING, newTreatment)
        elif patientDiagnosis.diagClassA == DiagClassA.SICK:
            newTreatment = patientTreatment._replace(rehab=False)
            return (CareTier.HOSP, newTreatment)
        elif patientDiagnosis.diagClassA == DiagClassA.VERYSICK:
            newTreatment = patientTreatment._replace(rehab=False)
            return (CareTier.ICU, newTreatment)
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSLTAC:
            newTreatment = patientTreatment._replace(rehab=False)
            return (CareTier.LTAC, newTreatment)
        elif patientDiagnosis.diagClassA == DiagClassA.DEATH:
            newTreatment = patientTreatment._replace(rehab=False)
            return (None, newTreatment)
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSVENT:
            newTreatment = patientTreatment._replace(rehab=False)
            return (CareTier.VENT, newTreatment)
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSSKILNRS:
            newTreatment = patientTreatment._replace(rehab=False)
            return (CareTier.SKILNRS, newTreatment)
        else:
            raise RuntimeError('Unknown DiagClassA %s' % str(patientDiagnosis.diagClassA))

    def getTransmissionFromMultiplier(self, careTier, **kwargs):
        """
        If the treatment elements in **kwargs have the given boolean values (e.g. rehab=True),
        return the scale factor by which the transmission coefficient tau is multiplied when
        the patient with this treatment is the source of the transmission.
        """
        raise RuntimeError('Base TreatmentPolicy was called.')        

    def getTransmissionToMultiplier(self, careTier, **kwargs):
        """
        If the treatment elements in **kwargs have the given boolean values (e.g. rehab=True),
        return the scale factor by which the transmission coefficient tau is multiplied when
        the patient with this treatment is the recipient of the transmission.
        """
        raise RuntimeError('Base TreatmentPolicy was called.')        

    @classmethod
    def getRelativeProb(cls, pthStatus, fromTier, toTier):
        """
        If the probability of transfer from fromTier to toTier of a patient at
        PthStatus.CLEAR is P, and the probability for a patient at the given PthStatus
        is kP, this routine returns the value of k.  Note that kP must still be less than 1.0,
        so there is an implied upper bound of P of 1.0/k.
        """
        return 1.0

    @classmethod
    def getEstimatedPrevalence(cls, pthStatus, abbrev, category, tier):
        """
        The return value provides an estimate prevalence (0.0 <= return val <= 1.0)
        of the pathogen for the given pathogen status at the facility named by abbrev,
        of the given category, at the given care tier.  One use of this value is to
        help keep patient flows in the correct range while rescaling the flow of a
        particular category of patient in response to colonization, etc.
        """
        if pthStatus == PthStatus.CLEAR:
            return 1.0
        else:
            return 0.0

    def setValue(self, key, val):
        """
        Setting values may be useful for changing phases in a scenario, for example. The
        values that can be set are treatment-specific; attempting to set an incorrect value
        is an error.
        
        The base class doesn't know how to set any values.
        """
        raise RuntimeError('Class %s does not know how to set the value %s'
                           % (type(self).__name__, key))


class TransferDestinationPolicy(Policy):
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

class PatientStats(object):
    def __init__(self):
        self.currentOccupancy = 0
        self.totalOccupancy = 0

    def addPatient(self):
        self.currentOccupancy += 1
        self.totalOccupancy += 1

    def remPatient(self):
        self.currentOccupancy -=1

class FacilityRegistry(object):
    def __init__(self,abbrev):
        ### there is a dict that will keep a list of patients on the registry of known something or other at this facility
        self.registryDict = {}
        self.facilityAbbrev = abbrev
        
    def startRegistry(self,key):
        ### Starting the registry with the same key will provide a no op
        if not self.isRegistryStarted(key):
            self.registryDict[key] = []
    
    def isRegistryStarted(self,key):
        return key in self.registryDict.keys()
    
    def registerPatient(self,key,patientName):
        if key not in self.registryDict.keys():
            raise RuntimeError("The Registry {1} at {0} has not been initialized".format(self.facilityAbbrev,key))
        
        if patientName not in self.registryDict[key]:
            self.registryDict[key].append(patientName)
    
    def transferRegistry(self, fromKey, toKey):
        if fromKey not in self.registryDict.keys():
            raise RuntimeError("In transferRegsitry: The Registry {1} at {0} has not been initialized".format(self.facilityAbbrev,fromKey))
        if toKey not in self.registryDict.keys():
            self.startRegistry(toKey)
        
        for name in self.registryDict[fromKey]:
            if name not in self.registryDict[toKey]:
                self.registryDict[toKey].append(name)
    
    def hasRegistry(self,key):
        return key in self.registryDict.keys()
    
    def isPatientInRegistry(self,key,patientName):
        if key not in self.registryDict.keys():
            raise RuntimeError("The Registry {1} at {0} has not been initialized".format(self.facilityAbbrev,key))
        return patientName in self.registryDict[key]
    
    def transferFromOtherRegistry(self,otherReg,fromKey,toKey):
        if fromKey not in otherReg.registryDict.keys():
            raise RuntimeError("In transferFromOtherRegsitry: The Registry {1} at {0} has not been initialized".format(self.facilityAbbrev,fromKey))
        if toKey not in self.registryDict.keys():
            self.startRegistry(toKey)
            
        for name in otherReg.registryDict[fromKey]:
            if name not in self.registryDict[toKey]:
                self.registryDict[toKey].append(name)
    
         
    def __str__(self):
        returnString = "<{0} Registry>\n".format(self.facilityAbbrev)
        for key,reg in self.registryDict.items():
            returnString += '   Key: {0}\n'.format(key)
            for pat in reg:
                returnString += '    - {0}\n'.format(pat)
        
        return returnString
    
        
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
        self.patientStats = PatientStats()
        self.registry = FacilityRegistry(self.abbrev)
        self.registry.startRegistry('knownCRECarrier')
        #self.registry.startRegistry('bundleCatches')
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
                    treatmentPolicyClass = pC
                if issubclass(pC, DiagnosticPolicy):
                    diagnosticPolicyClass = pC
        if not treatmentPolicyClasses:
            treatmentPolicyClasses = [TreatmentPolicy]
        self.transferDestinationPolicy = transferDestinationPolicyClass(patch, self.categoryNameMapper)
        self.treatmentPolicies = [treatmentPolicyClass(patch, self.categoryNameMapper)
                                  for treatmentPolicyClass in treatmentPolicyClasses]
        self.diagnosticPolicy = diagnosticPolicyClass(patch, self.categoryNameMapper)

    def __str__(self):
        return '<%s>' % self.name

    def setNoteHolder(self, noteHolder):
        self.noteHolder = noteHolder
        self.noteHolder.addNote({'name': self.name})

    def getNoteHolder(self):
        return self.noteHolder
    
    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        pass

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
            #print "myPayLoad = {0}".format(myPayload)
            #print "innerPayload = {0}".format(innerPayload)
            timeNow = super(Facility, self).handleIncomingMsg(msgType, innerPayload, timeNow)
            patientID, tier, isFrail, patientName = myPayload
            if patientID in self.patientDataDict:
                logger.debug('Patient %s has returned to %s' % (patientID, self.name))
                patientRec = pickle.loads(self.patientDataDict[patientID])
                if patientRec.departureDate is None:  # ie patient is being readmitted without previously leaving
                    self.patientStats.remPatient()
                patientRec.prevVisits += 1
                patientRec.arrivalDate = timeNow
                patientRec.departureDate = None
                self.patientDataDict[patientID] = pickle.dumps(patientRec)
            else:
#                patientRec = Facility.PatientRecord(patientID, timeNow, isFrail)
                patientRec = PatientRecord(patientID, timeNow, isFrail)
                self.patientDataDict[patientID] = pickle.dumps(patientRec,2)
                
            self.patientStats.addPatient()
            if timeNow != 0:  # Exclude initial populations
                nh = self.getNoteHolder()
                if nh:
                    nh.addNote({(CareTier.names[tier] + '_arrivals'): 1})
        elif issubclass(msgType, pyrheabase.DepartureMsg):
            myPayload, innerPayload = payload
            timeNow = super(Facility, self).handleIncomingMsg(msgType, innerPayload, timeNow)
            #print "innerPayload = {0}".format(innerPayload)
            patientID, tier, isFrail, patientName = myPayload
            if patientID not in self.patientDataDict:
                logger.error('%s has no record of patient %s' % (self.name, patientID))
            patientRec = pickle.loads(self.patientDataDict[patientID])
            patientRec.departureDate = timeNow
            self.patientDataDict[patientID] = pickle.dumps(patientRec,2)
            self.patientStats.remPatient()
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
        if senderAbbrev != self.abbrev: # Let's try not noting internal transfers
            self.getNoteHolder().addNote({(CareTier.names[tier] + '_found'): nSuccess,
                                          (CareTier.names[tier] + '_notfound'): nFail,
                                          bounceKey: nBounces,
                                          transferKey: nSuccess,
                                          fromToKey: nSuccess})

    def diagnose(self, patient, oldDiagnosis):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.
        """
        return self.diagnosticPolicy.diagnose(patient, oldDiagnosis, self)

    def prescribe(self, ward, patient, patientDiagnosis, patientTreatment):
        """
        This returns a tuple of either form (careTier, patientTreatment) or
        (careTier, patientTreatment, modifierList)
        """
        if patientDiagnosis.relocateFlag:
            modifierList = [pyrheabase.TierUpdateModFlag.FORCE_MOVE]
        else:
            modifierList = []
        careTier = None
        for tP in self.treatmentPolicies:
            newTier, patientTreatment = tP.prescribe(ward, patient, patientDiagnosis, patientTreatment,
                                                     modifierList)
            if careTier and (newTier != careTier):
                raise RuntimeError(('Treatment policies at %s prescribe different careTiers'
                                    % self.abbrev)
                                   (' for %s %s' % (patientDiagnosis, patientTreatment)))
            careTier = newTier
        return (careTier, patientTreatment, modifierList)

    def diagnosisFromCareTier(self, ward, timeNow):
        """
        If I look at a patient under a given care tier, what would I expect their diagnostic class
        to be?  This is used for patient initialization purposes.
        """
        return self.diagnosticPolicy.initializePatientDiagnosis(ward, timeNow)

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
        self._diagnosis = self.ward.fac.diagnosisFromCareTier(self.ward, timeNow)
        self._status = PatientStatus(PatientOverallHealth.HEALTHY,
                                     self._diagnosis.diagClassA, 0,
                                     self._diagnosis.pthStatus, 0,
                                     False, True, True, None)
        newTier, self._treatment = self.ward.fac.prescribe(self.ward,  # @UnusedVariable
                                                           self, self._diagnosis,
                                                           TREATMENT_DEFAULT)[0:2]

        self.lastUpdateTime = timeNow
        abbrev = self.ward.fac.abbrev
        #self.prevFac = None
        #self.cpReason = None
        self.id = (abbrev, PatientAgent.idCounters[abbrev])
        PatientAgent.idCounters[abbrev] += 1
        self.logger = logging.getLogger(__name__ + '.PatientAgent')

    def printSummary(self):
        print '%s as of %s' % (self.name, self.lastUpdateTime)
        print '  status: %s ' % str(self._status)
        print '  diagnosis: %s' % str(self._diagnosis)
        print '  treatment: %s' % self._treatment

    def getPthStatus(self):
        """Accessor for private status"""
        return self._status.pthStatus

    def getPthDiagnosis(self):
        """Accessor for private diagnosis"""
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
    
    def updateDiseaseState(self, treatment, facility, timeNow):
        """This should embody healing, community-acquired infection, etc."""
        dT = timeNow - self.lastUpdateTime
        previousStatus = self._status
        if dT > 0:  # moving from one ward to another can trigger two updates the same day
            treeL = []
            try:
                for src in [self.ward.fac, self.ward.iA]:
                    tree = src.getStatusChangeTree(self._status, self.ward, self._treatment,
                                                   self.lastUpdateTime, timeNow)
                    treeL.append(tree)
            except Exception, e:
                self.logger.critical('Got exception %s on patient %s (id %s) for tree %s'
                                     % (str(e), self.name, self.id, str(src)))
                raise
            try:
                treeL = self.ward.iA.filterStatusChangeTrees(treeL, self._status, self.ward, self._treatment,
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
                and self._status.pthStatus == PthStatus.COLONIZED):
                #print "New Infection at {0}".format(self.ward.fac.abbrev)
                self.ward.miscCounters['newColonizationsSinceLastChecked'] += 1
            if self.getTreatment('creBundle'):
                self.ward.miscCounters['creBundlesHandedOut'] += 1

    def updateEverything(self, timeNow):
        self.updateDiseaseState(self._treatment, self.ward.fac, timeNow)
        
        if self._status.diagClassA == DiagClassA.DEATH:
            return None, []
        self._diagnosis = self.ward.fac.diagnose(self, self._diagnosis)
        #print "prescibing in fac = {0}".format(self.ward.fac.abbrev)
        tpl = self.ward.fac.prescribe(self.ward, self, self._diagnosis, self._treatment)
        newTier, self._treatment = tpl[0:2]
        self._status = self._status._replace(justArrived=False)
        if len(tpl) == 3:
            modifierList = tpl[2]
        else:
            modifierList = []
        if self.debug:
            self.logger.debug('%s: status -> %s, diagnosis -> %s, treatment -> %s'
                              % (self.name, self._status, self._diagnosis,
                                 self._treatment))
            self.logger.debug('%s: modifiers are %s',
                              self.name,
                              [pyrheabase.TierUpdateModFlag.names[flg] for flg in modifierList])
    
        self.prevFac = self.ward.fac
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
            self.prevWardArr = self.locAddr
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
