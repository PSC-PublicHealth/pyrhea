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

import logging
from random import random
import types

import pyrheautils
from phacsl.utils.collections.phacollections import SingletonMetaClass
from facilitybase import PatientDiagnosis
from policybase import DiagnosticPolicy as BaseDiagnosticPolicy
from pathogenbase import PthStatus
from registry import Registry

_validator = None
_constants_values = '$(CONSTANTS)/$(PATHOGEN)/generic_diagnosis_constants.yaml'
_constants_schema = 'generic_diagnosis_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)

def parseConstantByFacilityCategory(fieldStr, key='category', innerKey='frac',
                                    constants=None):
    if constants is None:
        constants = _constants
    topD = {}
    for elt in constants[fieldStr]:
        cat = elt[key]
        topD[cat] = float(elt[innerKey]['value'])
    return topD

def _parseSameFacilityDiagnosisMemoryByCategory(fieldStr):
    return parseConstantByFacilityCategory(fieldStr)

def _parseCommunicateDiagnosisBetweenFacility(fieldStr, key='category'):
    return parseConstantByFacilityCategory(fieldStr, key=key)

def _parseRegistryAddCompliance(fieldStr):
    return parseConstantByFacilityCategory(fieldStr)

def _parseRegistrySearchCompliance(fieldStr):
    return parseConstantByFacilityCategory(fieldStr)

class GDPCore(object):
    """This is where we put things that are best shared across all instances"""
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        self.sameFacilityDiagnosisMemory = \
            _parseSameFacilityDiagnosisMemoryByCategory('sameFacilityDiagnosisMemory')
        self.rcvDiagnosisBetweenFacility = \
            _parseCommunicateDiagnosisBetweenFacility('receiveDiagnosisBetweenFacility',
                                                      key='categoryTo')
        self.sendDiagnosisBetweenFacility = \
            _parseCommunicateDiagnosisBetweenFacility('sendDiagnosisBetweenFacility',
                                                      key='categoryFrom')
        self.registryAddCompliance = \
            _parseRegistryAddCompliance('registryAddCompliance')
        self.registrySearchCompliance = \
            _parseRegistrySearchCompliance('registrySearchCompliance')


class GenericDiagnosticPolicy(BaseDiagnosticPolicy):
    def __init__(self, facility, patch, categoryNameMapper):
        #super(GenericDiagnosticPolicy, self).__init__(patch, categoryNameMapper)
        BaseDiagnosticPolicy.__init__(self, facility, patch, categoryNameMapper)
        self.effectiveness = _constants['pathogenDiagnosticEffectiveness']['value']
        self.falsePosRate = _constants['pathogenDiagnosticFalsePositiveRate']['value']
        self.core = GDPCore()
        self.increasedEffectivness = -1.0
        self.increasedFalsePosRate = -1.0
        self.useCentralRegistry = False

    def handlePatientArrival(self, ward, patient, transferInfoDict, timeNow):
        """
        This is called on patients when they arrive at a ward.
        """
        BaseDiagnosticPolicy.handlePatientArrival(self, ward, patient, transferInfoDict, timeNow)
        # Check for any info delivered with the transfer
        if 'carriesPth' in transferInfoDict:
            # Transfer probability was checked on the sending end
            rcvFacProb = self.core.rcvDiagnosisBetweenFacility[ward.fac.category]
            if random() <= rcvFacProb:
                with ward.fac.getPatientRecord(patient.id, timeNow=timeNow) as pRec:
                    pRec.carriesPth = True

    def handlePatientDeparture(self, ward, patient, timeNow):
        """
        This is called on patients when they depart from a ward.
        """
        BaseDiagnosticPolicy.handlePatientDeparture(self, ward, patient, timeNow)
        # Apparently there is a fair chance the patient record gets lost between visits
        if random() > self.core.sameFacilityDiagnosisMemory[ward.fac.category]:
            ward.fac.forgetPatientRecord(patient.id)

    def diagnose(self, ward, patientId, patientStatus, oldDiagnosis, timeNow=None):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.

        This version provides some awareness of pathogen status
        """
        if patientStatus.justArrived:
            with ward.fac.getPatientRecord(patientId, timeNow=timeNow) as pRec:

                if pRec.carriesPth:
                    diagnosedPthStatus = PthStatus.COLONIZED
                    if 'cpReason' not in pRec.noteD:
                        pRec.noteD['cpReason'] = 'passive'
                elif patientStatus.pthStatus == PthStatus.COLONIZED:
                    randVal = random()  # re-use this to get proper passive/xdro split
                    if randVal <= self.effectiveness:
                        diagnosedPthStatus = PthStatus.COLONIZED
                        pRec.noteD['cpReason'] = 'passive'
                    elif (self.useCentralRegistry and
                          (randVal <= self.increasedEffectiveness or
                           (random() <= self.core.registrySearchCompliance[ward.fac.category] and
                            Registry.getPatientStatus(str(ward.iA), patientId)))):
                        diagnosedPthStatus = PthStatus.COLONIZED
                        pRec.noteD['cpReason'] = 'xdro'
                    else:
                        diagnosedPthStatus = PthStatus.CLEAR  # Missed the diagnosis
                        pRec.noteD['cpReason'] = None
                else:
                    randVal = random()  # re-use this to get proper passive/xdro split
                    if randVal <= self.falsePosRate:
                        diagnosedPthStatus = PthStatus.COLONIZED
                        pRec.noteD['cpReason'] = 'passive'
                    elif self.useCentralRegistry and randVal <= self.increasedFalsePosRate:
                        diagnosedPthStatus = PthStatus.COLONIZED
                        pRec.noteD['cpReason'] = 'xdro'
                    else:
                        diagnosedPthStatus = PthStatus.CLEAR

                if diagnosedPthStatus in (PthStatus.COLONIZED, PthStatus.CHRONIC,
                                          PthStatus.INFECTED):
                    pRec.carriesPth = True  # though we might forget about this later
        else:
            diagnosedPthStatus = oldDiagnosis.pthStatus

        return PatientDiagnosis(patientStatus.overall,
                                patientStatus.diagClassA,
                                patientStatus.startDateA,
                                diagnosedPthStatus,
                                patientStatus.relocateFlag)

    def sendPatientTransferInfo(self, facility, patient, transferInfoDict):
        #transferInfoDict.update({'note': 'HUGE Hello from %s' % facility.name})
        # Maybe we remember to send known-carrier status with the patient
        pRec = facility.getPatientRecord(patient.id)
        if pRec.carriesPth:
            sendFacProb = self.core.sendDiagnosisBetweenFacility[facility.category]
            if random() <= sendFacProb:
                transferInfoDict['carriesPth'] = True

            if self.useCentralRegistry:
                if random() <= self.core.registryAddCompliance[facility.category]:
                    #print('here we are %s' % facility.category)
                    Registry.registerPatientStatus(patient.id,
                                                   str(patient.ward.iA),
                                                   patient._diagnosis,
                                                   facility.manager.patch)
        return transferInfoDict

    def setValue(self, key, val):
        """
        Setting values may be useful for changing phases in a scenario, for example. The
        values that can be set are treatment-specific; attempting to set an incorrect value
        is an error.

        This class supports setting pathogenDiagnosticEffectiveness, a fraction.
        """
        logger.info('%s setting %s to %s', type(self).__name__, key, val)
        if key == 'pathogenDiagnosticEffectiveness':
            if val:
                self.effectiveness = val
            else:
                # None means reset to default
                self.effectiveness = _constants['pathogenDiagnosticEffectiveness']['value']
            logger.info('effectiveness is now %s', self.effectiveness)
        elif key == 'pathogenDiagnosticEffectivenessIncreasedAwareness':
            if val:
                self.increasedEffectiveness = val
            else:
                # None means reset to initial setting
                self.increasedEffectiveness = -1.0
            logger.info('increasedEffectiveness is now %s', self.increasedEffectiveness)
        elif key == 'pathogenDiagnosticEffectivenessIncreasedFalsePositiveRate':
            if val:
                self.increasedFalsePosRate = val
            else:
                # None means reset to initial setting
                self.increasedFalsePosRate = -1.0
            logger.info('increasedFalsePosRate is now %s', self.increasedFalsePosRate)
        elif key == 'useCentralRegistry':
            assert type(val) == types.BooleanType, 'useCentralRegistry val should be boolean'
            self.useCentralRegistry = val
        else:
            super(GenericDiagnosticPolicy, self).setValue(key, val)


def getPolicyClasses():
    return [GenericDiagnosticPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(pyrheautils.pathTranslate(_constants_values),
                                         _constants_schema)
