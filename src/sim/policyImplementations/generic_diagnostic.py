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

import pyrheautils
from phacsl.utils.collections.phacollections import SingletonMetaClass
from facilitybase import PatientDiagnosis
from facilitybase import DiagnosticPolicy as BaseDiagnosticPolicy
from pathogenbase import PthStatus

_validator = None
_constants_values = '$(CONSTANTS)/generic_diagnosis_constants.yaml'
_constants_schema = 'generic_diagnosis_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)

def _parseSameFacilityDiagnosisMemoryByCategory(fieldStr):
    topD = {}
    for elt in _constants[fieldStr]:
        cat = elt['category']
        topD[cat] = float(elt['frac']['value'])
    
    return topD

def _parseCommunicateDiagnosisBetweenFacility(fieldStr):
    topD = {}
    for elt in _constants[fieldStr]:
        cat = elt['categoryFrom']
        topD[cat] = float(elt['frac']['value'])
    return topD    

class GDPCore(object):
    """This is where we put things that are best shared across all instances"""
    __metaclass__ = SingletonMetaClass
    
    def __init__(self):
        self.sameFacilityDiagnosisMemory = _parseSameFacilityDiagnosisMemoryByCategory('sameFacilityDiagnosisMemory')
        self.communicateDiagnosisBetweenFacility = _parseCommunicateDiagnosisBetweenFacility('communitcateDiagnosisBetweenFacility')
        
class GenericDiagnosticPolicy(BaseDiagnosticPolicy):
    def __init__(self,  patch, categoryNameMapper):
        #super(GenericDiagnosticPolicy, self).__init__(patch, categoryNameMapper)
        BaseDiagnosticPolicy.__init__(self, patch, categoryNameMapper)
        self.effectiveness = _constants['pathogenDiagnosticEffectiveness']['value']
        self.falsePosRate = _constants['pathogenDiagnosticFalsePositiveRate']['value']
        self.core = GDPCore()
        
    def diagnose(self, patient, oldDiagnosis, facility = None):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.
        
        This version provides some awareness of pathogen status
        """
        patientStatus = patient._status
        if patientStatus.justArrived:
            
            ### first, if this patient is transferred, did their status come with them
            ### They would have to be on the previous facility's registry first
            
            if patient.prevFac is not None and \
                patient.prevFac.category != "COMMUNITY" and \
                patient.prevFac.registry.isPatientInRegistry('knownCRECarrier',patient.name):
                ## This patient was on the registry at the previous location, lets do a dice role depending on
                ## If the previous location is different than the current
                if facility is not None: 
                    if patient.prevFac == facility:
                        pass
                        ### Patient being readmitted to the same facility
                        #sameFacProb = self.core.sameFacilityDiagnosisMemory[facility.category]
                        #if random() <= sameFacProb:
                        #    facility.registry.registerPatient('knownCRECarrier',patient.name)     
                    else:
                        #print "Previous Facililty category = {0}".format(patient.prevFac.category)
                        commFacProb = self.core.communicateDiagnosisBetweenFacility[patient.prevFac.category]
                        if random() <= commFacProb:
                            facility.registry.registerPatient('knownCRECarrier',patient.name)
            
                
            if facility is not None and facility.registry.isPatientInRegistry('knownCRECarrier',patient.name):
                diagnosedPthStatus = PthStatus.COLONIZED
                
            elif patientStatus.pthStatus == PthStatus.COLONIZED:
                diagnosedPthStatus = (PthStatus.COLONIZED if (random() <= self.effectiveness)
                                      else PthStatus.CLEAR)
            else:
                diagnosedPthStatus = (PthStatus.COLONIZED if (random() <= self.falsePosRate)
                                      else PthStatus.CLEAR)
            
            ### if they are known to be colonized at the facility, then we need ot see if they go on the 
            ### registry
            if diagnosedPthStatus == PthStatus.COLONIZED:
                sameFacProb = self.core.sameFacilityDiagnosisMemory[facility.category]
                if random() <= sameFacProb:
                    facility.registry.registerPatient('knownCRECarrier',patient.name)          
                    
        else:
            diagnosedPthStatus = oldDiagnosis.pthStatus
            
        return PatientDiagnosis(patientStatus.overall,
                                patientStatus.diagClassA,
                                patientStatus.startDateA,
                                diagnosedPthStatus,
                                patientStatus.relocateFlag)

    def setValue(self, key, val):
        """
        Setting values may be useful for changing phases in a scenario, for example. The
        values that can be set are treatment-specific; attempting to set an incorrect value
        is an error.
        
        This class supports setting pathogenDiagnosticEffectiveness, a fraction.
        """
        if key == 'pathogenDiagnosticEffectiveness':
            self.effectiveness = val
            print 'effectiveness is now %s' % self.effectiveness
        else:
            super(GenericDiagnosticPolicy, self).setValue(key, val)


def getPolicyClasses():
    return [GenericDiagnosticPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(pyrheautils.pathTranslate(_constants_values),
                                         _constants_schema)
