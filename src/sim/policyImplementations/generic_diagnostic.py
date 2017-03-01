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
from facilitybase import PatientDiagnosis
from facilitybase import DiagnosticPolicy as BaseDiagnosticPolicy
from pathogenbase import PthStatus

_validator = None
_constants_values = '$(MODELDIR)/constants/generic_diagnosis_constants.yaml'
_constants_schema = 'generic_diagnosis_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)

class GenericDiagnosticPolicy(BaseDiagnosticPolicy):
    def __init__(self,  patch, categoryNameMapper):
        super(GenericDiagnosticPolicy, self).__init__(patch, categoryNameMapper)
        self.effectiveness = _constants['pathogenDiagnosticEffectiveness']['value']
        
    def diagnose(self, patientStatus, oldDiagnosis):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.
        
        This version provides some awareness of pathogen status
        """
        
        if patientStatus.justArrived:
            if patientStatus.pthStatus == PthStatus.COLONIZED and (random() <= self.effectiveness):
                diagnosedPthStatus = PthStatus.COLONIZED
            else:
                diagnosedPthStatus = PthStatus.CLEAR
        else:
            diagnosedPthStatus = oldDiagnosis.getPthStatus()
            
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
