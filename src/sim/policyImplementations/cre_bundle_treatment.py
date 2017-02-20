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
import random

from phacsl.utils.collections.phacollections import SingletonMetaClass
import pyrheautils
from collections import deque
from facilitybase import CareTier
from facilitybase import TreatmentPolicy as BaseTreatmentPolicy

_validator = None
_constants_values = '$(POLICYDIR)/constants/cre_bundle_treatment_constants.yaml'
_constants_schema = 'cre_treatment_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)

class CREBCore(object):
    """This is where we put things that are best shared across all instances"""
    __metaclass__ = SingletonMetaClass
    
    def __init__(self):
        self.effectiveness = _constants['transmissibility']['value']
    

class CREBundleTreatmentPolicy(BaseTreatmentPolicy):
    def __init__(self,patch,categoryNameMapper):
        super(CREBundleTreatmentPolicy,self).__init__(patch,categoryNameMapper)
        
    def initializePatientTreatment(self, ward, patient):
        """
        This is called on patients at time zero, when they are first assigned to the
        ward in which they start the simulation.
        """
        try:
            patient.setTreatment(creBundle=True)
        except:
            msg = ('Cannot set CRE Bundle treatment for patient in {0}'.format(ward))
            logger.fatal(msg)
            raise RuntimeError(msg)
    
    def handlePatientArrival(self, ward, patient, timeNow):
        """
        This is called on patients when they arrive at a ward.
        """
        self.initializePatientTreatment(ward,patient)
        ward.missCounters['creBundlesHandedOut'] += 1

    def handlePatientDeparture(self, ward, patient, timeNow):
        """
        This is called on patients when they depart from a ward.
        """
        patient.setTreatment(creBundle=False)

    def getTransmissionMultiplier(self, careTier, **kwargs):
            """
            If the treatment elements in **kwargs have the given boolean values (e.g. rehab=True),
            return the scale factor by which the transmission coefficient tau is multiplied.
            """
            if 'creBundle' in kwargs and kwargs['creBundle']:
                return self.core.effectiveness
            else:
                return 1.0

def getPolicyClasses():
    return [CREBundleTreatmentPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(pyrheautils.pathTranslate(_constants_values),
                                         _constants_schema)
