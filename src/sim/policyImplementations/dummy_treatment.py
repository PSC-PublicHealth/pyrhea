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
from policybase import TreatmentPolicy as BaseTreatmentPolicy

_validator = None
_constants_values = '$(POLICYDIR)/transferatrandom_constants.yaml'
_constants_schema = 'transferatrandom_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)


class DummyTreatmentPolicy(BaseTreatmentPolicy):
    def initializePatientTreatment(self, ward, patient):
        """
        This is called on patients at time zero, when they are first assigned to the
        ward in which they start the simulation.
        """
        super(DummyTreatmentPolicy, self).initializePatientTreatment(ward, patient)

    def handlePatientArrival(self, ward, patient, transferInfoDict, timeNow):
        """
        This is called on patients when they arrive at a ward.
        """
        pass

    def handlePatientDeparture(self, ward, patient, timeNow):
        """
        This is called on patients when they depart from a ward.
        """
        pass

    def getTransmissionFromMultiplier(self, careTier, **kwargs):
        """
        If the treatment elements in **kwargs have the given boolean values (e.g. rehab=True),
        return the scale factor by which the transmission coefficient tau is multiplied when
        the patient with this treatment is the source of the transmission.
        """
        return 1.0

    def getTransmissionToMultiplier(self, careTier, **kwargs):
        """
        If the treatment elements in **kwargs have the given boolean values (e.g. rehab=True),
        return the scale factor by which the transmission coefficient tau is multiplied when
        the patient with this treatment is the recipient of the transmission.
        """
        return 1.0


def getPolicyClasses():
    return [DummyTreatmentPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(pyrheautils.pathTranslate(_constants_values),
                                         _constants_schema)
