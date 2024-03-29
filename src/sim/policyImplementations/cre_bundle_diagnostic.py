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

import pyrheautils
from phacsl.utils.collections.phacollections import SingletonMetaClass
from facilitybase import PatientDiagnosis
from generic_diagnostic import GenericDiagnosticPolicy
from pathogenbase import PthStatus
import labwork
from generic_diagnostic import parseConstantByFacilityCategory
from cre_bundle_treatment import TIERS_IMPLEMENTING_BUNDLE

_validator = None
_constants_values = '$(CONSTANTS)/cre_bundle_treatment_constants.yaml'
_constants_schema = 'cre_bundle_treatment_constants_schema.yaml'
_constants = None

LOGGER = logging.getLogger(__name__)

class SwabTestCore(object):
    """This is where we put things that are best shared across all instances"""
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        self.swabDelayDaysByCategory = \
            parseConstantByFacilityCategory('swabDelayDaysByCategory', innerKey='num',
                                            constants=_constants)


class SwabTest(labwork.LabWork):
    def __init__(self, category, debug=False):
        self.core = SwabTestCore()
        delay = self.core.swabDelayDaysByCategory[category]
        super(SwabTest, self).__init__(sensitivity=_constants['swabDiagnosticSensitivity']['value'],
                                       specificity=_constants['swabDiagnosticSpecificity']['value'],
                                       delayDays=delay,
                                       debug=debug)

    def trueTestFun(self, patientStatus):
        return patientStatus.pthStatus in (PthStatus.COLONIZED, PthStatus.CHRONIC,
                                           PthStatus.INFECTED)

    @classmethod
    def posAction(cls, patientRecord):
        patientRecord.carriesPth = True
        patientRecord.noteD['cpReason'] = 'swab'
        return patientRecord

    @classmethod
    def negAction(cls, patientRecord):
        patientRecord.carriesPth = False
        return patientRecord


class CREBundleDiagnosticPolicy(GenericDiagnosticPolicy):
    def __init__(self,  facility, patch, categoryNameMapper):
        #super(CREBundleDiagnosticPolicy, self).__init__(patch, categoryNameMapper)
        GenericDiagnosticPolicy.__init__(self, facility, patch, categoryNameMapper)
        self.swabTest = SwabTest(facility.category)
        self.active = False

    def diagnose(self, ward, patientId, patientStatus, oldDiagnosis, timeNow=None):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.

        This version provides some awareness of pathogen status
        """

        # Layer our diagnosis on top of any done by the base class
        parentDiagnosis = super(CREBundleDiagnosticPolicy, self).diagnose(ward,
                                                                          patientId,
                                                                          patientStatus,
                                                                          oldDiagnosis,
                                                                          timeNow=timeNow)

        with ward.fac.getPatientRecord(patientId, timeNow=timeNow) as pRec:
            diagnosedPthStatus = parentDiagnosis.pthStatus
            if self.active and patientStatus.justArrived:
                if ((not pRec.isContagious)
                    and parentDiagnosis.pthStatus in (PthStatus.CLEAR, PthStatus.RECOVERED)
                    and (ward.fac.implCategory, ward.tier) in TIERS_IMPLEMENTING_BUNDLE):
                    # Only test patients not already known to need isolation. Trust
                    # whatever method marked pRec to have set appropriate counters.
                    #
                    # This patient gets tested
                    self.swabTest.performLab(patientId, patientStatus, ward, timeNow)
                    ward.miscCounters['creSwabsPerformed'] += 1

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

        This class can be activated and deactivated via the 'active' flag (True/False)
        """
        if key == 'active':
            self.active = val
        else:
            super(CREBundleDiagnosticPolicy, self).setValue(key, val)


def getPolicyClasses():
    return [CREBundleDiagnosticPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(pyrheautils.pathTranslate(_constants_values),
                                         _constants_schema)
