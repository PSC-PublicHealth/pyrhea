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
from random import choice

import pyrheautils
from phacsl.utils.collections.phacollections import SingletonMetaClass
from typebase import CareTier
from facilitybase import PatientDiagnosis
from generic_diagnostic import GenericDiagnosticPolicy
from pathogenbase import PthStatus
import labwork

_validator = None
_constants_values = '$(CONSTANTS)/mrsa_control_treatment_constants.yaml'
_constants_schema = 'mrsa_control_treatment_constants_schema.yaml'
_constants = None

LOGGER = logging.getLogger(__name__)

TIERS_IMPLEMENTING = frozenset([('NURSINGHOME', CareTier.NURSING),
                                ('LTAC', CareTier.LTAC),
                                ('HOSPITAL', CareTier.ICU)])


def parseConstantByTier(fieldStr, key='tier', innerKey='frac', constants=None):
    backMap = {tierNm: tier for tier, tierNm in CareTier.names.items()}
    if constants is None:
        constants = _constants
    topD = {}
    for elt in constants[fieldStr]:
        tier = backMap[elt[key]]
        topD[tier] = elt[innerKey]['value']
    return topD


def _parseNumByTier(fieldStr):
    return parseConstantByTier(fieldStr, innerKey='num')


class SwabTestCore(object):
    """This is where we put things that are best shared across all instances"""
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        self.swabDelayDaysByTierMin = \
            _parseNumByTier('swabDelayDaysByTierMin')
        self.swabDelayDaysByTierMax = \
            _parseNumByTier('swabDelayDaysByTierMax')
        coveredTiersMin = set(self.swabDelayDaysByTierMin.keys())
        coveredTiersMax = set(self.swabDelayDaysByTierMin.keys())
        assert coveredTiersMin == coveredTiersMax, ('Max and min swab test delay data'
                                                    ' is not consistent')


class SwabTest(labwork.LabWork):
    def __init__(self, minDelayDays, maxDelayDays, debug=False):
        self.core = SwabTestCore()
        delay = choice(range(minDelayDays, maxDelayDays + 1))
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


class MRSAControlDiagnosticPolicy(GenericDiagnosticPolicy):
    def __init__(self,  facility, patch, categoryNameMapper):
        #super(MRSAControlDiagnosticPolicy, self).__init__(patch, categoryNameMapper)
        GenericDiagnosticPolicy.__init__(self, facility, patch, categoryNameMapper)
        self.swabCore = SwabTestCore()
        self.swabTestD = {}
        self.swabTestD = {tier: SwabTest(self.swabCore.swabDelayDaysByTierMin[tier],
                                         self.swabCore.swabDelayDaysByTierMax[tier])
                          for tier in self.swabCore.swabDelayDaysByTierMin}
        self.active = False

    def diagnose(self, ward, patientId, patientStatus, oldDiagnosis, timeNow=None):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.

        This version provides some awareness of pathogen status
        """

        # Layer our diagnosis on top of any done by the base class
        parentDiagnosis = super(MRSAControlDiagnosticPolicy, self).diagnose(ward,
                                                                            patientId,
                                                                            patientStatus,
                                                                            oldDiagnosis,
                                                                            timeNow=timeNow)

        with ward.fac.getPatientRecord(patientId, timeNow=timeNow) as pRec:
            diagnosedPthStatus = parentDiagnosis.pthStatus
            if patientStatus.justArrived:
                if (ward.tier == CareTier.ICU  # mandated by law.
                        or self.active):       # current scenario requires testing?
                    if ((not pRec.isContagious)
                        and parentDiagnosis.pthStatus in (PthStatus.CLEAR, PthStatus.RECOVERED)
                        and (ward.fac.implCategory, ward.tier) in TIERS_IMPLEMENTING):
                        # Only test patients not already known to need isolation. Trust
                        # whatever method marked pRec to have set appropriate counters.
                        #
                        # This patient gets tested
                        self.swabTestD[ward.tier].performLab(patientId, patientStatus, ward,
                                                             timeNow)
                        ward.miscCounters['swabsPerformed'] += 1

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
            super(MRSAControlDiagnosticPolicy, self).setValue(key, val)


def getPolicyClasses():
    return [MRSAControlDiagnosticPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(pyrheautils.pathTranslate(_constants_values),
                                         _constants_schema)
