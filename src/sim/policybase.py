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
from typebase import CareTier, PatientDiagnosis, PatientOverallHealth, DiagClassA
from pathogenbase import PthStatus, defaultPthStatus

from phacsl.utils.classutils import ClassIsInstanceMeta

logger = logging.getLogger(__name__)


class Policy(object):
    
    __metaclass__ = ClassIsInstanceMeta

    def __init__(self, patch, categoryNameMapper):
        self.patch = patch
        self.categoryNameMapper = categoryNameMapper


class DiagnosticPolicy(Policy):
    def handlePatientArrival(self, ward, patient, transferInfoDict, timeNow):
        """
        This is called on patients when they arrive at a ward.
        """
        logger.debug('%s arrives %s %s', '%s_%s'%patient.id, ward._name, timeNow)

    def diagnose(self, ward, patientId, patientStatus, oldDiagnosis, timeNow=None):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.
        """
        return PatientDiagnosis(patientStatus.overall,
                                patientStatus.diagClassA,
                                patientStatus.startDateA,
                                PthStatus.CLEAR,
                                patientStatus.relocateFlag)

    def initializePatientDiagnosis(self, careTier, overallHealth, timeNow):
        if careTier == CareTier.HOME:
            return PatientDiagnosis(overallHealth,
                                    DiagClassA.WELL, timeNow, defaultPthStatus, False)
        elif careTier == CareTier.NURSING:
            if overallHealth == PatientOverallHealth.FRAIL:
                return PatientDiagnosis(overallHealth, DiagClassA.WELL, timeNow,
                                        defaultPthStatus, False)
            else:
                # the patient must be here for rehab
                return PatientDiagnosis(overallHealth, DiagClassA.NEEDSREHAB, timeNow,
                                        defaultPthStatus, False)
        elif careTier == CareTier.LTAC:
            return PatientDiagnosis(overallHealth,
                                    DiagClassA.NEEDSLTAC, timeNow, defaultPthStatus, False)
        elif careTier == CareTier.HOSP:
            return PatientDiagnosis(overallHealth,
                                    DiagClassA.SICK, timeNow, defaultPthStatus, False)
        elif careTier == CareTier.ICU:
            return PatientDiagnosis(overallHealth,
                                    DiagClassA.VERYSICK, timeNow, defaultPthStatus, False)
        elif careTier == CareTier.VENT:
            return PatientDiagnosis(overallHealth,
                                    DiagClassA.NEEDSVENT, timeNow, defaultPthStatus, False)
        elif careTier == CareTier.SKILNRS:
            return PatientDiagnosis(overallHealth,
                                    DiagClassA.NEEDSSKILNRS, timeNow, defaultPthStatus, False)
        else:
            raise RuntimeError('Unknown care tier %s' % careTier)

    def sendPatientTransferInfo(self, facility, patient, transferInfoDict):
        """
        The information in the TransferInfo dictionary travels with the patient
        (actually with the BedRequest) from one facility to the next.  It may or may
        not contain useful information, for example whether or not the patient is contagious.
        This method provides an opportunity for the TreatmentPolicy to add info to
        the transferInfoDict.
        """
        return transferInfoDict

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
        if (ward.tier == CareTier.NURSING
                and patient.getDiagnosis().overall != PatientOverallHealth.FRAIL
                and patient.getDiagnosis().diagClassA == DiagClassA.NEEDSREHAB):
            patient.setTreatment(rehab=True)

    def handlePatientArrival(self, ward, patient, transferInfoDict, timeNow):
        """
        This is called on patients when they arrive at a ward.
        """
        raise RuntimeError('Base TreatmentPolicy was called for %s' % ward._name)

    def handlePatientDeparture(self, ward, patient, timeNow):
        """
        This is called on patients when they depart from a ward.
        """
        raise RuntimeError('Base TreatmentPolicy was called for %s' % ward._name)

    def sendPatientTransferInfo(self, facility, patient, transferInfoDict):
        """
        The information in the TransferInfo dictionary travels with the patient
        (actually with the BedRequest) from one facility to the next.  It may or may
        not contain useful information, for example whether or not the patient is contagious.
        This method provides an opportunity for the TreatmentPolicy to add info to
        the transferInfoDict.
        """
        return transferInfoDict

    def prescribe(self, ward, patientId, patientDiagnosis, patientTreatment, modifierList,
                  timeNow=None):
        """
        This returns a tuple of form (careTier, patientTreatment).
        modifierList is for functional modifiers, like pyrheabase.TierUpdateModFlag.FORCE_MOVE,
        and is not generally relevant to the decisions made by this method.
        """
        if patientDiagnosis.diagClassA == DiagClassA.WELL:
            if (patientDiagnosis.overall == PatientOverallHealth.HEALTHY
                or patientDiagnosis.overall == PatientOverallHealth.UNHEALTHY):
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
    def getOrderedCandidateFacList(self, facility, patientAgent, oldTier, newTier, timeNow):
        raise RuntimeError('Base TransferDestinationPolicy was called for %s' % facility.name)

class ScenarioPolicy(object):
    def __init__(self, name, patch):
        self.name = name
        self.patch = patch

    def begin(self, callingAgent, timeNow):
        logger.info('The scenario %s is beginning', self.name)
