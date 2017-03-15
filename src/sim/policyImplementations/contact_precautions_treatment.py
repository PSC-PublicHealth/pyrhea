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

from phacsl.utils.collections.phacollections import SingletonMetaClass
import pyrheautils
from facilitybase import CareTier
from facilitybase import TreatmentPolicy as BaseTreatmentPolicy
from pathogenbase import PthStatus

_validator = None
_constants_values = '$(CONSTANTS)/contact_precautions_constants.yaml'
_constants_schema = 'contact_precautions_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)

def _parseFracByStatusByTierByCategory(fieldStr):
    topD = {}
    tierBackMap = {val: key for key, val in CareTier.names.items()}
    statusBackMap = {val: key for key, val in PthStatus.names.items()}
    for elt in _constants[fieldStr]:
        cat = elt['category']
        tiers = elt['tiers']
        if cat not in topD:
            topD[cat] = {}
        for tierElt in tiers:
            tier = tierBackMap[tierElt['tier']]
            statuses = tierElt['pathogenStatuses']
            if tier not in topD[cat]:
                topD[cat][tier] = {}
            for statusElt in statuses:
                status = statusBackMap[statusElt['status']]
                val = statusElt['frac']['value']
                assert status not in topD[cat][tier], ('Redundant %s for %s %s %s' %
                                                       (fieldStr, cat, CareTier.names[tier],
                                                        PthStatus.names[status]))
                topD[cat][tier][status] = val
    return topD


class CPTPCore(object):
    """This is where we put things that are best shared across all instances"""
    __metaclass__ = SingletonMetaClass
    
    def __init__(self):
        self.baseFracTbl = _parseFracByStatusByTierByCategory('baseFractionUnderContactPrecautions')
        self.effectiveness = _constants['transmissibilityMultiplier']['value']


class ContactPrecautionsTreatmentPolicy(BaseTreatmentPolicy):
    """This policy implements contact precautions"""
    
    """
    If the presence of this treatment corresponds to a flag in TreatmentProtocol,
    the name of that flag is the treatmentKey.
    """
    treatmentKey = 'contactPrecautions'

    def __init__(self, patch, categoryNameMapper):
        super(ContactPrecautionsTreatmentPolicy, self).__init__(patch, categoryNameMapper)
        self.core = CPTPCore()

    def initializePatientTreatment(self, ward, patient):
        """
        This is called on patients at time zero, when they are first assigned to the
        ward in which they start the simulation.
        """
        pthStatus = patient.getPthStatus()
        cat = ward.fac.category
        tier = ward.tier
        try:
            frac = self.core.baseFracTbl[ward.fac.category][ward.tier][pthStatus]
            patient.setTreatment(contactPrecautions=(random() <= frac))
        except KeyError:
            if tier in CareTier.names:
                tier = CareTier.names[tier]
            if pthStatus in PthStatus.names:
                pthStatus = PthStatus.names[pthStatus]
            msg = ('No baseFractionUnderContactPrecautions entry for %s %s %s'
                   % (cat, tier, pthStatus))
            logger.fatal(msg)
            raise RuntimeError(msg)

    def handlePatientArrival(self, ward, patient, timeNow):
        """
        This is called on patients when they arrive at a ward.
        """
        self.initializePatientTreatment(ward, patient)

    def handlePatientDeparture(self, ward, patient, timeNow):
        """
        This is called on patients when they depart from a ward.
        """
        patient.setTreatment(contactPrecautions=False) # Forget any contact precautions
        
    def getTransmissionFromMultiplier(self, careTier, **kwargs):
        """
        If the treatment elements in **kwargs have the given boolean values (e.g. rehab=True),
        return the scale factor by which the transmission coefficient tau is multiplied when
        the patient with this treatment is the source of the transmission.
        """
        if 'contactPrecautions' in kwargs and kwargs['contactPrecautions']:
            return self.core.effectiveness
        else:
            return 1.0

    def getTransmissionToMultiplier(self, careTier, **kwargs):
        """
        If the treatment elements in **kwargs have the given boolean values (e.g. rehab=True),
        return the scale factor by which the transmission coefficient tau is multiplied when
        the patient with this treatment is the recipient of the transmission.
        """
        if 'contactPrecautions' in kwargs and kwargs['contactPrecautions']:
            return self.core.effectiveness
        else:
            return 1.0

    def prescribe(self, ward, patientDiagnosis, patientTreatment, modifierList):
        """
        This returns a tuple of form (careTier, patientTreatment).
        modifierList is for functional modifiers, like pyrheabase.TierUpdateModFlag.FORCE_MOVE,
        and is not generally relevant to the decisions made by this method.
        """
        newTier, newTreatment = \
            super(ContactPrecautionsTreatmentPolicy, self).prescribe(ward,
                                                                     patientDiagnosis,
                                                                     patientTreatment,
                                                                     modifierList)
        if patientDiagnosis.pthStatus not in (PthStatus.CLEAR, PthStatus.RECOVERED):
            newTreatment = newTreatment._replace(contactPrecautions=True)
            
        # Apparently no one stays on CP for more than 10 days in Nursing care tier
        if ward.tier == CareTier.NURSING:
            if ward.patch.loop.sequencer.getTimeNow() - patientDiagnosis.startDateA > 10:
                newTreatment = newTreatment._replace(contactPrecautions=False)
            
        # for accounting
        if newTreatment.contactPrecautions:
            #print "patient is on contact precautions"
            ward.miscCounters['patientDaysOnCP'] += ward.checkInterval

        return newTier, newTreatment


def getPolicyClasses():
    return [ContactPrecautionsTreatmentPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(pyrheautils.pathTranslate(_constants_values),
                                         _constants_schema)
