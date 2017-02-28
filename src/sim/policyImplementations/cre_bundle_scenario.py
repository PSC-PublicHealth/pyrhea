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
from pyrheabase import ScenarioPolicy as BaseScenarioPolicy
from cre_bundle_treatment import CREBundleTreatmentPolicy

_validator = None
_constants_values = '$(MODELDIR)/constants/cre_bundle_scenario_constants.yaml'
_constants_schema = 'scenario_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)

# These have to be in order of start-up time, because it's late and I'm tired.
# Format of tuples is (abbrev, startDate, endDate) where the dates start counting
# after burn-in.  For consistency with the calibration, we assume 11/28/2011 
# corresponds to post-burnin day number 150, and that date is the trigger date for
# this intervention.
# testLocations = [('THC_4058_L', 250, 830),
#                  ('THC_365_L', 320, 830),
#                  ('THC_6130_L', 383, 830),
#                  ('THC_2544_L', 452, 830)
#                  ]
testLocations = [('THC_4058_L', 4, 8),
                 ('ADVO_450_H', 4, 8),
                 ('ADVO_801_H', 4, 8),
                 ('ADVO_836_H', 4, 8),
                 ('THC_365_L', 6, 8),
                 ]

class CREBundleScenario(BaseScenarioPolicy):
    def __init__(self, name, patch):
        super(CREBundleScenario, self).__init__(name, patch)
        self.logThisString = _constants['stringToLogWhenStarting']
        
    def begin(self, callingAgent, timeNow):
        logger.warn(self.logThisString)
        assert hasattr(self.patch, 'allFacilities'), ('patch %s has no list of facilities!'
                                                      % self.patch.name)
        for abbrev, startDate, endDate in testLocations:
            if timeNow != startDate:
                assert(timeNow < startDate), 'It is too late to start intervention at %s' % abbrev
                timeNow = callingAgent.sleep(startDate - timeNow)
            for fac in self.patch.allFacilities:
                if fac.abbrev == abbrev:
                    fac.flushCaches()
                    for ward in fac.getWards():
                        ward.iA.flushCaches()
                    for tP in fac.treatmentPolicies:
                        if (type(tP).__name__ == CREBundleTreatmentPolicy.__name__):
                        #if isinstance(tP, CREBundleTreatmentPolicy):
                            tP.setValue('active', True)
                            print ('Activated CREBundleScenario at %s' % abbrev)
                            logger.info('Activated CREBundleScenario at %s' % abbrev)
                            break
                    else:
                        raise RuntimeError('%s does not have a CREBundleTreatmentPolicy' % abbrev)
                    break
            else:
                raise RuntimeError('Failed to find the facility %s' % abbrev)
        for abbrev, startDate, endDate in testLocations:
            if timeNow != endDate:
                assert(timeNow < endDate), 'It is too late to start intervention at %s' % abbrev
                timeNow = callingAgent.sleep(endDate - timeNow)
            for fac in self.patch.allFacilities:
                if fac.abbrev == abbrev:
                    fac.flushCaches()
                    for ward in fac.getWards():
                        ward.iA.flushCaches()
                    for tP in fac.treatmentPolicies:
                        if (type(tP).__name__ == CREBundleTreatmentPolicy.__name__):
                        #if isinstance(tP, CREBundleTreatmentPolicy):
                            tP.setValue('active', False)
                            print ('Deactivated CREBundleScenario at %s' % abbrev)
                            logger.info('Deactivated CREBundleScenario at %s' % abbrev)
                            break
                    else:
                        raise RuntimeError('%s does not have a CREBundleTreatmentPolicy' % abbrev)
                    break
            else:
                raise RuntimeError('Failed to find the facility %s' % abbrev)

def getPolicyClasses():
    return [CREBundleScenario]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(pyrheautils.pathTranslate(_constants_values),
                                         _constants_schema)
