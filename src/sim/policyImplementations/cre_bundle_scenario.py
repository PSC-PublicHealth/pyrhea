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
from policybase import ScenarioPolicy as BaseScenarioPolicy
from cre_bundle_treatment import CREBundleTreatmentPolicy
from cre_bundle_diagnostic import CREBundleDiagnosticPolicy

_validator = None
_constants_values = '$(CONSTANTS)/cre_bundle_scenario_constants.yaml'
_constants_schema = 'scenario_constants_schema.yaml'
_constants = None

LOGGER = logging.getLogger(__name__)


class CREBundleScenario(BaseScenarioPolicy):
    def __init__(self, name, patch):
        super(CREBundleScenario, self).__init__(name, patch)
        self.logThisString = _constants['stringToLogWhenStarting']
        self.locationImplementationInformation = _constants['locationsImplementingScenario']['facilities']
        self.evtList = []
        for elt in self.locationImplementationInformation:
            self.evtList.append((elt['times']['startDate'], elt['abbrev'], 'START'))
            self.evtList.append((elt['times']['endDate'], elt['abbrev'], 'END'))
        self.evtList.sort()


    def begin(self, callingAgent, timeNow):
        LOGGER.warn(self.logThisString)
        assert hasattr(self.patch, 'allFacilities'), ('patch %s has no list of facilities!'
                                                      % self.patch.name)
        baseTime = timeNow

        for when, abbrev, action in self.evtList:
            if baseTime + when != timeNow:
                assert(timeNow < baseTime + when), ('It is too late to %s intervention at %s'
                                                    % (action, abbrev))
                timeNow = callingAgent.sleep((baseTime + when) - timeNow)
            #print "{0}: {1} {2}".format(abbrev, when, action)

            for fac in self.patch.allFacilities:
                if fac.abbrev == abbrev:
                    fac.flushCaches()
                    for ward in fac.getWards():
                        ward.iA.flushCaches()
                    if action == 'START':
                        if type(fac.diagnosticPolicy).__name__ == CREBundleDiagnosticPolicy.__name__:
                            fac.diagnosticPolicy.setValue('active', True)
                        else:
                            raise RuntimeError('%s does not have a CREBundleDiagnosticPolicy' % abbrev)
                        for tP in fac.treatmentPolicies:
                            #if isinstance(tP, CREBundleTreatmentPolicy):
                            if type(tP).__name__ == CREBundleTreatmentPolicy.__name__:
                                tP.setValue('active', True)
                                for ward in fac.getWards():
                                    for patient in ward.getPatientList():
                                        tP.initializePatientTreatment(ward,patient)
                                LOGGER.info('Activated CREBundleScenario at %s' % abbrev)
                                break
                        else:
                            raise RuntimeError('%s does not have a CREBundleTreatmentPolicy' % abbrev)
                    elif action == 'END':
                        if type(fac.diagnosticPolicy).__name__ == CREBundleDiagnosticPolicy.__name__:
                            fac.diagnosticPolicy.setValue('active', False)
                        for tP in fac.treatmentPolicies:
                            if type(tP).__name__ == CREBundleTreatmentPolicy.__name__:
                                tP.setValue('active', False)
                                LOGGER.info('Deactivated CREBundleScenario at %s' % abbrev)
                                break
                    else:
                        raise RuntimeError('Nonsense action %s' % action)
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


