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
from generic_diagnostic import GenericDiagnosticPolicy
from registry import Registry

_validator = None
_constants_values = '$(CONSTANTS)/xdro_registry_scenario_constants.yaml'
_constants_schema = 'xdro_registry_scenario_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)

class XDRORegistryScenario(BaseScenarioPolicy):
    def __init__(self, name, patch):
        super(XDRORegistryScenario, self).__init__(name, patch)
        self.logThisString = _constants['stringToLogWhenStarting']
        self.facSet = frozenset(_constants['locationsImplementingScenario']['locAbbrevList'])
        self.newEffectiveness = _constants['enhancedDetectionFraction']['value']

    def begin(self, callingAgent, timeNow):
        logger.warn(self.logThisString)
        assert hasattr(self.patch, 'allFacilities'), ('patch %s has no list of facilities!'
                                                      % self.patch.name)
        # All we have to do is turn on the intervention in all locations implementing
        # the scenario.
        for fac in self.patch.allFacilities:
            if fac.abbrev in self.facSet:
                if fac.diagnosticPolicy.isinstance(GenericDiagnosticPolicy):
                    #print 'XDRO setting %s' % fac.abbrev
                    fac.diagnosticPolicy.setValue('useCentralRegistry', True)
                    fac.diagnosticPolicy.setValue('pathogenDiagnosticEffectivenessIncreasedAwareness',
                                                  self.newEffectiveness)




def getPolicyClasses():
    return [XDRORegistryScenario]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(pyrheautils.pathTranslate(_constants_values),
                                         _constants_schema)
