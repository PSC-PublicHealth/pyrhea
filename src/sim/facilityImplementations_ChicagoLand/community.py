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
logger = logging.getLogger(__name__)

from phacsl.utils.collections.phacollections import enum
import genericCommunity

PatientCategory = enum('BASE', 'COLONIZED')

class CommunityWard(genericCommunity.CommunityWard):
    def classify(self, agent):
        """Return the PatientCategory appropriate for this agent"""
        if agent._status.pthStatus == PthStatus.COLONIZED:
            return PatientCategory.COLONIZED
        elif agent._status.pthStatus == PthStatus.CLEAR:
            return PatientCategory.BASE
        else:
            raise CommunityWard.FreezerError('%s has unexpected PthStatus %s' %
                                             (agent.name, PthStatus.names[agent._status.pthStatus]))

        
class Community(genericCommunity.Community):
    def setCDFs(self, losModel):
        baseRate = losModel['parms'][0]
        colonizedRate = _constants['communityColonizedReadmissionRateScale']['value'] * baseRate
        self.cachedCDFs = {PatientCategory.BASE: CachedCDFGenerator(expon(scale=1.0/baseRate)),
                           PatientCategory.COLONIZED: CachedCDFGenerator(expon(scale=1.0/colonizedRate))}
        

genericCommunity.importCommunity(__name__)
