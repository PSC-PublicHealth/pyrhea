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
from facilitybase import TransferDestinationPolicy as BaseTransferDestinationPolicy
from facilitybase import CareTier, tierToQueueMap

_validator = None
_constants_values = '$(POLICYDIR)/transferatrandom_constants.yaml'
_constants_schema = 'transferatrandom_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)


class RandomTransferDestinationPolicy(BaseTransferDestinationPolicy):
    def __init__(self, patch, categoryNameMapper):
        super(RandomTransferDestinationPolicy, self).__init__(patch, categoryNameMapper)

    def getOrderedCandidateFacList(self, facility, oldTier, newTier, timeNow):
        queueClass = tierToQueueMap[newTier]
        facAddrList = [tpl[1] for tpl in self.patch.serviceLookup(queueClass.__name__)]
        random.shuffle(facAddrList)
        return facAddrList



def getPolicyClasses():
    return [RandomTransferDestinationPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(pyrheautils.pathTranslate(_constants_values),
                                         _constants_schema)
