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

_rhea_svn_id_ = "$Id$"

import logging

import os.path
import pyrheautils
from facilitybase import TransferDestinationPolicy as BaseTransferDestinationPolicy
from facilitybase import CareTier, tierToQueueMap

_validator = None
_constants_values = 'transferdestination_constants.yaml'
_constants_schema = 'transferdestination_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)


class MinTravelTimeCore(object):
    """This is where we put things that are best shared across all instances"""
    __metaclass__ = pyrheautils.SingletonMetaClass

    def __init__(self):
        transferMatrixFilePath = _constants['transitMatrixFilePath']
        logger.info('Importing the constant data file %s' % transferMatrixFilePath)
        rawTbl = pyrheautils.importConstants(transferMatrixFilePath,
                                             _constants['transitMatrixFileSchema'])
        self.tbl = {}
        for src, d in rawTbl.items():
            self.tbl[src] = {}
            for category, dd in d.items():  # @UnusedVariable
                for dest, entry in dd.items():
                    self.tbl[src][dest] = entry
        logger.info('Import complete.')


class MinTravelTimeTransferDestinationPolicy(BaseTransferDestinationPolicy):
    def __init__(self, patch):
        super(MinTravelTimeTransferDestinationPolicy, self).__init__(patch)
        self.core = MinTravelTimeCore()
        self.cache = {}

    def buildCacheEntry(self, oldFacility, newTier):
        srcAbbrev = oldFacility.abbrev
        queueClass = tierToQueueMap[newTier]
        tplList = []
        for info, addr in self.patch.serviceLookup(queueClass.__name__):
            innerInfo, destAbbrev = info  # @UnusedVariable
            travelTime = self.core.tbl[srcAbbrev][destAbbrev]['seconds']
            tplList.append((travelTime, info, addr))
        tplList.sort()
        return [c for a, b, c in tplList]  # @UnusedVariable

    def getOrderedCandidateFacList(self, oldFacility, oldTier, newTier, timeNow):
        if oldTier == CareTier.HOME or newTier == CareTier.HOME:
            # We have no travel time info for Community, so fall back
            return (super(MinTravelTimeTransferDestinationPolicy, self)
                    .getOrderedCandidateFacList(oldFacility, oldTier, newTier, timeNow))
        else:
            if newTier not in self.cache:
                self.cache[newTier] = self.buildCacheEntry(oldFacility, newTier)
#             print 'clause 3 %s: %s' % (oldFacility.abbrev, self.cache[newTier][:3])
            return self.cache[newTier][:]


def getPolicyClasses():
    return [MinTravelTimeTransferDestinationPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(os.path.join(os.path.dirname(__file__),
                                                      _constants_values),
                                         _constants_schema)
