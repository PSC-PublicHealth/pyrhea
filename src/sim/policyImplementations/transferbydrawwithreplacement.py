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

import os.path
from phacsl.utils.collections.phacollections import SingletonMetaClass
import pyrheautils
import random
from collections import deque
from policybase import TransferDestinationPolicy as BaseTransferDestinationPolicy
from facilitybase import CareTier, tierToQueueMap

_validator = None
_constants_values = 'transferbydrawwithreplacement_constants.yaml'
_constants_schema = 'transferbydrawwithreplacement_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)


class DWRCore(object):
    """This is where we put things that are best shared across all instances"""
    __metaclass__ = SingletonMetaClass

    def __init__(self, patch):
        self.patch = patch
        self.tierAddrMap = None
        self.tbl = None
        self.totTbl = None

    def _buildTierAddrMap(self):
        self.tierAddrMap = {}
        nmDict = CareTier.names
        logger.info('Building facility-to-address map')
        for tier in nmDict.keys():
            self.tierAddrMap[tier] = {}
            serviceTupleList = self.patch.serviceLookup(tierToQueueMap[tier].__name__)
            for info, addr in serviceTupleList:
                innerInfo, facAbbrev, facCoords = info  # @UnusedVariable
                self.tierAddrMap[tier][facAbbrev] = addr
        logger.info('map complete')

    def getTierAddrMap(self, tier):
        if self.tierAddrMap is None:
            self._buildTierAddrMap()
        return self.tierAddrMap[tier]
    
    def _buildWeightedLists(self):
        nmDict = CareTier.names
        tierFacSets = {tier: set(self.getTierAddrMap(tier).keys()) for tier in nmDict.keys()}
        pairsSeen = set()
        self.tbl = {}
        self.totTbl = {}
        for transferMatrixFilePath in _constants['transferFilePaths']:
            logger.info('Importing the weight data file %s', transferMatrixFilePath)
            rawTbl = pyrheautils.importConstants(transferMatrixFilePath,
                                                 _constants['transferFileSchema'])
            for srcName, rec in rawTbl.items():
                if srcName not in self.tbl:
                    self.tbl[srcName] = {}
                    self.totTbl[srcName] = {}
                for destName in rec.keys():
                    if (srcName, destName) in pairsSeen:
                        raise RuntimeError('Duplicate weight table entries for %s -> %s' %
                                           (srcName, destName))
                    else:
                        pairsSeen.add((srcName, destName))
                for tier in nmDict.keys():
                    if tier not in self.tbl[srcName]:
                        self.tbl[srcName][tier] = []
                        self.totTbl[srcName][tier] = 0.0
                    wtL = self.tbl[srcName][tier]
                    wtSum = self.totTbl[srcName][tier]
                    for destName, ct in rec.items():
                        if destName in tierFacSets[tier]:
                            wtL.append((ct, (destName, self.getTierAddrMap(tier)[destName])))
                            wtSum += ct
                    self.totTbl[srcName][tier] = wtSum
            logger.info('Import complete.')
        for srcName, subTbl in self.tbl.items():
            for wtL in subTbl.values():     
                wtL.sort(reverse=True)
        logger.info('Post-processing of weight data files is complete.')

    def getTierWeightedList(self, srcName, tier):
        """
        Returns a tuple containing:
        - a list of the form [(wt, (abbrev, addr)), ...] to be shuffled
        - a value wtSum equal to the sum of the weights in the list
        """
        if self.tbl is None:
            self._buildWeightedLists()
        return (self.tbl[srcName][tier][:], self.totTbl[srcName][tier])


class DrawWithReplacementTransferDestinationPolicy(BaseTransferDestinationPolicy):
    def __init__(self, patch, categoryNameMapper):
        super(DrawWithReplacementTransferDestinationPolicy, self).__init__(patch, categoryNameMapper)
        self.core = DWRCore(patch)

    def getOrderedCandidateFacList(self, oldFacility, oldTier, newTier, timeNow):
        pairList, tot = self.core.getTierWeightedList(oldFacility.abbrev, newTier)
#         print '%s %s -> %s: tot=%d' % (oldFacility.name, CareTier.names[oldTier], CareTier.names[newTier], tot)
#         print 'pairList: %s' % str([(a, b[0]) for a, b in pairList])
        pairList = deque(pairList)
#         print 'newTier: %s' % CareTier.names[newTier]
        facList = []
#        facList_sav = [destTpl for capacity, destTpl in pairList]
        try:
            while pairList:
                capSum = 0.0
                lim = random.random() * tot
                newPL = deque()
    #             print 'pairList: %s' % pairList
    #             print 'lim: %s' % lim
                while True:
                    if pairList:
                        capacity, destTpl = pairList.popleft()
    #                     print 'sum = %s tpl = (%s, %s)' % (capSum, capacity, destTpl)
                        capSum += capacity
                        if capSum >= lim:
                            facList.append(destTpl)
                            tot -= capacity
                            newPL.extend(pairList)
                            pairList = newPL
    #                         print 'breaking; facList = %s' % facList
                            break
                        else:
                            newPL.append((capacity, destTpl))
                    else:
                        break
        except IndexError, e:
            logger.error('Hit IndexError %s for %s %s -> %s at %s', e, oldFacility.abbrev,
                           oldTier, newTier, timeNow)
            logger.error('tot = %s, lim = %s, capSum = %s', tot, lim, capSum)
            logger.error('newPL is %s', str(newPL))
            pairList, tot = self.core.getTierWeightedList(oldFacility.abbrev, newTier)
            logger.error('initial tot = %s, pairList = %s', tot, pairList)
            raise
        return [b for a, b in facList]


def getPolicyClasses():
    return [DrawWithReplacementTransferDestinationPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(os.path.join(os.path.dirname(__file__),
                                                      _constants_values),
                                         _constants_schema)
