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
import random
from collections import deque
from facilitybase import TransferDestinationPolicy as BaseTransferDestinationPolicy
from facilitybase import CareTier, tierToQueueMap

_validator = None
_constants_values = 'transferbycapacity_constants.yaml'
_constants_schema = 'transferbycapacity_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)


class CapacityCore(object):
    """This is where we put things that are best shared across all instances"""
    __metaclass__ = pyrheautils.SingletonMetaClass

    def __init__(self, patch):
        self.patch = patch
        self.tierAddrMap = None
        transferMatrixFilePath = _constants['capacityFilePath']
        logger.info('Importing the constant data file %s' % transferMatrixFilePath)
        rawTbl = pyrheautils.importConstants(transferMatrixFilePath,
                                             _constants['capacityFileSchema'])
        nmDict = CareTier.names
        self.tbl = {tier: [] for tier in nmDict.keys()}
        for rec in rawTbl:
            abbrev = rec['abbrev']
            ctg = rec['category']
            if ctg == 'HOSPITAL':
                self.tbl[CareTier.HOSP].append((rec['capacity'], abbrev))
                if 'icuCapacity' in rec:
                    self.tbl[CareTier.ICU].append((rec['icuCapacity'], abbrev))
            elif ctg == 'LTAC':
                self.tbl[CareTier.LTAC].append((rec['capacity'], abbrev))
            elif ctg == 'NURSINGHOME':
                self.tbl[CareTier.NURSING].append((rec['capacity'], abbrev))
            elif ctg == 'COMMUNITY':
                self.tbl[CareTier.HOME].append((rec['capacity'], abbrev))
            else:
                raise RuntimeError('Capacity entry for %s has unknown category %s'
                                   % (abbrev, ctg))
        for k in nmDict.keys():
            self.tbl[k].sort()
        self.totTbl = {tier: sum([tpl[0] for tpl in self.tbl[tier]])
                       for tier in nmDict.keys()}
        logger.info('Import complete.')

    def buildTierAddrMap(self):
        self.tierAddrMap = {}
        nmDict = CareTier.names
        logger.info('Building facility-to-address map')
        for tier in nmDict.keys():
            self.tierAddrMap[tier] = {}
            serviceTupleList = self.patch.serviceLookup(tierToQueueMap[tier].__name__)
            for info, addr in serviceTupleList:
                innerInfo, facAbbrev = info  # @UnusedVariable
                self.tierAddrMap[tier][facAbbrev] = addr
        logger.info('map complete')

    def getTierAddrMap(self, tier):
        if self.tierAddrMap is None:
            self.buildTierAddrMap()
        return self.tierAddrMap[tier]


class CapacityTransferDestinationPolicy(BaseTransferDestinationPolicy):
    def __init__(self, patch):
        super(CapacityTransferDestinationPolicy, self).__init__(patch)
        self.core = CapacityCore(patch)

#     def buildCacheEntry(self, oldFacility, newTier):
#         srcAbbrev = oldFacility.abbrev
#         queueClass = tierToQueueMap[newTier]
#         tplList = []
#         for info, addr in self.patch.serviceLookup(queueClass.__name__):
#             innerInfo, destAbbrev = info  # @UnusedVariable
#             travelTime = self.core.tbl[srcAbbrev][destAbbrev]['seconds']
#             tplList.append((travelTime, info, addr))
#         tplList.sort()
#         return [c for a, b, c in tplList]  # @UnusedVariable

#     def getOrderedCandidateFacList(self, oldFacility, oldTier, newTier, timeNow):
#         pairList = self.core.tbl[newTier][:]
#         tot = self.core.totTbl[newTier]
# #         print 'newTier: %s' % CareTier.names[newTier]
#         facList = []
#         while pairList:
#             capSum = 0.0
#             idx = 0
#             lim = random.random() * tot
# #             print lim
#             while True:
# #                 print 'sum = %s idx = %s tpl = %s' % (capSum, idx, str(pairList[idx]))
#                 capSum += pairList[idx][0]
#                 if capSum >= lim:
#                     capacity, abbrev = pairList.pop(idx)  # @UnusedVariable
#                     facList.append(abbrev)
#                     tot -= capacity
# #                     print 'breaking; facList = %s' % facList
#                     break
#                 else:
#                     idx += 1
# #         print 'done! %s' % facList
#         tAM = self.core.getTierAddrMap(newTier)
#         return [tAM[facAbbrev] for facAbbrev in facList]

    def getOrderedCandidateFacList(self, oldFacility, oldTier, newTier, timeNow):
        pairList = deque(self.core.tbl[newTier])
        tot = self.core.totTbl[newTier]
#         print 'newTier: %s' % CareTier.names[newTier]
        facList = []
        facList_sav = [abbrev for capacity, abbrev in pairList]
        while pairList:
            capSum = 0.0
            lim = random.random() * tot
            newPL = deque()
#             print 'pairList: %s' % pairList
#             print 'lim: %s' % lim
            while True:
                capacity, abbrev = pairList.popleft()
#                 print 'sum = %s tpl = (%s, %s)' % (capSum, capacity, abbrev)
                capSum += capacity
                if capSum >= lim:
                    facList.append(abbrev)
                    tot -= capacity
                    newPL.extend(pairList)
                    pairList = newPL
#                     print 'breaking; facList = %s' % facList
                    break
                else:
                    newPL.append((capacity, abbrev))
#         print 'done! %s' % facList
#         facList = facList_sav
        tAM = self.core.getTierAddrMap(newTier)
        return [tAM[facAbbrev] for facAbbrev in facList]

#     def getOrderedCandidateFacList(self, oldFacility, oldTier, newTier, timeNow):
#         return (super(CapacityTransferDestinationPolicy, self)
#                 .getOrderedCandidateFacList(oldFacility, oldTier, newTier, timeNow))


def getPolicyClasses():
    return [CapacityTransferDestinationPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(os.path.join(os.path.dirname(__file__),
                                                      _constants_values),
                                         _constants_schema)
