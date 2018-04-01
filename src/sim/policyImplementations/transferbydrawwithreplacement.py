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
from collections import deque, defaultdict

from phacsl.utils.collections.phacollections import SingletonMetaClass
import pyrheautils
from policybase import TransferDestinationPolicy as BaseTransferDestinationPolicy
from facilitybase import CareTier, tierToQueueMap

_validator = None
_constants_values = '$(MODELDIR)/constants/transferbydrawwithreplacement_constants.yaml'
_constants_schema = 'transferbydrawwithreplacement_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)


def randomOrderByWt(pairList, tot, cull=None):
    """
    Given a list of the form [(weight, info), ...] in reverse sorted order
    and a value tot which is the sum of all the weights, return a list of
    info elements in weighted random order.  If cull is not None, exclude
    info[0]b==cull from the randomized list.
    """
    if cull is None:
        pairList = deque(pairList)
    else:
        tot -= sum([wt for wt, info in pairList if info[0] == cull])
        pairList = deque([(wt, info) for wt, info in pairList if info[0] != cull])
    shuffledList = []
#    shuffledList_sav = [info for weight, info in pairList]
    try:
        while pairList:
            wtSum = 0.0
            lim = random.random() * tot
            newPL = deque()
#             print 'pairList: %s' % pairList
#             print 'lim: %s' % lim
            while True:
                if pairList:
                    weight, info = pairList.popleft()
#                     print 'sum = %s tpl = (%s, %s)' % (wtSum, weight, info)
                    wtSum += weight
                    if wtSum >= lim:
                        shuffledList.append(info)
                        tot -= weight
                        newPL.extend(pairList)
                        pairList = newPL
#                         print 'breaking; shuffledList = %s' % shuffledList
                        break
                    else:
                        newPL.append((weight, info))
                else:
                    break
    except IndexError, e:
        logger.error('randomOrderByWt: %s: tot = %s, lim = %s, wtSum = %s', e, tot, lim, wtSum)
        logger.error('randomOrderByWt: newPL is %s', str(newPL))
        raise
    return shuffledList


class DWRCore(object):
    """This is where we put things that are best shared across all instances"""
    __metaclass__ = SingletonMetaClass

    def __init__(self, patch):
        self.patch = patch
        self.transferFilePaths = _constants['transferFilePaths']
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
        for transferMatrixFilePath in self.transferFilePaths:
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
                for tier in nmDict:
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

        # Let's try adding some options for rare facilities
        # - say if there are fewer than 3%
        totNumFacs = sum([len(fS) for tier, fS in tierFacSets.items()
                          if tier != CareTier.HOME])
        thresh = int(round(0.03 * totNumFacs))
        adjustThese = set()
        for tier, fS in tierFacSets.items():
            if len(fS) <= thresh:
                adjustThese.add(tier)
        for tier in adjustThese:
            fS = tierFacSets[tier]
            addrMap = self.getTierAddrMap(tier)
            for srcName in self.tbl:
                if tier in self.tbl[srcName]:
                    oldL = self.tbl[srcName][tier]
                    if oldL:
                        dct = defaultdict(lambda: 0)
                        for ct, (nm, addr) in oldL:
                            dct[nm] = ct
                        oldTot = sum(dct.values())
                        allNmL = [nm for nm in addrMap]
                        delta = (0.1 * oldTot)/len(allNmL)
                        for nm in allNmL:
                            dct[nm] += delta
                        newL = [(ct, (nm, addrMap[nm])) for nm, ct in dct.items()]
                        newTot = sum(dct.values())
                        self.tbl[srcName][tier] = newL
                        self.totTbl[srcName][tier] = newTot

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
#         super(DrawWithReplacementTransferDestinationPolicy, self).__init__(patch,
#                                                                            categoryNameMapper)
        BaseTransferDestinationPolicy.__init__(self, patch, categoryNameMapper)
        self.core = DWRCore(patch)

    def getOrderedCandidateFacList(self, oldFacility, patientAgent, oldTier, newTier, timeNow):
        pairList, tot = self.core.getTierWeightedList(oldFacility.abbrev, newTier)
#             print 'newTier: %s' % CareTier.names[newTier]
        try:
            return [b for a, b in randomOrderByWt(pairList, tot, cull=oldFacility.abbrev)]
        except IndexError, e:
            logger.error('Hit IndexError %s for %s %s -> %s at %s', e, oldFacility.abbrev,
                           oldTier, newTier, timeNow)
            pairList, tot = self.core.getTierWeightedList(oldFacility.abbrev, newTier)
            logger.error('initial tot = %s, pairList = %s', tot, pairList)
            raise

def getPolicyClasses():
    return [DrawWithReplacementTransferDestinationPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(_constants_values,
                                         _constants_schema)
