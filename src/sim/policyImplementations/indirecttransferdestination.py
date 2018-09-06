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
from collections import defaultdict

import pyrheautils
from phacsl.utils.collections.phacollections import SingletonMetaClass, enum
from policybase import TransferDestinationPolicy as BaseTransferDestinationPolicy
from facilitybase import CareTier, tierToQueueMap
import facilitybase
from pyrheabase import TierUpdateModKey
import transferbydrawwithreplacement
import categorydrawwithreplacement
from genericCommunity import BYPASS_KEY

_validator = None
_constants_values = '$(MODELDIR)/constants/indirecttransferdestination_constants.yaml'
_constants_schema = 'indirecttransferdestination_constants_schema.yaml'
_constants = None

LOGGER = logging.getLogger(__name__)

TransType = enum('HOSP_INDIRECT', 'NH_READMIT', 'OTHER')

class IndirectTransferCore(object):
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
        LOGGER.info('Building facility-to-address map')
        for tier in nmDict.keys():
            self.tierAddrMap[tier] = {}
            serviceTupleList = self.patch.serviceLookup(tierToQueueMap[tier].__name__)
            for info, addr in serviceTupleList:
                innerInfo, facAbbrev, facCoords = info  # @UnusedVariable
                self.tierAddrMap[tier][facAbbrev] = addr
        LOGGER.info('map complete')

    def getTierAddrMap(self, tier):
        if self.tierAddrMap is None:
            self._buildTierAddrMap()
        return self.tierAddrMap[tier]

    def _buildWtLists(self, tmpTbl):
        for srcName in tmpTbl:
            for tier in tmpTbl[srcName]:
                wtL = []
                for destName, ct in tmpTbl[srcName][tier].items():
                    wtL.append((ct, (destName, self.getTierAddrMap(tier)[destName])))
                self.tbl[srcName][tier] = wtL

    def _addRecToTbl(self, srcName, rec, tmpTbl, tierFacSets):
        if srcName not in self.tbl:
            self.tbl[srcName] = {}
            tmpTbl[srcName] = {}
            self.totTbl[srcName] = {}
        for tier in CareTier.names.keys():
            if tier not in self.tbl[srcName]:
                self.tbl[srcName][tier] = {}
                tmpTbl[srcName][tier] = defaultdict(lambda: 0)
                self.totTbl[srcName][tier] = 0.0
            wtSum = self.totTbl[srcName][tier]
            for destName, ct in rec.items():
                if destName in tierFacSets[tier]:
                    tmpTbl[srcName][tier][destName] += ct
                    wtSum += ct
            self.totTbl[srcName][tier] = wtSum

    def _buildWeightedLists(self):
        tierFacSets = {tier: set(self.getTierAddrMap(tier).keys())
                       for tier in CareTier.names.keys()}
        self.tbl = {}
        tmpTbl = {}
        self.totTbl = {}
        for transferFilePath in _constants['transferFilePaths']:
            LOGGER.info('Importing the weight data file %s', transferFilePath)
            rawTbl = pyrheautils.importConstants(transferFilePath,
                                                 _constants['transferFileSchema'])
            for srcName, rec in rawTbl.items():
                self._addRecToTbl(srcName, rec, tmpTbl, tierFacSets)
            LOGGER.info('Import complete.')
        for transferFilePath in _constants['bypassTransferFilePaths']:
            LOGGER.info('Importing the weight data file %s', transferFilePath)
            rawTbl = pyrheautils.importConstants(transferFilePath,
                                                 _constants['transferFileSchema'])
            assert rawTbl.keys() == ['COMMUNITY'], ("Bypass transfer file {0}".format(transferFilePath)
                                                    + " has an unexpected format")
            for rec in rawTbl.values():
                self._addRecToTbl(BYPASS_KEY, rec, tmpTbl, tierFacSets)
            LOGGER.info('Import complete.')
            
        self._buildWtLists(tmpTbl)
        for srcName, subTbl in self.tbl.items():
            for wtL in subTbl.values():
                wtL.sort(reverse=True)
        LOGGER.info('Post-processing of weight data files is complete.')

    def getWeightedList(self, srcName, tier):
        """
        Returns a tuple containing:
        - a list of the form [(wt, (abbrev, addr)), ...] to be shuffled
        - a value wtSum equal to the sum of the weights in the list
        Because the srcName is specific to the Hospital Indirect or NH Readmit case, we
        can store both types of entry in the same table.
        """
        if self.tbl is None:
            self._buildWeightedLists()
        return (self.tbl[srcName][tier][:], self.totTbl[srcName][tier])


class IndirectTransferDestinationPolicy(BaseTransferDestinationPolicy):
    def __init__(self, patch, categoryNameMapper):
        super(IndirectTransferDestinationPolicy, self).__init__(patch, categoryNameMapper)
        self.core = IndirectTransferCore(patch)
        self.cache = {}
        self.fallbackPolicy = (categorydrawwithreplacement.
                               CategoryDrawWithReplacementTransferDestinationPolicy(patch,
                                                                            categoryNameMapper))

    def getOrderedCandidateFacList(self, thisFacility, patientAgent,  # @UnusedVariable
                                   oldTier, newTier, modifierDct, timeNow):  # @UnusedVariable
        assert TierUpdateModKey.FLOW_KEY in modifierDct, 'FLOW_KEY has not been set'
        flowKey = modifierDct[TierUpdateModKey.FLOW_KEY]
        if flowKey is None:
            return self.fallbackPolicy.getOrderedCandidateFacList(thisFacility, patientAgent,
                                                                  oldTier, newTier,
                                                                  modifierDct, timeNow)
        else:
            try:
                pairList, tot = self.core.getWeightedList(flowKey, newTier)
                return [addr for destNm, addr                               # @UnusedVariable
                        in transferbydrawwithreplacement.randomOrderByWt(pairList, tot,
                                                                         cull=thisFacility.abbrev)]
            except IndexError, e:
                LOGGER.error('Hit IndexError %s for %s %s -> %s at %s', e, thisFacility.abbrev,
                             oldTier, newTier, timeNow)
                raise


def getPolicyClasses():
    return [IndirectTransferDestinationPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(_constants_values,
                                         _constants_schema)
