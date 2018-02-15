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
from phacsl.utils.collections.phacollections import SingletonMetaClass, enum
from policybase import TransferDestinationPolicy as BaseTransferDestinationPolicy
from facilitybase import CareTier, tierToQueueMap
import facilitybase
import transferbydrawwithreplacement

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

    def _buildWeightedLists(self):
        nmDict = CareTier.names
        tierFacSets = {tier: set(self.getTierAddrMap(tier).keys()) for tier in nmDict.keys()}
        pairsSeen = set()
        self.tbl = {}
        self.totTbl = {}
        for transferMatrixFilePath in _constants['transferFilePaths']:
            LOGGER.info('Importing the weight data file %s', transferMatrixFilePath)
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
            LOGGER.info('Import complete.')
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
        self.fallbackPolicy = (transferbydrawwithreplacement.
                               DrawWithReplacementTransferDestinationPolicy(patch,
                                                                            categoryNameMapper))

    def classifyTrans(self, newTier, timeTupleL):
        """
        Given a list returned by buildTimeTupleList, what type of transfer is this?
        Returns a TransType enum, the relevant starting facility, and the days since
        discharge from that facility
        """
        hospCategoryL = ['HOSPITAL', 'LTAC', 'LTACH']
        curDays, curFac, curCategory, curTier = timeTupleL[0]
        if curCategory in hospCategoryL:
            return TransType.OTHER, curFac, 0
        else:
            dayTot = 0
            firstHOffset = None
            firstHDays = None
            firstNHOffset = None
            firstNHDays = None
            firstComOffset = None
            firstComDays = None
            for offset in xrange(len(timeTupleL)):
                days, fac, cat, tier = timeTupleL[offset]  # @UnusedVariable
                if cat in hospCategoryL:
                    firstHOffset = offset
                    firstHDays = dayTot
                    break
                elif cat == 'COMMUNITY' and firstComOffset is None:
                    firstComOffset = offset
                    firstComDays = dayTot
                    firstNHOffset = None  # Only count NH's which come before home
                    firstNHDays = None
                elif cat != 'COMMUNITY' and firstNHOffset is None:
                    firstNHOffset = offset
                dayTot += days
                if dayTot > 365:
                    return TransType.OTHER, curFac, 0
            else:
                # This patient has never been to a hospital
                return TransType.OTHER, curFac, 0
            if firstComOffset is None:
                # Patient hasn't been home since last hosp visit, so this cannot be
                # an indirect transfer
                return TransType.OTHER, curFac, 0
            elif firstNHOffset is None:
                assert firstHOffset is not None, 'logic error'
                return TransType.HOSP_INDIRECT, timeTupleL[firstHOffset][1], firstHDays
            else:
                # This is believed to count as NH-Readmit from the last NH before going home.
                return TransType.NH_READMIT, timeTupleL[firstNHOffset][1], firstNHDays

    def getOrderedCandidateFacList(self, thisFacility, patientAgent,  # @UnusedVariable
                                   oldTier, newTier, timeNow):  # @UnusedVariable
        timeTupleL = facilitybase.buildTimeTupleList(patientAgent, timeNow)
        transType, rootFac, daysSince = self.classifyTrans(newTier, timeTupleL)
        if transType == TransType.OTHER:
            return self.fallbackPolicy.getOrderedCandidateFacList(thisFacility, patientAgent,
                                                                  oldTier, newTier, timeNow)
        else:
            if transType in (TransType.HOSP_INDIRECT, TransType.NH_READMIT):
                pairList, tot = self.core.getWeightedList(rootFac, newTier)
            else:
                raise RuntimeError('Unknown transfer type %s' % transType)
#             print '%s %s -> %s: tot=%d' % (thisFacility.name, CareTier.names[oldTier],
#                                            CareTier.names[newTier], tot)
#             print 'pairList: %s' % str([(wt, tpl[0]) for wt, tpl in pairList])
#             print 'newTier: %s' % CareTier.names[newTier]

            try:
                return [addr for destNm, addr                                   # @UnusedVariable
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
