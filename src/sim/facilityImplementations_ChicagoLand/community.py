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

from scipy.stats import expon
from collections import defaultdict
from phacsl.utils.collections.phacollections import enum, SingletonMetaClass
from stats import CachedCDFGenerator
from pathogenbase import PthStatus
from typebase import PatientOverallHealth, CareTier
import facilitybase
import genericCommunity
import pyrheautils

# Setting _constants_schema here does not work because of order of operations
# _constants_schema = 'community_ChicagoLand_constants_schema.yaml'

class CommunityManagerCore(object):
    """This is a place to put infrastructure we must share between community managers"""
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        self.rateScale = 1.0
        self.p = 1.0
        self.Q = _constants['kalmanQ']['value']
        self.H = _constants['kalmanH']['value']

    def kalmanUpdate(self, totPop, meanPop, callerAbbrev):
        """
        Perform a Kalman update of the rate scaling factor.  The nomenclature is
        from Welch & Bishop, "An Introduction to the Kalman Filter" and
        http://scipy-cookbook.readthedocs.io/items/KalmanFiltering.html .

        x is the scale factor for the rate constant (to be estimated),
          expected to be near 1.0
        z is the population delta: z = totPop - meanPop
        H is dz/dx, an input constant
        So our model equation is z = H(x - 1.0)
        p is the estimate error covariance- to be estimated
        A is dxNew/dxOld = 1.0
        Q is the process variance, also an input constant
        R is the measurement variance, which we expect scales as 1/totPop
          since the standard deviation scales as 1/sqrt(N)
        W and V are 1.0
        """

        x = self.rateScale
        P = self.p
        z = totPop - meanPop
        Q = self.Q
        A = 1.0
        H = self.H
        R = 1.0/totPop
        xHatMinus = x
        PMinus = A*P*A + Q
        K = (PMinus * H) / (H * PMinus * H + R)
        xHat = xHatMinus + (K * (z - H*(xHatMinus - 1.0)))
        P = (1.0 - K * H) * PMinus
#         print ('Kalman update for %s: x= %s P=%s z=%s R=%s -> x=%s P=%s' %
#                (callerAbbrev, x, self.p, z, R, xHat, P))
        self.rateScale = xHat
        self.p = P


class CommunityManager(genericCommunity.CommunityManager):
    def __init__(self, name, patch, facility):
        super(CommunityManager, self).__init__(name, patch, facility)
        self.core = CommunityManagerCore()
        self.meanPop = self.fac.meanPop
        self.lastKalmanUpdateTime = 0

    def getProbThaw(self, patientCategory, dT):
        return self.core.rateScale * self.fac.cachedCDFs[patientCategory].intervalProb(0, dT)

    def perTickActions(self, timeNow):
        super(CommunityManager, self).perTickActions(timeNow)
        if timeNow != self.lastKalmanUpdateTime:
            totPop = sum([sum([len(frz.frozenAgentList) for frz in ward.freezers.values()])
                          for ward in self.fac.getWards()])
            self.core.kalmanUpdate(totPop, self.meanPop, self.fac.name)
            self.lastKalmanUpdateTime = timeNow


class CommunityWard(genericCommunity.CommunityWard):
    def classify(self, agent, timeNow):
        """Return the PatientCategory appropriate for this agent"""
        health = PatientOverallHealth.names[agent.getStatus().overall]
        if agent.getStatus().pthStatus == PthStatus.COLONIZED:
            return '%s_colonized' % health
        elif agent.getStatus().pthStatus == PthStatus.CLEAR:
            return '%s_base' % health
        elif agent.getStatus().pthStatus == PthStatus.UNDETCOLONIZED:
            return '%s_undetcolonized' % health
        else:
            raise genericCommunity.FreezerError('%s has unexpected PthStatus %s' %
                                                (agent.name,
                                                 PthStatus.names[agent._status.pthStatus]))


class CommunityCore(object):
    """This is a place to put infrastructure we must share between communities"""
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        self.tbl = None

    def _buildWeightedLists(self, ctgNameMapper):
        self.tbl = {}
        for transferFilePath in _constants['srcToCategoryMapFilePaths']:
            logger.info('Importing the weight data file %s', transferFilePath)
            rawTbl = pyrheautils.importConstants(transferFilePath,
                                                 _constants['srcToCategoryMapFileSchema'])
            # normalize on the fly
            for srcName, rec in rawTbl.items():
                tot = sum(rec.values())
                if srcName not in self.tbl:
                    self.tbl[srcName] = {}
                for ctg, ct in rec.items():
                    self.tbl[srcName][ctgNameMapper(ctg)] = ct/tot
            logger.info('Import complete.')

        # We need to add an entry for 'no source' representing the net weights
        # marginalized across sources
        dct = defaultdict(lambda: 0.0)
        for rec in self.tbl.values():
            for ctg, ct in rec.items():
                dct[ctg] += ct
        tot = sum(dct.values())
        self.tbl[None] = {ctg: ct / tot for ctg, ct in dct.items()}
        logger.info('Post-processing of weight data is complete')

    def getCategoryWeights(self, pFAbbrev, ctgNameMapper):
        """
        Return net fraction of patients transferring to HOSPITAL, NURSINGHOME/SNF,
        VSNF, and LTAC in that order as a tuple.  The sum of the return values is
        expected to be 1.0.
        """
        if self.tbl is None:
            self._buildWeightedLists(ctgNameMapper)
        rec = self.tbl[pFAbbrev]
        rslt = [rec[key] if key in rec else 0.0
                for key in ['HOSPITAL', 'NURSINGHOME', 'VSNF', 'LTAC']]
        return tuple(rslt)

class Community(genericCommunity.Community):
    def __init__(self, descr, patch, policyClasses=None, categoryNameMapper=None):
        self.meanPop = descr['meanPop']['value']
        super(Community, self).__init__(descr, patch, policyClasses=policyClasses,
                                        categoryNameMapper=categoryNameMapper,
                                        managerClass=CommunityManager,
                                        wardClass=CommunityWard)
        self.core = CommunityCore()

    def calcTierRateConstants(self, pFAbbrev):
        """
        These constants represent the fraction of the newly-sick that die, that
        need rehab as a (v)SNF, and that need care at the HOSP, ICU, LTAC, SKILNRS,
        and VENT tiers respectively.  The return values should sum to 1.0.
        """
        (deathRate, needsRehabRate, needsHospRate, needsICURate, needsLTACRate,
         needsSkilNrsRate, needsVentRate) = (super(Community, self)
                                             .calcTierRateConstants(pFAbbrev))

        # un-weight to adjust for deathRate
        liveRate = 1.0 - deathRate
        assert liveRate > 0.0, 'Surely not *everyone* dies!'
        needsRehabRate /= liveRate
        needsHospRate /= liveRate
        needsICURate /= liveRate
        needsLTACRate /= liveRate
        needsSkilNrsRate /= liveRate
        needsVentRate /= liveRate

        # Adjust for this particular patient's history
        wtHOSP, wtSNF, wtVSNF, wtLTACH = self.core.getCategoryWeights(pFAbbrev,
                                                                      self.categoryNameMapper)
        needsVSNFTotRate = needsRehabRate + needsSkilNrsRate + needsVentRate
        if needsVSNFTotRate > 0.0:
            fracVSNFBedsRehab = needsRehabRate / needsVSNFTotRate
            fracVSNFBedsVent = needsVentRate / needsVSNFTotRate
        else:
            fracVSNFBedsRehab = 0.0
            fracVSNFBedsVent = 0.0
        fracVSNFBedsSkil = 1.0 - (fracVSNFBedsRehab + fracVSNFBedsVent)
        fracHospBedsICU = needsICURate/(needsHospRate + needsICURate)
        fracHospBedsHosp = 1.0 - fracHospBedsICU

        newSkilNrsRate = fracVSNFBedsSkil * wtVSNF
        newVentRate = fracVSNFBedsVent * wtVSNF
        newRehabRate = wtSNF + fracVSNFBedsRehab*wtVSNF
        newLTACRate = wtLTACH
        newHospRate = fracHospBedsHosp * wtHOSP
        newICURate = fracHospBedsICU * wtHOSP

        # re-weight to adjust for deathRate
        needsRehabRate = newRehabRate * liveRate
        needsHospRate = newHospRate * liveRate
        needsICURate = newICURate * liveRate
        needsLTACRate = newLTACRate * liveRate
        needsSkilNrsRate = newSkilNrsRate * liveRate
        needsVentRate = newVentRate * liveRate

        return (deathRate, needsRehabRate, needsHospRate, needsICURate, needsLTACRate,
                needsSkilNrsRate, needsVentRate)


def generateFull(facilityDescr, patch, policyClasses=None, categoryNameMapper=None):
    return genericCommunity.generateFull(facilityDescr, patch, policyClasses=policyClasses,
                                         categoryNameMapper=categoryNameMapper,
                                         communityClass=Community)


genericCommunity.importCommunity(__name__)
