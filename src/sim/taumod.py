#! /usr/bin/env python

###################################################################################
# Copyright   2018, Pittsburgh Supercomputing Center (PSC).  All Rights Reserved. #
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

from __future__ import print_function
import sys
import time
from math import sqrt
from collections import defaultdict
import six.moves.cPickle as pickle
import logging
import optparse
import pandas as pd
import yaml

from lockserver import Lock, SharedLock
from phacsl.utils.collections.phacollections import SingletonMetaClass


CONFIG_FNAME = 'taumod_config.yaml'

LOGGER = logging.getLogger(__name__)

# Golden ratio, for use during iterative optimization
GOLDEN_PHI = (1 + sqrt(5.0))/2.0
GOLDEN_RESPHI = 2.0 - GOLDEN_PHI

class Config(object):
    """Seems like a long way to go to avoid globals, but I'll play along"""
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        with open(CONFIG_FNAME, 'rU') as f:
            self.config = yaml.load(f)

    def __getitem__(self, key):
        return self.config[key]


def writeDateFile(s):
    with Lock('taumod'):
        with open(Config()['DateFile'], "w") as f:
            f.write(s)


def appendDateFile(s):
    with Lock('taumod'):
        with open(Config()['DateFile'], "a") as f:
            f.write(s)


def readDateFile():
    # give time for taumod to have started up and dropped the date file
    for i in xrange(60):
        try:
            with SharedLock('taumod'):
                with open(Config()['DateFile'], "r") as f:
                    return f.readlines()
        except:
            pass
        time.sleep(5)
        LOGGER.info("retry read datefile")
    raise RuntimeError("Can't read date file")

# utilities for the pyrhea workers
def getNewTauDict(lastDate):
    LOGGER.info('waiting')
    while True:
        lines = readDateFile()
        nextDate = int(lines[0])
        if nextDate <= lastDate:
            time.sleep(1)
            continue

        LOGGER.info('done waiting')
        tauFileName = lines[1].strip()
        with SharedLock(tauFileName):
            with open(tauFileName, "rb") as f:
                tauDict = pickle.load(f)
        return nextDate, tauDict

# it's really dumb having the taufile and the expected prevalences being passed in this way but
# it's easier to get them from a running instance of pyrhea than to generate them in a separate process.
def updateInfo(uid, prevStatsFile, tauFile, expectedFile):
    s = "%s %s %s %s\n"%(uid, prevStatsFile, tauFile, expectedFile)
    appendDateFile(s)


class TauMod(object):
    def __init__(self, workerCount=None):
        """
        workerCount is the expected number of workers, defaulting to the config file value.
        All other adjustable parameters are drawn from the config file
        """

        # number of processes we're waiting for responses from
        if workerCount is None:
            self.workerCount = Config()['WorkerCount']
        else:
            self.workerCount = workerCount

        # how frequently we want data from those processes
        self.updatePeriod=Config()['UpdatePeriod']

        # list of days in the past that should be used for generating statistics
        # ie if timeNow were 50 and dayList were [0,7,14], we'd use stats from
        # days [50, 43, 36]
        self.dayList = Config()['DayList']

        self.nextDay = 0
        self.endDay = Config()['EndDay']

        # factors for calculating tau adjustment
        self.minRatio = Config()['MinRatio']
        self.maxRatio = Config()['MaxRatio']
        self.adjFactor = Config()['AdjFactor']

        self.createWorkFile({})

        self.overridePrevalence = {(ent['fac'], ent['tier']): ent['value']
                                   for ent in Config()['OverridePrevalence']}
        self.isAdjustableDict = {(ent['fac'], ent['tier']): ent['value']
                                 for ent in Config()['IsAdjustable']}
        self.pastData = None
        self.current_best_tau = None
        self.nUpdates = 0

    def getExpected(self, fac, tier):
        if (fac,tier) in self.overridePrevalence:
            return self.overridePrevalence[(fac,tier)]
        if ('ALL', tier) in self.overridePrevalence:
            return self.overridePrevalence[('ALL',tier)]

        return self.expectedPrevalence[(fac, tier)]


    def createWorkFile(self, tauFile):
        self.nextDay += self.updatePeriod

        # create initial dummy taudict
        tFN = Config()['TauFileName']
        with Lock(tFN):
            with open(tFN, "w") as f:
                pickle.dump(tauFile, f)

        s = "%s\n%s\n"%(self.nextDay, tFN)
        writeDateFile(s)

    def waitForWorkers(self):
        while True:
            lines = readDateFile()
            if len(lines) - 2 == self.workerCount:
                return
            time.sleep(1)

    def getWorkerStats(self):
        lines = readDateFile()
        prevStats = []
        tauDicts = []
        expectedDicts = []
        days = [self.nextDay - d for d in self.dayList]
        for line in lines[2:]:
            uid, prevStatsFileName, tauFileName, expectedFileName = line.strip().split(" ") # @UnusedVariable

            pS = pd.read_msgpack(prevStatsFileName)
            prevStats.append(pS[pS.day.isin(days)])
            with SharedLock(tauFileName):
                with open(tauFileName, "rb") as f:
                    tauDicts.append(pickle.load(f))

            with SharedLock(expectedFileName):
                with open(expectedFileName, "rb") as f:
                    expectedDicts.append(pickle.load(f))

        return prevStats, tauDicts, expectedDicts

    def isAdjustable(self, fac, tier):
        if (fac, tier) in self.isAdjustableDict:
            return self.isAdjustableDict[(fac, tier)]
        elif ('ALL', tier) in self.isAdjustableDict:
            return self.isAdjustableDict[('ALL', tier)]
        else:
            return False

    def collatePrev(self, prevStats):

        s = pd.concat(prevStats, ignore_index=True)

        # print overall stats
        tierGroups = s.groupby(['tier'])
        tierSums = tierGroups.sum()
        for tier in tierGroups.groups.keys():
            colRatio = float(tierSums['COLONIZED'][(tier)]) / tierSums['TOTAL'][(tier)]
            print("%s Prevalence: %s, (colonized %s, total %s)") % (tier,
                                                                    colRatio,
                                                                    tierSums['COLONIZED'][(tier)],
                                                                    tierSums['TOTAL'][(tier)])

        return s

    def adjustTaus(self, facDF, tauDicts):
        # assume all tauDicts are identical (maybe we should verify this?)
        tauDict = tauDicts[0]
        facGrps = facDF.groupby(['tier', 'fac'])
        facSums = facGrps.sum()

        facDF['prev_sample'] = facDF['COLONIZED'].astype(float)/facDF['TOTAL'].astype(float)
        facSampGrps = facDF.groupby(['tier', 'fac'])
        facSampMax = facSampGrps.max()['prev_sample']
        facSampMin = facSampGrps.min()['prev_sample']
        facSampMean = facSampGrps.mean()['prev_sample']
        facSampMedian = facSampGrps.median()['prev_sample']
        facSampStdv = facSampGrps.std()['prev_sample']
        facSampQ1 = facSampGrps.quantile(q=0.25)['prev_sample']
        facSampQ3 = facSampGrps.quantile(q=0.75)['prev_sample']

        #algorithm = 'gr_independent'
        algorithm = 'running_average'

        if algorithm == 'gr_independent':
            if self.current_best_tau is None:
                self.current_best_tau = defaultdict(float)

            for key, colonized in facSums['COLONIZED'].to_dict().items():
                tier, fac = key
                if not self.isAdjustable(fac, tier):
                    continue

                total = facSums['TOTAL'][key]
                if total:
                    prevalence = float(colonized) / total
                else:
                    prevalence = 0.0

                expected = self.getExpected(fac, tier)

                if prevalence == 0.0:
                    ratio = self.maxRatio
                else:
                    ratio = expected / prevalence
                    if ratio > self.maxRatio:
                        ratio = self.maxRatio
                    elif ratio < self.minRatio:
                        ratio = self.minRatio

                    ctr = 0.5 * (ratio + 1.0)
                    ratio = ctr + GOLDEN_RESPHI * (ratio - ctr)

                #newTau = ((ratio - 1) * self.adjFactor + 1) * tauDict[(fac, tier)]
                newTau = min(ratio*tauDict[(fac, tier)], 0.9999)
                self.current_best_tau[(fac, tier)] = (self.current_best_tau[(fac, tier)] + newTau) / 2

                print (("%s %s: total: %s, prevalence %s, expected %s," 
                       " prev_mean %s, prev_stdv %s prev_min %s, prev_q1 %s, prev_q2 %s, prev_q3 %s, prev_max %s,"
                       " ratio %s, tau %s, newTau %s, current_best_tau %s")
                        % (fac, tier, total, prevalence, expected,
                           facSampMean[key], facSampStdv[key], facSampMin[key],
                           facSampQ1[key], facSampMedian[key], facSampQ3[key],
                           facSampMax[key], ratio, tauDict[(fac,tier)], newTau,
                           self.current_best_tau[(fac, tier)]))

                tauDict[(fac, tier)] = newTau

        elif algorithm == 'running_average':
            if self.pastData is None:
                self.pastData = defaultdict(float)

            if self.current_best_tau is None:
                self.current_best_tau = defaultdict(float)

            icuTauMultiplier = 1.10

            for key, colonized in facSums['COLONIZED'].to_dict().items():
                tier, fac = key
                if not self.isAdjustable(fac, tier):
                    continue
                total = facSums['TOTAL'][key]

                if total:
                    prevalence = float(colonized) / total
                else:
                    prevalence = 0.0

                expected = self.getExpected(fac, tier)

                if self.nUpdates == 0:
                    prevalence = float(colonized)/total
                elif tier == 'ICU':
                    prevalence = (1.0/3.0) * (self.pastData[(fac, tier)]
                                              + min(1.0,
                                                    icuTauMultiplier*self.pastData[(fac,'HOSP')])
                                              + prevalence)
                else:
                    prevalence = 0.5 * (prevalence + self.pastData[(fac, tier)])
                self.pastData[(fac, tier)] = prevalence

                if prevalence == 0.0:
                    ratio = self.maxRatio
                else:
                    ratio = expected / prevalence
                    if ratio > self.maxRatio:
                        ratio = self.maxRatio
                    elif ratio < self.minRatio:
                        ratio = self.minRatio

                ctr = 0.5 * (ratio + 1.0)
                ratio = ctr + GOLDEN_RESPHI * (ratio - ctr)

                newTau = min(ratio*tauDict[(fac, tier)], 0.9999)
                self.current_best_tau[(fac, tier)] = (self.current_best_tau[(fac, tier)] + newTau) / 2

                print (("%s %s: total: %s, prevalence %s, expected %s," 
                       " prev_mean %s, prev_stdv %s prev_min %s, prev_q1 %s, prev_q2 %s, prev_q3 %s, prev_max %s,"
                       " ratio %s, tau %s, newTau %s, current_best_tau %s")
                        % (fac, tier, total, prevalence, expected,
                           facSampMean[key], facSampStdv[key], facSampMin[key],
                           facSampQ1[key], facSampMedian[key], facSampQ3[key],
                           facSampMax[key], ratio, tauDict[(fac,tier)], newTau,
                           self.current_best_tau[(fac, tier)]))

                tauDict[(fac, tier)] = newTau

        else:
            raise RuntimeError('unknown tau update algorithm %s' % algorithm)

        self.nUpdates += 1
        return tauDict


    def process(self):
        LOGGER.info("waiting!")
        self.waitForWorkers()
        LOGGER.info("processing!")
        prevStats, tauDicts, expectedDicts = self.getWorkerStats()
        self.expectedPrevalence = expectedDicts[0]

        facDF = self.collatePrev(prevStats)
        tauDict = self.adjustTaus(facDF, tauDicts)
        self.createWorkFile(tauDict)


def main():
    parser = optparse.OptionParser(usage="""
    %prog [--debug] [-N nWorkers]
    """)
    parser.add_option("-d", "--debug", action="store_true",
                      help="debugging output")
    parser.add_option("-N", "--workers", action="store", type="int",
                      help="expected number of workers (optional)")
    opts, args = parser.parse_args()
    if args:
        parser.error('Invalid arguments %s' % args)
    if opts.debug:
        logLvl = 'debug'
    else:
        logLvl = 'info'
    logging.basicConfig(level=getattr(logging, logLvl.upper(), None))

    if opts.workers:
        nWorkers = opts.workers
    else:
        nWorkers = None
    tm = TauMod(workerCount=nWorkers)

    while tm.nextDay <= tm.endDay:
        tm.process()
        


if __name__ == "__main__":
    main()
