import numpy as np
import bcolz as bz
import time
import cPickle as pickle
import pandas as pd

if 0:
    from fcntl import lockf, LOCK_UN, LOCK_EX

    LockFile = "date.lock"
    DateFile = "date.info"

    def initLock():
        with open(LockFile, "w") as f:
            pass

    class Lock(object):
        def __init__(self):
            self.f = open(LockFile, "w")
            lockf(self.f.fileno(), LOCK_EX)

        def __enter__(self):
            return self

        def __exit__(self, type, value, traceback):
            lockf(self.f.fileno(), LOCK_UN)
            self.f.close()

            
from lockserver import Lock


            
def writeDateFile(s):
    with Lock('taumod'):
        with open(DateFile, "w") as f:
            f.write(s)

def appendDateFile(s):
    with Lock('taumod'):
        with open(DateFile, "a") as f:
            f.write(s)

def readDateFile():
    with Lock('taumod'):
        with open(DateFile, "r") as f:
            return f.readlines()



# utilities for the pyrhea workers
def getNewTauDict(lastDate):
    while(True):
        lines = readDateFile()
        nextDate = int(lines[0])
        if nextDate <= lastDate:
            time.sleep(1)
            continue

        tauFileName = lines[1].strip()
        with open(tauFileName, "rb") as f:
            tauDict = pickle.load(f)
        return nextDate, tauDict

# it's really dumb having the taufile and the expected prevalences being passed in this way but
# it's easier to get them from a running instance of pyrhea than to generate them in a separate process.
def updateInfo(uid, prevStatsFile, tauFile, expectedFile):
    s = "%s %s %s %s\n"%(uid, prevStatsFile, tauFile, expectedFile)
    appendDateFile(s)


TauFileName = "taudict.pkl"
    
class TauMod(object):
    def __init__(self, workerCount, updatePeriod, dayList):
        self.workerCount = workerCount    # number of processes we're waiting for responses from
        self.updatePeriod = updatePeriod  # how frequently we want data from those processes
        self.dayList = dayList            # list of days in the past that should be used for generating statistics
                                          # ie if timeNow were 50 and dayList were [0,7,14], we'd use stats from
                                          # days [50, 43, 36]

        self.nextDay = 0

        # factors for calculating tau adjustment
        # probably should be passed in.
        self.minRatio = 0.5
        self.maxRatio = 2.0
        self.adjFactor = 0.05


        self.createWorkFile({})


        self.overridePrevalence = {
            ('ALL','HOSP'):0.005,
            ('ALL','ICU'):0.03,
            ('ALL', 'NURSING'): 1.5,
            ('ALL', 'SKILNRS'): 1.5,
        }
        

        
    def getExpected(self, fac, tier):
        if (fac,tier) in self.overridePrevalence:
            return self.overridePrevalence[(fac,tier)]
        if ('ALL', tier) in self.overridePrevalence:
            return self.overridePrevalence[('ALL',tier)]

        return self.expectedPrevalence[(fac, tier)]


    def createWorkFile(self, tauFile):
        self.nextDay += self.updatePeriod
        initLock()

        # create initial dummy taudict
        with open(TauFileName, "w") as f:
            pickle.dump(tauFile, f)

        s = "%s\n%s\n"%(self.nextDay, TauFileName)
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
        for line in lines[2:]:
            uid, prevStatsFileName, tauFileName, expectedFileName = line.strip().split(" ")

            prevStats.append(pd.read_msgpack(prevStatsFileName))
            with open(tauFileName, "rb") as f:
                tauDicts.append(pickle.load(f))

            with open(expectedFileName, "rb") as f:
                expectedDicts.append(pickle.load(f))

        return prevStats, tauDicts, expectedDicts

                
    def isAdjustable(self, fac, tier):
        # this is something that should be passed in but for now:
        if tier == 'LTAC':
            return True
        if tier == 'ICU':
            return True
        return False


    def collatePrev(self, prevStats):
        days = [self.nextDay - d for d in self.dayList]
        
        s = pd.concat(prevStats, ignore_index=True)

        # print overall stats
        tierSums = s[s.day.isin(days)].groupby(['tier']).sum()
        for tier in ['HOSP', 'ICU', 'LTAC', 'NURSING', 'SKILNRS', 'VENT']:
            print "%s Prevalence: %s, (colonized %s, total %s)"%(tier,
                                                                 float(tierSums['COLONIZED'][(tier)]) / tierSums['TOTAL'][(tier)],
                                                                 tierSums['COLONIZED'][(tier)],
                                                                 tierSums['TOTAL'][(tier)])

        


        facSums = s[s.day.isin(days)].groupby(['tier', 'fac']).sum()
        return facSums

    def adjustTaus(self, facSums, tauDicts):
        # assume all tauDicts are identical (maybe we should verify this?)
        tauDict = tauDicts[0]
        
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

            newTau = ((ratio - 1) * self.adjFactor + 1) * tauDict[(fac, tier)]


            print "%s %s: prevalence %s, expected %s, ratio %s, tau %s, newTau %s"%(fac, tier, prevalence, expected, ratio, tauDict[(fac,tier)], newTau)
            tauDict[(fac,tier)] = newTau

        return tauDict
            
            

    def process(self):
        print "waiting!"
        self.waitForWorkers()
        print "processing!"
        prevStats, tauDicts, expectedDicts = self.getWorkerStats()
        self.expectedPrevalence = expectedDicts[0]
        
        facSums = self.collatePrev(prevStats)
        tauDict = self.adjustTaus(facSums, tauDicts)
        self.createWorkFile(tauDict)


def main():
    tm = TauMod(workerCount=5, updatePeriod=50, dayList=[0,5,10,15,20,25,30,35])

    while(1):
        tm.process()
    
        
if __name__ == "__main__":
    main()
