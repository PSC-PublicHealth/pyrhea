from pathogenbase import PthStatus
from facilitybase import CareTier

import numpy as np
import bcolz as bz
import time
import os


class Monitor(object):
    def __init__(self, patch, totalRunDays, stopFn = None, nextStop = None):
        self.patch = patch
        self.nextStop = nextStop
        self.stopFn = stopFn
        dtype = [('fac', 'S20'), ('tier', 'S10'), ('ward', 'uint32'), ('day', 'uint32')]
        for i in xrange(len(PthStatus.names)):
            dtype.append((PthStatus.names[i], 'uint32'))
        dtype.append(('TOTAL', 'uint32'))
        ra = np.recarray((0,), dtype = dtype)
        self.pthData = bz.ctable(ra)
        self.pthDataDF = None    # pthData as a pandas dataframe
        self.uniqueID = str(os.getpid()) + "_" + str(patch.patchId) + "_" + str(time.time())

    def XXXgetPthData(self, timeNow):
        """
        traverses each ward of each facility and updates pthData
        """

        for fac in self.patch.allFacilities:
            for ward in fac.getWards():
                # for now let's drop any "HOME" wards
                if ward.tier == CareTier.HOME:
                    continue
                pPC = ward.iA.getPatientPthCounts(timeNow)
                for pthStatusEnum, count in pPC.items():
                    pthStatus = PthStatus.names[pthStatusEnum]
                    self.pthData[fac.abbrev][pthStatus][timeNow] += count
                    self.pthData[fac.abbrev]['total'][timeNow] += count

    def collectPthData(self, timeNow):
        """
        traverses each ward of each facility and updates pthData
        """

        for fac in self.patch.allFacilities:
            for ward in fac.getWards():
                # for now let's drop any "HOME" wards
                if ward.tier == CareTier.HOME:
                    continue
                pPC = ward.iA.getPatientPthCounts(timeNow)
                counts = []
                total = 0
                for i in xrange(len(PthStatus.names)):
                    counts.append(pPC[i])
                    total += pPC[i]
                counts.append(total)
                row = [fac.abbrev, CareTier.names[ward.tier], ward.wardNum, timeNow] + counts
                self.pthData.append(row)

        self.pthDataDF = None

    def getPthData(self):
        if self.pthDataDF is None:
            self.pthDataDF = self.pthData.todataframe()
        return self.pthDataDF


    def setStopTimeFn(self, when, fn):
        """
        set a time on the simulation when we should stop and wait for an update
        """
        self.nextStop = when
        self.stopFn = fn

    def daily(self, timeNow):
        self.collectPthData(timeNow)
        print "next stop: %s, timeNow: %s"%(self.nextStop, timeNow)
        if self.nextStop is not None:
            if timeNow >= self.nextStop:
                self.nextStop = self.stopFn(timeNow, self)

    def createDailyCallbackFn(self):
        def fn(loop, timeNow):
            self.daily(timeNow)
        return fn

    def writeData(self, fileName):
        pthData = self.getPthData()
        fileName = "%s_%s.bcz"%(fileName, self.uniqueID)
        pthData.to_msgpack(fileName, compress="zlib")
