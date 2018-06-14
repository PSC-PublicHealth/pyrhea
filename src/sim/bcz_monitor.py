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

from pathogenbase import PthStatus
from facilitybase import CareTier

import numpy as np
import bcolz as bz
import time
import os
import logging

LOGGER = logging.getLogger(__name__)


# this is a dictionary linking anything 
TrackableValues = {}
PreregisteredKeys = set()

def addTrackable(key, fn, dtypes):
    global TrackableValues
    if key in TrackableValues:
        pass
#        logger.warning("TrackableValue keyword collision: %s"%key)
    TrackableValues[key] = (fn, dtypes)

def preregisterKey(key):
    global PreregisteredKeys
    PreregisteredKeys.add(key)

class Monitor(object):
    def __init__(self, filename, patch, stopFn=None, nextStop=None):
        self.patch = patch
        self.nextStop = nextStop
        self.stopFn = stopFn
        self.filename = filename
        self.dtype = [('fac', 'S20'), ('tier', 'S10'), ('ward', 'uint32'), ('day', 'uint32')]
        for i in xrange(len(PthStatus.names)):
            self.dtype.append((PthStatus.names[i], 'uint32'))
        self.dtype.append(('TOTAL', 'uint32'))
        self.pthDataDF = None    # pthData as a pandas dataframe
        self.uniqueID = str(os.getpid()) + "_" + str(patch.patchId) + "_" + str(time.time())
        self.registeredFns = []
        self.registeredKeys = set()

    def registerFields(self, dtypes, fn):
        self.registeredFns.append(fn)
        self.dtype.extend(dtypes)

    def registerTrackable(self, key):
        global TrackableValues
        
        # check if this key is already registered
        if key in self.registeredKeys:
            return
        
        fn, dtypes = TrackableValues[key]
        self.registerFields(dtypes, fn)
        self.registeredKeys.add(key)

    def addPreregistered(self):
        global PreregisteredKeys
        for key in PreregisteredKeys:
            self.registerTrackable(key)

    def solidifyFields(self):
        ra = np.recarray((0,), dtype = self.dtype)
        self.pthData = bz.ctable(ra, rootdir=self.filename, mode="w")
    
    def collectData(self, timeNow):
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
                for fn in self.registeredFns:
                    row.extend(fn(ward, timeNow))
                self.pthData.append(row)

        self.pthDataDF = None
        resetMiscCounters(self.patch, timeNow)

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
        self.collectData(timeNow)
        LOGGER.debug("daily: next stop: %s, timeNow: %s", self.nextStop, timeNow)
        if self.nextStop is not None:
            if timeNow >= self.nextStop:
                self.nextStop = self.stopFn(timeNow, self)

    def createDailyCallbackFn(self):
        def fn(loop, timeNow):
            self.daily(timeNow)
        return fn

    def writeData(self):
        self.pthData.flush()
        
    def XXXwriteData(self, fileName):
        pthData = self.getPthData()
        fileName = "%s_%s.mpz"%(fileName, self.uniqueID)
        pthData.to_msgpack(fileName, compress="zlib")

        
class oldMonitor(object):
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
        LOGGER.debug("daily: next stop: %s, timeNow: %s", self.nextStop, timeNow)
        if self.nextStop is not None:
            if timeNow >= self.nextStop:
                self.nextStop = self.stopFn(timeNow, self)

    def createDailyCallbackFn(self):
        def fn(loop, timeNow):
            self.daily(timeNow)
        return fn

    def writeData(self, fileName):
        pthData = self.getPthData()
        fileName = "%s_%s.mpz"%(fileName, self.uniqueID)
        pthData.to_msgpack(fileName, compress="zlib")



_resetCountNeeded = 1
_resetRequestCount = 0
_resetRequestTime = None

def resetMiscCounters(patch, timeNow):
    """
    stupid hack to reset the miscCounters only after both per day callbacks have been run that need to look at them
    """
    global _resetCountNeeded, _resetRequestCount, _resetRequestTime
    if _resetRequestTime != timeNow:
        _resetRequestCount = 0
    _resetRequestCount += 1
    if _resetRequestCount != _resetCountNeeded:
        return

    # looks like we should reset everything
    for fac in patch.allFacilities:
        for ward in fac.getWards():
            ward.miscCounters = defaultdict(lambda: 0)
    
