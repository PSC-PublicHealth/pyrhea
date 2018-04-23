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

'''
Created on Apr 20, 2018

@author: xuth, welling
'''

import logging
import cPickle as pickle
import taumod
from typebase import CareTier

LOGGER = logging.getLogger(__name__)

def getTauDict(patch):
    ret = {}
    for fac in patch.allFacilities:
        for ward in fac.getWards():
            # for now let's drop any "HOME" wards
            if ward.tier == CareTier.HOME:
                continue
            assert hasattr(ward.iA, 'tau'), 'Ward %s has no tau' % ward._name
            ret[(fac.abbrev, CareTier.names[ward.tier])] = ward.iA.tau
    return ret


def overrideTaus(patch, tauDict):
    for fac in patch.allFacilities:
        flushedFac = False
        for ward in fac.getWards():
            assert hasattr(ward.iA, 'tau'), 'Ward %s has no tau' % ward._name
            key = (fac.abbrev, CareTier.names[ward.tier])
            if key in tauDict:
                if not flushedFac:
                    fac.flushCaches()
                    flushedFac = True
                ward.iA.flushCaches()
                ward.iA.tau = tauDict[key]


class TauAdjuster(object):
    '''
    Handle the interface to the process which is updating run parameters
    '''
    def __init__(self, monitor):
        '''
        monitor is a Monitor object, which collects data and periodically pauses for updates.
        '''
        self.patch = monitor.patch
        self.monitor = monitor
        self.expectedPrevalence = self.getColonizedTargets()
        self.tauHistory = {}

        self.nextDate, tauDict = taumod.getNewTauDict(-1)
        overrideTaus(self.patch, tauDict)

    def getColonizedTargets(self):
        ret = {}
        for fac in self.patch.allFacilities:
            for ward in fac.getWards():
                # for now let's drop any "HOME" wards
                if ward.tier == CareTier.HOME:
                    continue
                assert hasattr(ward.iA, 'initialFracColonized'), (('ward %s'
                                                                  ' has no initialFracColonized')
                                                                  % ward._name)
                ret[(fac.abbrev, CareTier.names[ward.tier])] = ward.iA.initialFracColonized
        return ret

    def processPrevData(self):
        pthData = self.monitor.getPthData()
        pthDataName = "cre_prev_%s.mpk"%self.monitor.uniqueID
        pthData.to_msgpack(pthDataName, compress="zlib")

        tauDataName = "tau_data_%s.pkl"%self.monitor.uniqueID
        tauDict = getTauDict(self.patch)
        with open(tauDataName, "wb") as f:
            pickle.dump(tauDict, f, 2)

        expectedName = "expected_data_%s.pkl"%self.monitor.uniqueID
        with open(expectedName, "wb") as f:
            pickle.dump(self.expectedPrevalence, f, 2)

        taumod.updateInfo(self.monitor.uniqueID, pthDataName, tauDataName, expectedName)
        self.nextDate, tauDict = taumod.getNewTauDict(self.nextDate)
        overrideTaus(self.patch, tauDict)

    def daily(self, timeNow):
        if timeNow < self.nextDate:
            return self.nextDate

        self.processPrevData()

        return self.nextDate

    def createCallbackFn(self):
        def fn(timeNow, mon):
            return self.daily(timeNow)
        return fn
        