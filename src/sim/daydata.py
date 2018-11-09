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
A DayData singleton maintains a cache of dicts, each of which is produced by a generating
function and valid for one 'day'.  When data for a new day is requested, the DayDataGroup
uses the generating function to regenerate the dicts.

Created on Nov 9, 2018

@author: welling
'''

from phacsl.utils.collections.phacollections import SingletonMetaClass

class DayData(object):
    def __init__(self, patch, genFun):
        self.patch = patch
        self.genFun = genFun
        self.curDay = None
        self.curDct = None

    def get(self, timeNow):
        if self.curDay != timeNow:
            self.curDct = self.genFun(self.patch, timeNow)
            self.curDay = timeNow
        return self.curDct

class DayDataGroup(object):
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        self.readyDct = {}
        self.rawDct = {}

    def add(self, key, genFun):
        self.rawDct[key] = genFun

    def get(self, patch, key, timeNow):
        readyKey = (patch.patchId, key)
        if readyKey not in self.readyDct:
            assert key in self.rawDct, 'DayData object %s has not been defined' % key
            self.readyDct[readyKey] = DayData(patch, self.rawDct[key])
        return self.readyDct[readyKey].get(timeNow)


