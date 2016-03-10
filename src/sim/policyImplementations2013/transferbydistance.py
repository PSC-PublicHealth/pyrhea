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

_rhea_svn_id_ = "$Id$"

import logging

import os.path
import math
import pyrheautils
from facilitybase import TransferDestinationPolicy as BaseTransferDestinationPolicy
from facilitybase import CareTier, tierToQueueMap

_validator = None
_constants_values = 'transferbydistance_constants.yaml'
_constants_schema = 'transferbydistance_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)


def longitudeLatitudeSep(lon1, lat1, lon2, lat2):
    "Inputs are in floating point degrees.  Returns separation in kilometers"
    scale = lat1r = lat2r = lon1r = lon2r = a = b = 0.0
    try:
        scale = math.pi / 180.
        lat1r = lat1*scale
        lon1r = lon1*scale
        lat2r = lat2*scale
        lon2r = lon2*scale

        a = math.sin(lat1r)*math.sin(lat2r)
        b = math.cos(lat1r)*math.cos(lat2r) * math.cos(lon2r-lon1r)
        apb = a + b
        if apb > 1.0:
            apb = 1.0  # avoid rounding error
        if apb < -1.0:
            apb = -1.0  # avoid rounding error
        c = math.acos(apb)
    except Exception, e:
        print ('longitudeLatitudeSep: <%s> <%s> <%s> <%s> -> %s %s -> %s %s'
               % (lon1, lat1, lon2, lat2, lat1r, lat2r, a, b))
        raise e
    R = 6378.  # radius of earth in km; in miles it's 3963.189
    return R*c


class MinDistanceCore(object):
    """This is where we put things that are best shared across all instances"""
    __metaclass__ = pyrheautils.SingletonMetaClass

    def __init__(self):
        pass


class MinDistanceTransferDestinationPolicy(BaseTransferDestinationPolicy):
    def __init__(self, patch):
        super(MinDistanceTransferDestinationPolicy, self).__init__(patch)
        self.core = MinDistanceCore()
        self.cache = {}

    def buildCacheEntry(self, oldFacility, newTier):
        srcLon, srcLat = oldFacility.coords
        queueClass = tierToQueueMap[newTier]
        tplList = []
        for info, addr in self.patch.serviceLookup(queueClass.__name__):
            innerInfo, destAbbrev, destCoords = info  # @UnusedVariable
            dstLon, dstLat = destCoords
            sep = longitudeLatitudeSep(srcLon, srcLat, dstLon, dstLat)
            tplList.append((sep, info, addr))
        tplList.sort()
        return [c for a, b, c in tplList]  # @UnusedVariable

    def getOrderedCandidateFacList(self, oldFacility, oldTier, newTier, timeNow):
        if newTier not in self.cache:
            self.cache[newTier] = self.buildCacheEntry(oldFacility, newTier)
        return self.cache[newTier][:]


def getPolicyClasses():
    return [MinDistanceTransferDestinationPolicy]


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(os.path.join(os.path.dirname(__file__),
                                                      _constants_values),
                                         _constants_schema)
