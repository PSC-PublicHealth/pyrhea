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

import math
import logging
import json
import re
import types
import matplotlib.pyplot as plt
from descartes import PolygonPatch

# http://stackoverflow.com/questions/18229628/python-profiling-using-line-profiler-clever-way-to-remove-profile-statements
try:
    profile = __builtins__.profile
except AttributeError:
    # No line profiler, provide a pass-through version
    def profile(func):
        """Dummy profile in case it is not present in __builtin__"""
        return func

logger = logging.getLogger(__name__)


def dot(v1, v2):
    return (v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2])


def cross(v1, v2):
    v1x, v1y, v1z = v1
    v2x, v2y, v2z = v2
    v3x = v1y * v2z - v1z * v2y
    v3y = v1z * v2x - v1x * v2z
    v3z = v1x * v2y - v1y * v2x
    return (v3x, v3y, v3z)


def scale(v, fac):
    return (fac*v[0], fac*v[1], fac*v[2])


def normalize(v):
    scl = math.sqrt(dot(v, v))
    assert scl > 0.0, "Tried to normalize a zero-length vector"
    return (v[0]/scl, v[1]/scl, v[2]/scl)


def mapTo3D(longitude, latitude, worldRadius):
    """Project from latitude and longitude to (x,y,z)"""
    lonRad = math.pi * longitude / 180.0
    latRad = math.pi * latitude / 180.0
    z = worldRadius * math.sin(latRad)
    x = worldRadius * math.cos(latRad) * math.sin(lonRad)
    y = worldRadius * math.cos(latRad) * math.cos(lonRad)
    return (x, y, z)


def vecPlus(v1, v2):
    return (v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2])


def vecMinus(v1, v2):
    return (v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])


class MapProjection(object):
    def __init__(self, mapCtr, worldRadius):
        """
        mapCtr should be a (lon, lat) tuple.  worldRadius scales the final projection.
        """
        self.worldRadius = worldRadius
        self.mapCtr3D = mapTo3D(mapCtr[0], mapCtr[1], worldRadius)

        ctrLatRad = math.asin(self.mapCtr3D[2] / worldRadius)
        ctrLonRad = math.atan2(self.mapCtr3D[0], self.mapCtr3D[1])
        print 'Input map center longitude and latitude: %s' % str(mapCtr)
        print 'Map Center: %s' % str(self.mapCtr3D)
        print 'Map Center longitude and latitude: %s %s' % (180.*ctrLonRad/math.pi,
                                                            180.*ctrLatRad/math.pi)

        # X and Y vectors are normalized vectors in the longitude and latitude directions
        self.mapXVec = (math.cos(ctrLatRad) * math.cos(ctrLonRad),
                        -math.cos(ctrLatRad) * math.sin(ctrLonRad),
                        0.0)
        self.mapYVec = (-math.sin(ctrLatRad) * math.sin(ctrLonRad),
                        -math.sin(ctrLatRad) * math.cos(ctrLonRad),
                        math.cos(ctrLatRad))
        print 'mapXVec: %s' % str(self.mapXVec)
        print 'mapYVec: %s' % str(self.mapYVec)

    @profile
    def project(self, lon, lat):
        sep = vecMinus(mapTo3D(lon, lat, self.worldRadius), self.mapCtr3D)
        x = dot(sep, self.mapXVec)
        y = dot(sep, self.mapYVec)
        return (x, y)

    @profile
    def projCoordList(self, coordL):
        newCoords = []
        lonSum = 0.0
        latSum = 0.0
        count = 0
        for l1 in coordL:
            nl1 = []
            for l2 in l1:
                lon, lat = self.project(l2[0], l2[1])
                nl1.append([lon, lat])
                lonSum += lon
                latSum += lat
                count += 1
            newCoords.append(nl1)
        return newCoords, lonSum/float(count), latSum/float(count)

    @profile
    def projPoly(self, poly):
        if poly['type'] == 'Polygon':
            newCoords, lonMean, latMean = self.projCoordList(poly['coordinates'])
        elif poly['type'] == 'MultiPolygon':
            newCoords = []
            lonVals = []
            latVals = []
            for coordList in poly['coordinates']:
                nC, lonM, latM = self.projCoordList(coordList)
                newCoords.append(nC)
                lonVals.append(lonM)
                latVals.append(latM)
            lonMean = sum(lonVals)/float(len(lonVals))
            latMean = sum(latVals)/float(len(latVals))
        else:
            logger.fatal('cannot handle poly type %s', poly['type'])
        newPoly = poly.copy()
        newPoly['coordinates'] = newCoords
        return newPoly, (lonMean, latMean)


class Map(object):
    mrkDefaultSz = 10
#     mrkSzDict = {'*': 60, '+': 60, 'o': 60}
#     mrkDefaultSz = 30
    tractLabelDefaultSz = 10

    def __init__(self, geoDataPathList, stateCodeStr, countyCodeStr, ctrLon, ctrLat,
                 annotate=True,
                 annotateTracts=True,
                 nameMap=None,
                 regexCodes=False,
                 mrkSzDict=None):
        """
        Annotation of tracts requires that both annotate and annotateTracts be True.
        Annotation of markers requires only that annotate be True.
        """
        if isinstance(geoDataPathList, types.StringTypes):
            geoDataPathList = [geoDataPathList]
        if nameMap is None:
            nameMap = {'STATE': 'STATE', 'COUNTY': 'COUNTY', 'TRACT': 'TRACT', 'NAME': 'NAME',
                       'GEOKEY': 'GEO_ID'}
        self.nameMap = nameMap.copy()
        if mrkSzDict is None:
            mrkSzDict = {'*': 15, '+': 15, 'o': 10}
        self.mrkSzDict = mrkSzDict.copy()
        self.tractPolyDict = {}
        self.tractPropertyDict = {}
        self.mrkCoordDict = {}
        stateKey = self.nameMap['STATE']
        countyKey = self.nameMap['COUNTY']
        geoKey = self.nameMap['GEOKEY']
        if regexCodes:
            stateRe = re.compile(stateCodeStr)
            countyRe = re.compile(countyCodeStr)
        else:
            stateRe = None
            countyRe = None
        for geoDataPath in geoDataPathList:
            with open(geoDataPath, 'r') as f:
                json_data = json.load(f)
            if regexCodes:
                for feature in json_data['features']:
                    if (stateRe.match(feature['properties'][stateKey]) and
                            countyRe.match(feature['properties'][countyKey])):
                        geoID = feature['properties'][geoKey]
                        self.tractPolyDict[geoID] = feature['geometry']
                        self.tractPropertyDict[geoID] = feature['properties']
            else:
                for feature in json_data['features']:
                    if (feature['properties'][stateKey] == stateCodeStr and
                            feature['properties'][countyKey] == countyCodeStr):
                        geoID = feature['properties'][geoKey]
                        if geoID in self.tractPolyDict:
                            logger.fatal('Tract ID collision for ID %s', geoID)
                        self.tractPolyDict[geoID] = feature['geometry']
                        self.tractPropertyDict[geoID] = feature['properties']

        if not (ctrLon and ctrLat):
            if self.tractPolyDict:
                ctrLon = 0.0
                ctrLat = 0.0
                count = 0
                coordListList = []
                for poly in self.tractPolyDict.values():
                    if poly['type'] == 'Polygon':
                        coordListList.append(poly['coordinates'][0])
                    elif poly['type'] == 'MultiPolygon':
                        coordListList.extend([elt[0] for elt in poly['coordinates']])
                    else:
                        logger.fatal('cannot handle poly type %s', poly['type'])
                    for coordList in coordListList:
                        sumLon = sum([lon for lon, lat in coordList])  # UnusedVariable
                        sumLat = sum([lat for lon, lat in coordList])  # UnusedVariable
                        lonMean = float(sumLon) / len(coordList)
                        latMean = float(sumLat) / len(coordList)
                    ctrLat += latMean
                    ctrLon += lonMean
                    count += 1
                ctrLon /= count
                ctrLat /= count
            else:
                logger.fatal('No projection center provided and I cannot calculate one')

        self.mapProj = MapProjection((ctrLon, ctrLat), 100.0)
        self.fig = plt.figure()
        self.ax = self.fig.gca()
        self.annotateMarkers = annotate
        self.annotateTracts = annotateTracts and annotate

    def tractIDList(self):
        return self.tractPolyDict.keys()

    def getTractProperties(self, geoID):
        return self.tractPropertyDict[geoID]

    def plotTract(self, geoID, clr):
        poly = self.tractPolyDict[geoID]
        poly, (ctrLon, ctrLat) = self.mapProj.projPoly(poly)
        if poly['type'] == 'Polygon':
            self.ax.add_patch(PolygonPatch(poly, fc=clr, ec=clr, alpha=0.5, zorder=2))
        elif poly['type'] == 'MultiPolygon':
            for cL in poly['coordinates']:
                subPoly = {'type': 'Polygon', 'coordinates': cL}
                self.ax.add_patch(PolygonPatch(subPoly, fc=clr, ec=clr, alpha=0.5, zorder=2))
        else:
            logger.fatal('cannot handle poly type %s', poly['type'])

        if self.annotateTracts:
            self.ax.annotate(self.tractPropertyDict[geoID][self.nameMap['NAME']],
                             (ctrLon, ctrLat), (ctrLon, ctrLat),
                             fontsize=Map.tractLabelDefaultSz,
                             ha='center', va='center')

    def plotMarker(self, lon, lat, mark, label, clr, sz=None):
        if sz is None:
            sz = 1.0
        coords = (self.mapProj.project(lon, lat), label)
        key = (mark, clr)
        if key in self.mrkCoordDict:
            self.mrkCoordDict[key].append((coords, sz))
        else:
            self.mrkCoordDict[key] = [(coords, sz)]
        

    def draw(self, mrkLegendDict=None):
        legendL = []
        artistL = []
        for (mrk, clr), coordList in self.mrkCoordDict.items():
            if mrk in self.mrkSzDict:
                mrkSz = self.mrkSzDict[mrk]
            else:
                mrkSz = Map.mrkDefaultSz
            artist = self.ax.scatter([x for ((x, y), abbrev), sz in coordList],
                                     [y for ((x, y), abbrev), sz in coordList],
                                     c=clr,
                                     s=[(sz * mrkSz) for ((x, y), abbrev), sz in coordList],
                                     marker=mrk,
                                     alpha=0.5)
            if mrkLegendDict:
                artistL.append(artist)
                if mrk in mrkLegendDict:
                    legendL.append(mrkLegendDict[mrk])
                else:
                    legendL.append('???')
            if self.annotateMarkers:
                for ((x, y), abbrev), sz in coordList:
                    self.ax.annotate(abbrev, (x, y), (x, y),
                                     ha='left', va='bottom', fontsize=sz*mrkSz)
        if mrkLegendDict:
            self.ax.legend(artistL, legendL)
        self.ax.axis('scaled')
        plt.show()