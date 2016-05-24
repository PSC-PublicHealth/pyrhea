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

_rhea_svn_id_ = "$Id: map_facilities.py 133 2016-03-07 20:12:44Z welling $"

import os.path
import math
import optparse
import logging
import signal
import yaml
import numpy as np
import matplotlib.pyplot as plt
from map_transfer_matrix import parseFacilityData
from notes_plotter import importNotes
from pyrheautils import loadModulesFromDir
import schemautils
from descartes import PolygonPatch
import json
import __builtin__
# http://stackoverflow.com/questions/18229628/python-profiling-using-line-profiler-clever-way-to-remove-profile-statements
try:
    profile = __builtin__.profile
except AttributeError:
    # No line profiler, provide a pass-through version
    def profile(func): return func


schemaDir = '../schemata'
inputSchema = 'rhea_input_schema.yaml'

logger = None


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
    scale = math.sqrt(dot(v, v))
    assert scale > 0.0, "Tried to normalize a zero-length vector"
    return (v[0]/scale, v[1]/scale, v[2]/scale)


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


def project(v, xVec, yVec, center):
    """xVec and yVec are presumed to be normalized"""
    vCtr = vecMinus(v, center)
    vx = scale(xVec, dot(vCtr, xVec))
    vy = scale(yVec, dot(vCtr, yVec))
    return vecPlus(vx, vy)


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
def projPoly(poly, mapProj):
    newCoords = []
    lonSum = 0.0
    latSum = 0.0
    count = 0
    for l1 in poly['coordinates']:
        nl1 = []
        for l2 in l1:
            lon, lat = mapProj.project(l2[0], l2[1])
            nl1.append([lon, lat])
            lonSum += lon
            latSum += lat
            count += 1
        newCoords.append(nl1)
    newPoly = poly.copy()
    newPoly['coordinates'] = newCoords
    return newPoly, (lonSum/float(count), latSum/float(count))


def categoryToColor(cat):
    return {'HOSPITAL': 'red', 'LTAC': 'pink', 'NURSINGHOME': 'blue', 'COMMUNITY': 'pink'}[cat]


def buildLclWorkVec(orderedFacList, facDict, implementationDir):
    implDict = {}
    for mod in loadModulesFromDir(implementationDir, ['category', 'estimateWork']):
        implDict[mod.category] = mod
    wtDict = {}
    for abbrev, descr in facDict.items():
        wtDict[abbrev] = implDict[descr['category']].estimateWork(descr)
    lclWorkVec = np.array([wtDict[k] for k in orderedFacList])
    return lclWorkVec


def buildTransferMatrix(orderedFacList, facDict, notesPath):
    idx = {nm: i for i, nm in enumerate(orderedFacList)}
    notes = importNotes(notesPath)
    transferMatrix = np.zeros((len(orderedFacList), len(orderedFacList)))
    for fmAbbrev in orderedFacList:
        fmKey = '%s_%s' % (facDict[fmAbbrev]['category'].upper(), fmAbbrev)
        if fmKey in notes:
            for k, v in notes[fmKey].items():
                if k.endswith('_transfer'):
                    toAbbrev = k[:-len('_transfer')]
                    transferMatrix[idx[fmAbbrev], idx[toAbbrev]] += v
        else:
            print 'no notes for %s' % fmKey
    return transferMatrix


def checkInputFileSchema(fname, schemaFname):
    if logger is None:
        myLogger = logging.getLogger(__name__)
    else:
        myLogger = logger
    try:
        with open(fname, 'rU') as f:
            inputJSON = yaml.safe_load(f)
        validator = schemautils.getValidator(os.path.join(os.path.dirname(__file__), schemaFname))
        nErrors = sum([1 for e in validator.iter_errors(inputJSON)])  # @UnusedVariable
        if nErrors:
            myLogger.error('Input file violates schema:')
            for e in validator.iter_errors(inputJSON):
                myLogger.error('Schema violation: %s: %s' %
                               (' '.join([str(word) for word in e.path]), e.message))
            raise RuntimeError('Input file violates its schema')
        else:
            return inputJSON
    except Exception, e:
        myLogger.error('Error checking input against its schema: %s' % e)
        raise RuntimeError('Input file violates its schema')


def main():
    global logger

    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)

    logging.basicConfig(format="%%(levelname)s:%%(name)s:%%(message)s")

    parser = optparse.OptionParser(usage="""
    %prog [-v][-d][-n notes1.pkl [-n notes2.pkl [...]] input.yaml
    """)
    parser.add_option("-v", "--verbose", action="store_true",
                      help="verbose output")
    parser.add_option("-d", "--debug", action="store_true",
                      help="debugging output")
    parser.add_option("-n", "--notes", action="append",
                      help="pickled notes file to provide transfer data")

    opts, args = parser.parse_args()

    if opts.debug:
        logging.getLogger('root').setLevel(logging.DEBUG)
    elif opts.verbose:
        logging.getLogger('root').setLevel(logging.INFO)
    else:
        logging.getLogger('root').setLevel(logging.WARNING)
    logger = logging.getLogger(__name__)

    if opts.notes:
        notesFileList = opts.notes
    else:
        notesFileList = []

    if len(args) == 1:
        inputDict = checkInputFileSchema(args[0], os.path.join(schemaDir, inputSchema))
    else:
        parser.error("A YAML-format file specifying prototype model parameters must be specified.")
    parser.destroy()

    geoData = '/home/welling/geo/tiger/tigr_2010_06.json'

    thisDir = os.path.dirname(__file__)
    modelDir = os.path.join(thisDir, '../../models/OrangeCounty2013/')
    facDirList = inputDict['facilityDirs']
    defaultNotesPath = os.path.join(thisDir, '../sim/notes.pkl')
    implementationDir = inputDict['facilityImplementationDir']

    facDict = parseFacilityData(facDirList)

    schemautils.setSchemaBasePath(schemaDir)

    orderedFacList = facDict.keys()[:]
    orderedFacList.sort()

    transferMatrix = np.zeros((len(orderedFacList), len(orderedFacList)))
    for notesPath in notesFileList:
        transferMatrix += buildTransferMatrix(orderedFacList, facDict, notesPath)
#         diagPairs = zip(orderedFacList, transferMatrix.diagonal())
#         for abbrev, val in diagPairs:
#             assert val == 0 or facDict[abbrev]['category'].upper() in ['HOSPITAL'], 'fail %s %s' % (abbrev, val)

    lclWorkVec = buildLclWorkVec(orderedFacList, facDict, implementationDir)
    print transferMatrix
    print lclWorkVec

# 
#     mapProj = MapProjection((sum([r['longitude'] for r in facDict.values()]) / len(facDict),
#                              sum([r['latitude'] for r in facDict.values()]) / len(facDict)),
#                             100.0)
# 
#     BLUE = '#6699cc'
#     fig = plt.figure()
#     ax = fig.gca()
# 
#     with open(geoData, 'r') as f:
#         json_data = json.load(f)
#     for feature in json_data['features']:
#         if feature['properties']['COUNTY'] == '059':
#         #if True:
#             # print feature['properties']['TRACT']
#             # print feature['properties']
#             poly = feature['geometry']
#             if poly['type'] == 'Polygon':
#                 poly, (ctrLon, ctrLat) = projPoly(poly, mapProj)
#                 ax.add_patch(PolygonPatch(poly, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2))
#                 ax.annotate(feature['properties']['NAME'], (ctrLon, ctrLat), (ctrLon, ctrLat),
#                             fontsize=5, ha='center', va='center')
#             else:
#                 pass
# 
#     for cg, mrk, sz in {('HOSPITAL', '*', 15), ('LTAC', '+', 15), ('NURSINGHOME', 'o', 10)}:
#         coordList = [(mapProj.project(rec['longitude'], rec['latitude']),
#                      rec['abbrev'])
#                      for rec in facDict.values() if rec['category'] == cg]
#         coordList = [(x, y, abbrev) for (x, y), abbrev in coordList]
#         ax.scatter([x for x, y, abbrev in coordList], [y for x, y, abbrev in coordList],
#                    c=categoryToColor(cg), marker=mrk)
#         for x, y, abbrev in coordList:
#             ax.annotate(abbrev, (x, y), (x, y), ha='left', va='bottom', fontsize=sz)
# 
# 
# #     for x, y, clr, abbrev in coordList:
# #         ax.annotate(abbrev, xy=(x, y), xytext=(x, y))
# #     for abbrev, offset in indexDict.items():
# #         xy = (xVals[offset], yVals[offset])
# #         xytext = (xVals[offset]+0.0125, yVals[offset]+0.0125)
# #         scatterAx.annotate(abbrev, xy=xy, xytext=xytext)
# 
#     ax.axis('scaled')
#     plt.show()


if __name__ == "__main__":
    main()
