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

###############################
# A good command line for the resulting .dot file is:
#
#   neato -Gspline=true -Nfontsize=5 -Nheight=0.1 -Efontsize=5 -n2 -Tsvg -ograph.svg graph.dot
#
###############################

import sys
import os.path

import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools
import math
import types

import numpy as np
if __name__ == "__main__":
    import matplotlib.pyplot as plt

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))
import pyrheautils
import schemautils
from pyrhea import checkInputFileSchema
import yaml


class Graph(object):

    def __init__(self, fileName, title=None):
        self.ofile = open(fileName, 'w')
        self.title = title

    def __enter__(self):
        if self.title is None:
            self.ofile.write('digraph myGraph {\n')
        else:
            self.ofile.write('digraph myGraph {\n overlap=false; label="%s";' % self.title)
        return self

    def __exit__(self, tp, value, traceback):
        self.ofile.write('}\n')
        self.ofile.close()

    def addNode(self, nodeName, attrDict={}):
        attrStr = 'label=%s' % nodeName
        for k1, k2 in [('label', 'label'),
                       ('width', 'width'),
                       ('color', 'color'),
                       ('shape', 'shape'),
                       ('style', 'style'),
                       ('height', 'height'),
                       ('fillcolor', 'fillcolor'),
                       ('pos', 'pos')]:
            if k1 in attrDict:
                attrStr += ', %s=%s' % (k2, attrDict[k1])
#         attrStr = attrStr[1:]  # drop leading comma
        self.ofile.write('%s [%s];\n' % (nodeName, attrStr))

    def addEdge(self, fromNodeName, toNodeName, attrDict={}):
        attrStr = ''
        if 'weight' in attrDict:
            attrDict['_layoutWt'] = int(math.ceil(10.0*attrDict['weight']))
            attrDict['_penWt'] = max(0.3, 3.0*math.log(attrDict['weight']))
        for k1, k2 in [('_penWt', 'penwidth'),
                       ('_layoutWt', 'weight'),
                       ('label', 'label'),
                       ('color', 'color'),
                       ('dir', 'dir')]:
            if k1 in attrDict:
                attrStr += ', %s=%s' % (k2, attrDict[k1])
        attrStr = attrStr[1:]  # drop leading comma
        self.ofile.write('%s -> %s [%s];\n' % (fromNodeName, toNodeName, attrStr))


catchLocX = []
catchLocY = []
catchSize = []
catchColor = []
worldRadius = None
mapCtr = None
mapXVec = None
mapYVec = None


def facToNodeAttrs(facRec):
    global catchLocX, catchLocY, catchSz, catchColor
    attrDict = {'label': facRec['abbrev']}

    if 'nBeds' in facRec:
        try:
            if isinstance(facRec['nBeds'], (types.IntType, types.FloatType)):
                szFac = facRec['nBeds']
            else:
                szFac = float(facRec['nBeds']['value'])
        except ValueError:
            print '%s has a bad nBeds value' % facRec['abbrev']
            szFac = None
    elif 'meanPop' in facRec:
        if isinstance(facRec['meanPop'], (types.IntType, types.FloatType)):
            szFac = float(facRec['meanPop'])
        else:
            try:
                szFac = float(facRec['meanPop']['value'])
            except ValueError:
                print ('%s meanPop invalid value <%s>' %
                       (facRec['abbrev'], facRec['meanPop']['value']))
                szFac = None
    else:
        print '%s has no nBeds or meanPop' % facRec['abbrev']
        szFac = None

    if facRec['category'] == 'HOSPITAL':
        attrDict['shape'] = 'box'
        catchColor.append('blue')
    elif facRec['category'] == 'NURSINGHOME':
        attrDict['shape'] = 'ellipse'
        catchColor.append('red')
    elif facRec['category'] == 'LTAC':
        attrDict['shape'] = 'triangle'
        catchColor.append('blue')
    elif facRec['category'] == 'COMMUNITY':
        attrDict['shape'] = 'star'
        szFac *= 0.01
        catchColor.append('green')
    else:
        attrDict['shape'] = 'octagon'
        catchColor.append('violet')

    if szFac is None:
        attrDict['color'] = 'blue'
        attrDict['style'] = 'dashed'
    else:
        attrDict['width'] = 0.01 * float(szFac)

    x, y = mapProject(facRec['latitude'], facRec['longitude'],
                      worldRadius, mapCtr, mapXVec, mapYVec)
    # catchLocX.append(facRec['longitude'])
    # catchLocY.append(facRec['latitude'])
    catchLocX.append(x)
    catchLocY.append(y)
    catchSize.append(szFac)

    attrDict['pos'] = '"%f,%f?"' % (x, y)

    return attrDict


def importTransferTable(fname, abbrevPrefix='To_', keyLabel=''):
    with open(fname, 'r') as f:
        keys, transferRecs = csv_tools.parseCSV(f)  # @UnusedVariable

    transferDict = {}
    for r in transferRecs:
        newR = {}
        for k, v in r.items():
            if abbrevPrefix:
                if k.startswith(abbrevPrefix):
                    newR[k[len(abbrevPrefix):]] = v
            else:
                if len(k):
                    newR[k] = v
        transferDict[r[keyLabel]] = newR
    return transferDict


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


def mapTo3D(latitude, longitude, worldRadius):
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


def mapProject(lat, lon, r, ctr, xVec, yVec):
    """Returns (x, y) where x and y are the projected map coordinates"""
    sep = vecMinus(mapTo3D(lat, lon, r), ctr)
#     print 'sep: %s' % str(sep)
    return (dot(sep, xVec), dot(sep, yVec))


def parseFacilityData(fnameOrNameList):
    if not isinstance(fnameOrNameList, types.ListType):
        fnameOrNameList = [fnameOrNameList]
    facDict = {}
    for fn in fnameOrNameList:
        junkKeys, facRecs = yaml_tools.parse_all(fn)  # @UnusedVariable
        for r in facRecs:
            assert r['abbrev'] not in facDict, 'Redundant definitions for %s' % r['abbrev']
            facDict[r['abbrev']] = r
    return facDict


def initializeMapCoordinates(facRecs):
    global worldRadius, mapCtr, mapXVec, mapYVec

    # Scale this to get good node spacing
    worldRadius = 250000.0

    # Find the center vector of all the data- this will map to (0,0)
    ctr = (0.0, 0.0, 0.0)
    n = 0
    for rec in facRecs:
        ctr = vecPlus(ctr, mapTo3D(rec['latitude'], rec['longitude'], worldRadius))
        n += 1
    mapCtr = scale(ctr, 1.0/float(n))
    ctrLatRad = math.asin(mapCtr[2] / worldRadius)
    ctrLonRad = math.atan2(mapCtr[0], mapCtr[1])
    print 'Map Center: %s' % str(mapCtr)
    print 'Map Center latitude and longitude: %s %s' % (180.*ctrLatRad/math.pi,
                                                        180.*ctrLonRad/math.pi)

    # X and Y vectors are normalized vectors in the longitude and latitude directions
    mapXVec = (math.cos(ctrLatRad) * math.cos(ctrLonRad),
               -math.cos(ctrLatRad) * math.sin(ctrLonRad),
               0.0)
    mapYVec = (-math.sin(ctrLatRad) * math.sin(ctrLonRad),
               -math.sin(ctrLatRad) * math.cos(ctrLonRad),
               math.cos(ctrLatRad))
    print 'mapXVec: %s' % str(mapXVec)
    print 'mapYVec: %s' % str(mapYVec)


def writeDotGraph(fname, title, facDict, transferDict, inclusionSet=None):
    initializeMapCoordinates(facDict.values())
    facDict = {k.lower(): v for k, v in facDict.items()}  # convert keys to lower case
    transOutDict = {}
    transInDict = {}
    invertDict = {}
    for src, r in transferDict.items():
        src = src.lower()
        totTransfersOut = float(sum([v for v in r.values()]))
        if totTransfersOut == 0.0:
            print '%s has no outgoing transfers' % src
        transOutDict[src] = totTransfersOut
        for dst, v in r.items():
            dst = dst.lower()
            if dst not in transInDict:
                transInDict[dst] = 0.0
            transInDict[dst] += float(v)
            if dst not in invertDict:
                invertDict[dst] = {}
            invertDict[dst][src] = v

    minWtCutoff = 0.1
    minCntCutoff = 2
    wtScale = 0.01

    createdSet = set()
    with Graph(fname, title=title) as g:
        for src in facDict.keys():
            src = src.lower()
#             if (facDict[src]['longitude'] < -87.6 or facDict[src]['longitude'] > -87.2
#                     or facDict[src]['latitude'] < 41.41 or facDict[src]['latitude'] > 41.61):
#                 continue
            if src not in createdSet:
                g.addNode(src, facToNodeAttrs(facDict[src]))
                createdSet.add(src)

        for src, r in transferDict.items():
            src = src.lower()
            if (src not in facDict 
#                 or
#                 facDict[src]['longitude'] < -87.6 or facDict[src]['longitude'] > -87.2
#                     or facDict[src]['latitude'] < 41.41 or facDict[src]['latitude'] > 41.61
                    ):
                continue
            for dst, v in r.items():
                dst = dst.lower()
                if (dst not in facDict 
#                     or
#                     facDict[dst]['longitude'] < -87.6 or facDict[dst]['longitude'] > -87.2
#                         or facDict[dst]['latitude'] < 41.41 or facDict[dst]['latitude'] > 41.61
                        ):
                    continue
                if transOutDict[src] > 0.0 and transInDict[dst] > 0.0:
                    if src in invertDict and dst in invertDict[src]:
                        reverseV = invertDict[src][dst]
                    else:
                        reverseV = 0.0
                    visWt = float(v)/min(transOutDict[src], transInDict[dst])
                    if (visWt >= minWtCutoff
                            and src in facDict and dst in facDict
                            and ((inclusionSet is None)
                                 or ((facDict[src]['category'], facDict[dst]['category'])
                                     in inclusionSet))):
                        if (src, dst) not in createdSet and src != dst:
                            if v >= reverseV:
                                symV = reverseV
                                asymV = v - reverseV
                                if symV >= minCntCutoff:
                                    g.addEdge(src, dst, {'weight': wtScale * symV, 'dir': 'both',
                                                         'label': symV, 'color': 'black'})
                                if asymV >= minCntCutoff:
                                    g.addEdge(src, dst, {'weight': wtScale * asymV, 'color': 'red',
                                                         'label': asymV})
                            else:
                                symV = v
                                asymV = reverseV - v
                                if symV >= minCntCutoff:
                                    g.addEdge(src, dst, {'weight': wtScale * symV, 'dir': 'both',
                                                         'label': symV, 'color': 'black'})
                                if asymV >= minCntCutoff:
                                    g.addEdge(dst, src, {'weight': wtScale * asymV, 'color': 'red',
                                                         'label': asymV})
                            createdSet.add((src, dst))
                            createdSet.add((dst, src))

    print 'wrote %s' % fname
    oname = '%s.svg' % os.path.splitext(os.path.basename(fname))[0]
    print 'A good post-processing command might be:'
    print ('neato -Gspline=true -Nfontsize=5 -Nheight=0.1 -Efontsize=5'
           ' -n2 -Tsvg -o%s %s' % (oname, fname))
    return transInDict, transOutDict


def addFacilityTypeSynonyms(inclusionSet):
    inclusionSet = set(inclusionSet)  # in case we get passed a list
    synDict = {'NURSINGHOME': ['SNF', 'VSNF'],
               'LTAC': ['LTACH'],
               'HOSPITAL': []}
    newSet = inclusionSet.copy()
    for tpl in inclusionSet:
        fm, to = tpl
        for newFm in synDict[fm] + [fm]:
            for newTo in synDict[to] + [to]:
                newSet.add((newFm, newTo))
    return newSet

def main():
    SCHEMA_DIR = '../schemata'
    INPUT_SCHEMA = 'rhea_input_schema.yaml'

    runDesc = sys.argv[1]
    schemautils.setSchemaBasePath(SCHEMA_DIR)
    inputDict = checkInputFileSchema(sys.argv[1],
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    pyrheautils.prepPathTranslations(inputDict)
    facDir = pyrheautils.pathTranslate('$(MODELDIR)/facilityfacts')
    facDict = parseFacilityData(facDir)
    modelDir = inputDict['modelDir']
    
    
    
    #     directTransferDict = importTransferTable('transfer_matrix_direct_normalized.csv')
    #     readmitTransferDict = importTransferTable('transfer_matrix_readmit_normalized.csv')
    #     directTransferDict = importTransferTable('Transferring_matrix_abbrev_9-18-2014_with-update-10-2-2014_copy_RHEA_Direct_fix_HSOU+silos+SCLE.csv',
    #                                              abbrevPrefix=None)
    #    directTransferDict = importTransferTable('/home/welling/git/pyrhea/models/ChicagoLand/'
    #                                             'Matrices_LOS_08092016_Transfer3day.csv',
    #                                             abbrevPrefix='To_', keyLabel='UNIQUE_ID')

    transferYaml = pyrheautils.pathTranslate('$(MODELDIR)/transfer_counts.yaml')
    with open(transferYaml) as f:
        directTransferDict = yaml.load(f)


    readmitTransferDict = {}
    transferDict = {}
    for k, rec in directTransferDict.items():
        transferDict[k] = rec.copy()
    for k, r2 in readmitTransferDict.items():
        r1 = transferDict[k]
        for k2, v in r2.items():
            r1[k2] += v

    initializeMapCoordinates(facDict.values())

    inclusionSet = [('NURSINGHOME', 'HOSPITAL'),
                    ('NURSINGHOME', 'LTAC'),
                    ('LTAC', 'NURSINGHOME'),
                    ('HOSPITAL', 'NURSINGHOME'),
                    ('HOSPITAL', 'HOSPITAL'),
                    ('NURSINGHOME', 'NURSINGHOME'),
                    ('LTAC', 'LTAC'),
                    ('LTAC', 'HOSPITAL'),
                    ('HOSPITAL', 'LTAC')
                    ]
    inclusionSet = addFacilityTypeSynonyms(inclusionSet)
    title = "HOSPITAL + LTAC internal direct transfers"
    # title = 'NURSINGHOME internal direct transfers'
    # title = "NURSINGHOME to/from HOSPITAL+LTAC direct + indirect transfers"

    transInDict, transOutDict = writeDotGraph('graph.dot', title,
                                              facDict, transferDict, inclusionSet)

    orderedSources = facDict.keys()[:]
    orderedSources.sort()
    totalIn = 0
    totalOut = 0
    totalInHosp = 0
    totalOutHosp = 0
    totalInNH = 0
    totalOutNH = 0
    totBeds = 0
    totBedsHosp = 0
    totBedsNH = 0
    totBedsHospFail = False
    totBedsFail = False
    totBedsNHFail = False
    for src in orderedSources:
        key = None
        bedCt = None
        print ('%s (%s): %d %d;' %
               (src, facDict[src]['category'], transInDict[src.lower()],
                transOutDict[src.lower()])),
        totalIn += transInDict[src.lower()]
        totalOut += transOutDict[src.lower()]
        if facDict[src]['category'] in ['NURSINGHOME', 'SNF', 'VSNF']:
            totalInNH += transInDict[src.lower()]
            totalOutNH += transOutDict[src.lower()]
            if 'nBeds' in facDict[src]:
                key = 'nBeds'
            elif 'meanPop' in facDict[src]:
                key = 'meanPop'
            else:
                print ' no nBeds and no meanPop'
            if key:
                if isinstance(facDict[src][key], types.IntType):
                    bedCt = facDict[src][key]
                else:
                    try:
                        bedCt = int(facDict[src][key]['value'])
                    except ValueError:
                        pass
            if bedCt:
                totBeds += bedCt
                totBedsNH += bedCt
                print '%s beds' % bedCt
            else:
                print ' no nBeds and no meanPop'
                totBedsNHFail = True
        else:
            totalInHosp += transInDict[src.lower()]
            totalOutHosp += transOutDict[src.lower()]
            if 'meanPop' in facDict[src]:
                key = 'meanPop'
            elif 'nBeds' in facDict[src]:
                key = 'nBeds'
            else:
                print ' no nBeds and no meanPop'
            if key:
                if isinstance(facDict[src][key], types.FloatType):
                    bedCt = facDict[src][key]
                else:
                    try:
                        bedCt = float(facDict[src][key]['value'])
                    except ValueError:
                        pass
            if bedCt:
                totBeds += bedCt
                totBedsHosp += bedCt
                print '%s beds' % bedCt
            else:
                print ' no nBeds and no meanPop'
                totBedsHospFail = True


    print 'Total incoming transfers: %s (%s hospitals, %s NH)' % (totalIn, totalInHosp, totalInNH)
    print 'Total outgoing transfers: %s (%s hospitals, %s NH)' % (totalOut, totalOutHosp,
                                                                  totalOutNH)
    if totBedsHospFail:
        if totBedsNHFail:
            print 'Incomplete data; beds were not counted'
        else:
            print 'Tot NH beds: %s (others not counted)' % totBedsNH
    elif totBedsNHFail:
        print 'Tot HOSP beds: %s (others not counted)' % totBedsHosp
    else:
        print 'Total beds: %s (%s hospital beds, %s NH beds)' % (totBeds, totBedsHosp, totBedsNH)

    # fig, ax = plt.subplots()
    # ax.scatter(catchLocX, catchLocY, s=catchSize, c=catchColor, alpha=0.3)
    # ax.grid(True)
    # fig.tight_layout()
    # plt.show()

if __name__ == "__main__":
    main()
