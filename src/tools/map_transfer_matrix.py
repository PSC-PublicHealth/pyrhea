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

import os.path

import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools
import math

import numpy as np
import matplotlib.pyplot as plt


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
                       ('pos', 'pos')]:
            if k1 in attrDict:
                attrStr += ', %s=%s' % (k2, attrDict[k1])
        attrStr = attrStr[1:]  # drop leading comma
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
        szFac = facRec['nBeds']
    elif 'meanPop' in facRec:
        szFac = facRec['meanPop']
    else:
        print '%s has no nBeds or meanPop'
        szFac = None
    if szFac is None:
        attrDict['color'] = 'blue'
        attrDict['style'] = 'dashed'
    else:
        scaledSz = 0.01 * float(szFac)
        attrDict['width'] = scaledSz

    if facRec['category'] == 'HOSPITAL':
        attrDict['shape'] = 'box'
        catchSize.append(facRec['meanPop'])
        catchColor.append('blue')
    elif facRec['category'] == 'NURSINGHOME':
        attrDict['shape'] = 'ellipse'
        catchSize.append(facRec['nBeds'])
        catchColor.append('red')
    else:
        attrDict['shape'] = 'triangle'
        catchSize.append(facRec['meanPop'])
        catchColor.append('blue')

    x, y = mapProject(facRec['latitude'], facRec['longitude'],
                      worldRadius, mapCtr, mapXVec, mapYVec)
    # catchLocX.append(facRec['longitude'])
    # catchLocY.append(facRec['latitude'])
    catchLocX.append(x)
    catchLocY.append(y)

    attrDict['pos'] = '"%f,%f?"' % (x, y)

    return attrDict


def importTransferTable(fname):
    with open(fname, 'r') as f:
        keys, transferRecs = csv_tools.parseCSV(f)  # @UnusedVariable

    transferDict = {}
    for r in transferRecs:
        newR = {}
        for k, v in r.items():
            if k.startswith('To_'):
                newR[k[3:]] = v
        transferDict[r['']] = newR
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


def parseFacilityData(fname):
    junkKeys, facRecs = yaml_tools.parse_all(fname)  # @UnusedVariable
    facDict = {r['abbrev']: r for r in facRecs}
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


def writeDotGraph(fname, title, facDict, transferDict, inclusionSet):
    facDict = {k.lower(): v for k, v in facDict.items()}  # convert keys to lower case
    transOutDict = {}
    transInDict = {}
    invertDict = {}
    for src, r in transferDict.items():
        totTransfersOut = float(sum([v for v in r.values()]))
        if totTransfersOut == 0.0:
            print '%s has no outgoing transfers' % src
        transOutDict[src] = totTransfersOut
        for dst, v in r.items():
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
        for src, r in transferDict.items():
            for dst, v in r.items():
                if transOutDict[src] > 0.0 and transInDict[dst] > 0.0:
                    if src in invertDict and dst in invertDict[src]:
                        reverseV = invertDict[src][dst]
                    else:
                        reverseV = 0.0
                    visWt = float(v)/min(transOutDict[src], transInDict[dst])
                    if (visWt >= minWtCutoff
                            and src in facDict and dst in facDict
                            and ((facDict[src]['category'], facDict[dst]['category'])
                                 in inclusionSet)):
                        if src not in createdSet:
                            g.addNode(src, facToNodeAttrs(facDict[src]))
                            createdSet.add(src)
                        if dst not in createdSet:
                            g.addNode(dst, facToNodeAttrs(facDict[dst]))
                            createdSet.add(dst)
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


def main():
    facDict = parseFacilityData('/home/welling/workspace/pyRHEA/models/OrangeCounty/'
                                'facilityfactsCurrent')

    directTransferDict = importTransferTable('transfer_matrix_direct_normalized.csv')
    readmitTransferDict = importTransferTable('transfer_matrix_readmit_normalized.csv')
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
    # title = "HOSPITAL + LTAC internal direct transfers"
    # title = 'NURSINGHOME internal direct transfers'
    title = "NURSINGHOME to/from HOSPITAL+LTAC direct + indirect transfers"

    transInDict, transOutDict = writeDotGraph('graph.dot', title,
                                              facDict, transferDict, inclusionSet)

    orderedSources = transferDict.keys()[:]
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
    for src in orderedSources:
        print '%s (%s): %d %d' % (src, facDict[src]['category'], transInDict[src],
                                  transOutDict[src])
        totalIn += transInDict[src]
        totalOut += transOutDict[src]
        if facDict[src]['category'] == 'NURSINGHOME':
            totalInNH += transInDict[src]
            totalOutNH += transOutDict[src]
            if 'nBeds' in facDict[src]:
                totBeds += facDict[src]['nBeds']
                totBedsNH += facDict[src]['nBeds']
            else:
                print '%s has no nBeds' % src
        else:
            totalInHosp += transInDict[src]
            totalOutHosp += transOutDict[src]
            if 'meanPop' in facDict[src]:
                totBeds += facDict[src]['meanPop']
                totBedsHosp += facDict[src]['meanPop']
            else:
                print '%s has no meanPop' % src

    print 'Total incoming transfers: %s (%s hospitals, %s NH)' % (totalIn, totalInHosp, totalInNH)
    print 'Total outgoing transfers: %s (%s hospitals, %s NH)' % (totalOut, totalOutHosp,
                                                                  totalOutNH)
    print 'Total beds: %s (%s hospitals, %s NH)' % (totBeds, totBedsHosp, totBedsNH)

    # fig, ax = plt.subplots()
    # ax.scatter(catchLocX, catchLocY, s=catchSize, c=catchColor, alpha=0.3)
    # ax.grid(True)
    # fig.tight_layout()
    # plt.show()

if __name__ == "__main__":
    main()
