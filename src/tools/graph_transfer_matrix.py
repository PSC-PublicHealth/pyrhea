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

import csv_tools
import yaml_tools
import sys
import math


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
                       ('height', 'height')]:
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


def facToNodeAttrs(facRec):
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
#         scaledSz = 0.17 * math.sqrt(float(szFac))
#         attrDict['width'] = scaledSz
#         attrDict['height'] = scaledSz
        scaledSz = 0.03 * float(szFac)
        attrDict['width'] = scaledSz

    if facRec['category'] == 'HOSPITAL':
        attrDict['shape'] = 'box'
    elif facRec['category'] == 'NURSINGHOME':
        attrDict['shape'] = 'ellipse'
    else:
        attrDict['shape'] = 'triangle'
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

directTransferDict = importTransferTable('transfer_matrix_direct_normalized.csv')
readmitTransferDict = importTransferTable('transfer_matrix_readmit_normalized.csv')

transferDict = {}
for k, rec in directTransferDict.items():
    transferDict[k] = rec.copy()
for k, r2 in readmitTransferDict.items():
    r1 = transferDict[k]
    for k2, v in r2.items():
        r1[k2] += v

junkKeys, facRecs = yaml_tools.parse_all('/home/welling/workspace/pyRHEA/models/LemonCounty/facilityfacts7')
facDict = {r['abbrev']: r for r in facRecs}

transOutDict = {}
transInDict = {}
invertDict = {}
for src, r in transferDict.items():
    totTransfersOut = float(sum([v for k, v in r.items()]))
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

inclusionSet = [('NURSINGHOME', 'HOSPITAL'),
                ('NURSINGHOME', 'LTAC'),
                ('LTAC', 'NURSINGHOME'),
                ('HOSPITAL', 'NURSINGHOME'),
                ('HOSPITAL', 'HOSPITAL'),
                ('NURSINGHOME', 'NURSINGHOME'),
                ('LTAC','LTAC'),
                ('LTAC','HOSPITAL'),
                ('HOSPITAL', 'LTAC')]
#title = "HOSPITAL + LTAC internal direct transfers"
#title = 'NURSINGHOME internal direct transfers'
title = "NURSINGHOME to/from HOSPITAL+LTAC direct transfers"
minWtCutoff = 0.1
minCntCutoff = 2
wtScale = 0.1

createdSet = set()
with Graph('graph.dot', title=title) as g:
    for src, r in transferDict.items():
        for dst, v in r.items():
            if transOutDict[src] > 0.0 and transInDict[dst] > 0.0:
                reverseV = invertDict[src][dst]
                visWt = float(v)/min(transOutDict[src], transInDict[dst])
                if (visWt >= minWtCutoff
                        and ((facDict[src]['category'], facDict[dst]['category'])
                             in inclusionSet)):
                    if src not in createdSet:
                        g.addNode(src, facToNodeAttrs(facDict[src]))
                        createdSet.add(src)
                    if dst not in createdSet:
                        g.addNode(dst, facToNodeAttrs(facDict[dst]))
                        createdSet.add(dst)
                    if (src, dst) not in createdSet:
                        if v >= reverseV:
                            symV = reverseV
                            asymV = v - reverseV
                            if symV >= minCntCutoff:
                                g.addEdge(src, dst, {'weight': wtScale * symV, 'dir': 'both',
                                                     'label': symV, 'color':'black'})
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
    print '%s (%s): %d %d' % (src, facDict[src]['category'], transInDict[src], transOutDict[src])
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
print 'Total outgoing transfers: %s (%s hospitals, %s NH)' % (totalOut, totalOutHosp, totalOutNH)
print 'Total beds: %s (%s hospitals, %s NH)' % (totBeds, totBedsHosp, totBedsNH)
