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
        width = None
        shape = None
        color = None
        style = None
        if attrDict['category'] == 'HOSPITAL':
            shape = 'box'
            try:
                width = 0.01*attrDict['meanPop']
            except:
                print '%s has no meanPop' % attrDict['abbrev']
                color = 'blue'
                style = 'dashed'
        elif attrDict['category'] == 'NURSINGHOME':
            shape = 'ellipse'
            try:
                width = 0.01*attrDict['nBeds']
            except:
                print '%s has no nBeds' % attrDict['abbrev']
                color = 'blue'
                style = 'dashed'
        else:
            shape = 'triangle'
            try:
                width = 0.01*attrDict['meanPop']
            except:
                print '%s has no meanPop' % attrDict['abbrev']
                color = 'blue'
                style = 'dashed'

        attrStr = 'label=%s' % nodeName
        for var, key in [(width, 'width'), (color, 'color'), (shape, 'shape'), (style, 'style')]:
            if var is not None:
                attrStr += ', %s=%s' % (key, var)
        self.ofile.write('%s [%s];\n' % (nodeName, attrStr))

    def addEdge(self, fromNodeName, toNodeName, attrDict={}):
        attrStr = ''
        if 'weight' in attrDict:
            attrDict['_intWt'] = int(math.floor(attrDict['weight']))
        for k1, k2 in [('weight', 'penwidth'),
                       ('_intWt', 'weight'),
                       ('label', 'label'),
                       ('color', 'color'),
                       ('dir', 'dir')]:
            if k1 in attrDict:
                attrStr += ', %s=%s' % (k2, attrDict[k1])
        attrStr = attrStr[1:]  # drop leading comma
        self.ofile.write('%s -> %s [%s];\n' % (fromNodeName, toNodeName, attrStr))

with open('transfer_matrix_direct_normalized.csv', 'r') as f:
    keys, transferRecs = csv_tools.parseCSV(f)

transferDict = {}
for r in transferRecs:
    newR = {}
    for k, v in r.items():
        if k.startswith('To_'):
            newR[k[3:]] = v
    transferDict[r['']] = newR

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
            transInDict[k] = 0.0
        transInDict[k] += float(v)
        if dst not in invertDict:
            invertDict[dst] = {}
        invertDict[dst][src] = v

inclusionSet = [('NURSINGHOME', 'HOSPITAL'),
                ('NURSINGHOME', 'LTAC'),
                ('LTAC', 'NURSINGHOME'),
                ('HOSPITAL', 'NURSINGHOME')]
#title = "HOSPITAL + LTAC internal direct transfers"
#title = 'NURSINGHOME internal direct transfers'
title = "NURSINGHOME to/from HOSPITAL+LTAC direct transfers"
minWtCutoff = 0.05
minCntCutoff = 2

createdSet = set()
with Graph('graph.dot', title=title) as g:
    for src, r in transferDict.items():
        if transOutDict[src] > 0.0:
            for dst, v in r.items():
                reverseV = invertDict[src][dst]
                if v >= reverseV:
                    symV = reverseV
                    asymV = v - reverseV
                    symWt = float(symV)/transOutDict[src]
                    asymWt = float(asymV)/transOutDict[src]
                    totWt = float(v)/transOutDict[src]
                    if (totWt >= minWtCutoff
                            and ((facDict[src]['category'], facDict[dst]['category'])
                                 in inclusionSet)):
                        if src not in createdSet:
                            g.addNode(src, facDict[src])
                            createdSet.add(src)
                        if dst not in createdSet:
                            g.addNode(dst, facDict[dst])
                            createdSet.add(dst)
                        if (src, dst) not in createdSet:
                            if symWt > 0.0 and symV >= minCntCutoff:
                                g.addEdge(src, dst, {'weight': 10.0*symWt, 'dir': 'both',
                                                     'label': symV})
                            if asymWt > 0.0 and asymV >= minCntCutoff:
                                g.addEdge(src, dst, {'weight': 10.0*asymWt, 'color': 'red',
                                                     'label': asymV})
                            createdSet.add((src, dst))
                            createdSet.add((dst, src))
                else:
                    pass  # handle it when we get to dst
