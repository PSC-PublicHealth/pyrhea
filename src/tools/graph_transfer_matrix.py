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
                attrStr += ',%s=%s' % (key, var)
        self.ofile.write('%s [%s];\n' % (nodeName, attrStr))

    def addEdge(self, fromNodeName, toNodeName, attrDict={}):
        attrStr = ''
        if 'penwidth' in attrDict:
            attrStr += ',penwidth=%f, weight=%d' % (attrDict['penwidth'],int(math.floor(attrDict['penwidth'])))
        attrStr = attrStr[1:]  # drop leading comma
        self.ofile.write('%s -> %s [%s];\n' % (fromNodeName, toNodeName, attrStr))

with open('transfer_matrix_direct_normalized.csv', 'r') as f:
    keys, transferRecs = csv_tools.parseCSV(f)
transferDict = {r['']: r for r in transferRecs}

junkKeys, facRecs = yaml_tools.parse_all('/home/welling/workspace/pyRHEA/models/LemonCounty/facilityfacts7')
facDict = {r['abbrev']: r for r in facRecs}

culledTransferDict = {}
for src, r in transferDict.items():
    culledTransferDict[src] = []
    totTransfers = float(sum([v for k, v in r.items() if k.startswith('To_')]))
    if totTransfers == 0.0:
        print '%s has no outgoing transfers' % src
        continue
    for k, v in r.items():
        if k.startswith('To_'):
            dst = k[3:]
            assert dst in facDict
            fracTransfers = float(v)/totTransfers
            culledTransferDict[src].append((dst, fracTransfers))

#inclusionSet = ['HOSPITAL', 'LTAC']
#title = "HOSPITAL + LTAC internal direct transfers"
inclusionSet = ['NURSINGHOME']
#title = 'NURSINGHOME internal direct transfers'
title = "NURSINGHOME to/from HOSPITAL+LTAC direct transfers"

with Graph('graph.dot', title=title) as g:
    for k, attrDict in facDict.items():
        #if attrDict['category'] in inclusionSet:
        if True:
            g.addNode(k, attrDict)

    for src, tplVec in culledTransferDict.items():
        for dst, wt in tplVec:
#             if (wt > 0.01 and facDict[dst]['category'] in inclusionSet
#                     and facDict[src]['category'] in inclusionSet):
            if (wt > 0.05 and ((facDict[dst]['category'] not in inclusionSet
                                and facDict[src]['category'] in inclusionSet)
                               or (facDict[dst]['category'] in inclusionSet
                                and facDict[src]['category'] not in inclusionSet))):
                g.addEdge(src, dst, {'penwidth': 10.0*wt})
