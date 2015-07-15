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

import sys
import csv_tools
import yaml_tools
import math

import numpy as np
import matplotlib.pyplot as plt

from scipy.cluster.vq import whiten, kmeans, vq
from scipy.stats import lognorm

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


def clrCurves(xPts, yPts, clrPts, ctPts, clr):
    xxPts = []
    yyPts = []
    ccPts = []
    cctPts = []
    for x, y, c, ct in zip(xPts, yPts, clrPts, ctPts):
        if c == clr:
            xxPts.append(x)
            yyPts.append(y)
            ccPts.append(c)
            cctPts.append(ct)
    # Scale the y values
    s = sum(yyPts)
    yyPts = [y/s for y in yyPts]
    return xxPts, yyPts, ccPts, cctPts


junkKeys, facRecs = yaml_tools.parse_all('/home/welling/workspace/pyRHEA/models/OrangeCounty/facilityfacts10')
facDict = {r['abbrev']: r for r in facRecs}

directTransferDict = importTransferTable('transfer_matrix_direct_normalized.csv')

try:
    with open('/home/welling/workspace/pyRHEA/models/OrangeCounty/transitmatrix.csv', 'rU') as f:
        transitKeys, transitRecs = csv_tools.parseCSV(f)  #
except Exception, e:
    print 'Could not parse input transitmatrix'
    transitKeys = ['From', 'To', 'Seconds', 'Meters']
    transitRecs = []
transitDict = {}
for r in transitRecs:
    if r['From'] not in transitDict:
        transitDict[r['From']] = {}
    transitDict[r['From']][r['To']] = {'Seconds': r['Seconds'], 'Meters': r['Meters']}

clrMap = {'HOSPITAL:HOSPITAL': 'red',
          'HOSPITAL:LTAC': 'red',
          'LTAC:HOSPITAL': 'red',
          'LTAC:LTAC': 'red',
          'HOSPITAL:NURSINGHOME': 'blue',
          'LTAC:NURSINGHOME': 'blue',
          'NURSINGHOME:NURSINGHOME': 'yellow',
          'NURSINGHOME:HOSPITAL': 'green',
          'NURSINGHOME:LTAC': 'green'}
xPts = []
yPts = []
clrPts = []
ctPts = []

if len(sys.argv) > 1:
    facVec = [f for f in sys.argv[1:] if f in directTransferDict]
else:
    facVec = facDict.keys()
for fmLoc in facVec:
    if fmLoc in directTransferDict:
        for toLoc in directTransferDict[fmLoc]:
            if toLoc == fmLoc:
                continue
            if toLoc in facDict:
                ttoRec = facDict[fmLoc]['totalTransfersOut']
                if fmLoc in transitDict and toLoc in transitDict[fmLoc]:
                    x, y = transitDict[fmLoc][toLoc]['Seconds'], directTransferDict[fmLoc][toLoc]
                elif toLoc in transitDict and fmLoc in transitDict[toLoc]:
                    x, y = transitDict[toLoc][fmLoc]['Seconds'], directTransferDict[fmLoc][toLoc]
                else:
                    print 'No transit data for %s -> %s' % (fmLoc, toLoc)
                xPts.append(x)
                yPts.append(float(y))
                clrPts.append(clrMap['%s:%s' %
                                     (facDict[fmLoc]['category'],
                                      facDict[toLoc]['category'])])
                ctPts.append(y)
            else:
                continue
    else:
        print 'No information for %s' % fmLoc


fig1, axes = plt.subplots(nrows=2, ncols=2)

labels = ['Travel Time (seconds)', 'Fraction of Direct Transfers']

for yId, xId, clr, lbl in [(0, 0, 'red', 'Hospital to Hospital'),
                           (0, 1, 'blue', 'Hospital to Nursing Home'),
                           (1, 0, 'green', 'Nursing Home to Hospital'),
                           (1, 1, 'yellow', 'Nursing Home to Nursing Home')
                           ]:
    xxPts, yyPts, ccPts, cctPts = clrCurves(xPts, yPts, clrPts, ctPts, clr)
    axes[xId, yId].scatter(xxPts, yyPts, c=ccPts)
    axes[xId, yId].grid(True)
    axes[xId, yId].set_title("%s (%d transfers)" % (lbl, sum(cctPts)))
    axes[xId, yId].set_xlabel(labels[0])
    axes[xId, yId].set_ylabel(labels[1])

fig1.tight_layout()
if len(sys.argv) > 1:
    fig1.canvas.set_window_title(', '.join(facVec))
else:
    fig1.canvas.set_window_title("Regression")
plt.show()
