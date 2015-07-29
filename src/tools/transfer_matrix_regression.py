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
from collections import defaultdict
import types

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


def regressionModel(rec):
    try:
        # simplest model
        return rec['timeSec']
        # NH -> Hospital regression model
    #     return -(-6.759e-05 * rec['timeSec']
    #              + 3.706e-05 * rec['toTotTransfersIn']
    #              + -1.140e-03 * rec['toMeanLOS'])
        # Hospital -> NH regression model
#         return -(-2.185e-05 * rec['timeSec']
#                  + 9.792e-05 * rec['toMeanPop'] * rec['toTransientPatientFrac'])
    except Exception, e:
        print 'Error on %s -> %s: %s' % (rec['from'], rec['to'], e)
        return None
    
def regressionLabel():
     return 'Travel Time (seconds)'
#     return "NH -> Hosp regression model"
#    return 'Hosp -> NH regression model'


junkKeys, facRecs = yaml_tools.parse_all('/home/welling/workspace/pyRHEA/models/OrangeCounty/facilityfacts11')
facDict = {r['abbrev']: defaultdict(lambda: None, r) for r in facRecs}

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

hhRecs = []
nnRecs = []
hnRecs = []
nhRecs = []
recCollections = {'HOSPITAL:HOSPITAL': hhRecs,
                  'HOSPITAL:LTAC': hhRecs,
                  'LTAC:HOSPITAL': hhRecs,
                  'LTAC:LTAC': hhRecs,
                  'HOSPITAL:NURSINGHOME': hnRecs,
                  'LTAC:NURSINGHOME': hnRecs,
                  'NURSINGHOME:NURSINGHOME': nnRecs,
                  'NURSINGHOME:HOSPITAL': nhRecs,
                  'NURSINGHOME:LTAC': nhRecs}

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
yPts2 = []

totalsDict = {}

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
                if fmLoc in transitDict and toLoc in transitDict[fmLoc]:
                    timeSec = transitDict[fmLoc][toLoc]['Seconds']
                elif toLoc in transitDict and fmLoc in transitDict[toLoc]:
                    timeSec = transitDict[toLoc][fmLoc]['Seconds']
                else:
                    print 'No transit data for %s -> %s' % (fmLoc, toLoc)
                kStr = '%s:%s' % (facDict[fmLoc]['category'], facDict[toLoc]['category'])
                fmTotTransfersOut = sum([v['count']['value'] for v in facDict[fmLoc]['totalTransfersOut']])
                newRec = {'from': fmLoc, 'to': toLoc,
                          'timeSec': timeSec, 'nTransfered': directTransferDict[fmLoc][toLoc],
                          'fmMeanPop': facDict[fmLoc]['meanPop'],
                          'toMeanPop': facDict[toLoc]['meanPop'],
                          'fmNBeds': facDict[fmLoc]['nBeds'],
                          'toNBeds': facDict[toLoc]['nBeds'],
                          'fmMeanLOS': facDict[fmLoc]['meanLOS'],
                          'toMeanLOS': facDict[toLoc]['meanLOS'],
                          'fmTotTransfersIn': facDict[fmLoc]['totalTransfersIn'],
                          'toTotTransfersIn': facDict[toLoc]['totalTransfersIn'],
                          'fmTotTransfersOut': fmTotTransfersOut,
                          }
                for loc, key in [(fmLoc, 'fmTransientPatientFrac'),
                                 (toLoc, 'toTransientPatientFrac')]:
                    if facDict[loc]['category'] in ['HOSPITAL', 'LTAC']:
                        v = 1.0
                    elif facDict[loc]['losModel'] is None:
                        v = None
                    else:
                        v = facDict[loc]['losModel']['parms'][0]
                    newRec[key] = v
                for k in newRec.keys():
                    if isinstance(newRec[k], types.DictType):
                        newRec[k] = newRec[k]['value']
                recCollections[kStr].append(newRec)
                if kStr not in totalsDict:
                    totalsDict[kStr] = 0
                xPts.append(regressionModel(newRec))
                yPts.append(float(newRec['nTransfered']))
                clrPts.append(clrMap[kStr])
                ctPts.append(newRec['nTransfered'])
                totalsDict[kStr] += newRec['nTransfered']
            else:
                continue
    else:
        print 'No information for %s' % fmLoc

for k, v in totalsDict.items():
    print '%s: %s' % (k, v)

for fName, recSet in [('hosp_hosp_transfer_data.csv', hhRecs),
                      ('hosp_nh_transfer_data.csv', hnRecs),
                      ('nh_hosp_transfer_data.csv', nhRecs),
                      ('nh_nh_transfer_data.csv', nnRecs)]:
    with open(fName, 'w') as f:
        csv_tools.writeCSV(f, recSet[0].keys(), recSet)

fig1, axes = plt.subplots(nrows=2, ncols=2)

labels = [regressionLabel(), 'Fraction of Direct Transfers']

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
