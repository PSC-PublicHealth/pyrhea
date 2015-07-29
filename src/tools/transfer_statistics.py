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
import math

import numpy as np
import matplotlib.pyplot as plt

from scipy.cluster.vq import whiten, kmeans, vq



indexDict = {}
valVec = []
offset = 0
facKeys, facRecs= yaml_tools.parse_all('/home/welling/workspace/pyRHEA/models/OrangeCounty/facilityfacts10')
facDict = {r['abbrev']: r for r in facRecs}


tblRecs = []
for abbrev in facDict.keys():
    try:
        bSVec = [(facDict[abbrev]['meanPop'] 
                  if 'meanPop' in facDict[abbrev] 
                  else float(facDict[abbrev]['nBeds'])),
                 float(facDict[abbrev]['totalTransfersIn']['value'])/facDict[abbrev]['totalDischarges']['value'],
                 sum([float(thing['count']['value']) for thing in facDict[abbrev]['totalTransfersOut']])/facDict[abbrev]['totalDischarges']['value']
                 ]
        valVec.append(bSVec)
        indexDict[abbrev] = offset
        offset += 1
        tblRecs.append({'abbrev': abbrev
                        })
    except Exception, e:
        print 'Cannot map %s: %s' % (abbrev, e)

reverseMap = {v: k for k, v in indexDict.items()}

features = whiten(valVec)
print 'finished whiten'
book, distortion = kmeans(features, 3)
print 'finished kmeans'
# print book
# print distortion
code, dist = vq(features, book)
# print code

categories = ['HOSPITAL', 'LTAC', 'NURSINGHOME']

#trackers = ['EDNA', 'ELIZ', 'HSOU', 'NNRC', 'ORRH', 'PALM', 'SCRT']
trackers = ['TUSH', 'KINW', 'KINB', 'HSOU', 'COLL']
clrs = ['red', 'blue', 'green']
fig1, axes = plt.subplots()
scatterAx = axes
fig2, histoAxes = plt.subplots(nrows=1, ncols=3)

xIndex = 1
yIndex = 2

labels = ['approx pop', 'transfers in / discharges', 'transfers out / discharges']

xVals = [v[xIndex] for v in valVec]

yVals = [v[yIndex] for v in valVec]

#cVals = [clrs[code[i]] for i in xrange(len(valVec))]
cVals = [clrs[categories.index(facDict[reverseMap[i]]['category'])] for i in xrange(len(xVals))]
scatterAx.scatter(xVals, yVals, c=cVals)

for abbrev, offset in indexDict.items():
    xy = (xVals[offset], yVals[offset])
    #xytext = (xVals[offset]+0.5, yVals[offset]+0.05)
    xytext = (xVals[offset]+0.02, yVals[offset]+0.02)
    scatterAx.annotate(abbrev, xy=xy, xytext=xytext)

xT = []
yT = []
cT = []
for i, x, y, c in zip(xrange(len(valVec)), xVals, yVals, cVals):
    if reverseMap[i] in trackers:
        xT.append(x)
        yT.append(y)
        cT.append(c)
scatterAx.scatter(xT, yT, c=cT, marker='+', s=200)

unwhitenedBook = book * np.std(valVec, axis=0)[None, :]
xCtr = unwhitenedBook[:, xIndex]
yCtr = unwhitenedBook[:, yIndex]
clrCtr = [clrs[i] for i in xrange(xCtr.shape[0])]
scatterAx.scatter(xCtr, yCtr, marker='^', c=clrCtr, s=200)
scatterAx.grid(True)
scatterAx.set_title("LOS distribution characteristics by facility")
scatterAx.set_xlabel(labels[xIndex])
scatterAx.set_ylabel(labels[yIndex])

for i in xrange(xCtr.shape[0]):
    samples = []
    for abbrev, offset in indexDict.items():
        if code[offset] == i:
            samples.append(valVec[indexDict[abbrev]])
    histoAxes[i].hist(samples, bins=100, range=(0.0, 400.0), stacked=True)
    histoAxes[i].set_title('locations in ' + clrs[i])

fig1.tight_layout()
fig1.canvas.set_window_title("Clustering")
fig2.tight_layout()
fig2.canvas.set_window_title("Cluster LOS Histograms")
plt.show()
