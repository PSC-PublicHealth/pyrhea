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


import csv_tools
import yaml_tools
import math

import numpy as np
import matplotlib.pyplot as plt

from scipy.cluster.vq import whiten, kmeans, vq


def importLOSTable(fname):
    with open(fname, 'r') as f:
        keys, losRecs = csv_tools.parseCSV(f)  # @UnusedVariable
        
    losListDict = {}
    for r in losRecs:
        k = r['# Abbreviation']
        v = r['LOS']
        if k not in losListDict:
            losListDict[k] = []
        losListDict[k].append(v)
#     for k,v in losListDict.items():
#         print '%s: %d' % (k,len(v))
    return losListDict

def fiveNumberSummary(vec):
    return [ len(vec),
            np.mean(vec),
            np.std(vec),
            np.percentile(vec,25.0),
            np.median(vec),
            np.percentile(vec,75.0),
            ]
    
def betterSummary(vec):
    fNS = fiveNumberSummary(vec)
    N = fNS[0]
    mean = fNS[1]
    stdv = fNS[2]
    q1 = fNS[3]
    median = fNS[4]
    q3 = fNS[5]
    #return [median/mean, ((q3-q1)/stdv), q1/median, q3/median]
    #return [median, median/mean, (q3-median)/(median-q1)]
    return [median, median/mean, (q3-q1)/median]

indexDict = {}
valVec = []
offset = 0
losListDict = importLOSTable('/home/welling/git/rhea-dante/test/nursing_home_CI_decolonization_2014/Length_of_Stay_2007_to_2009_OC_Nursing_Homes-12-18-11_SMB_with_abbrev_RHEA.csv')
tblRecs = []
for abbrev, losList in losListDict.items():
    if len(losList) >= 10:
        indexDict[abbrev] = offset
        bSVec = betterSummary(losList)
        valVec.append(bSVec)
        offset += 1
        tblRecs.append({'abbrev': abbrev,
                        'median': bSVec[0],
                        'medianOverMean': bSVec[1],
                        'quartileBalance': bSVec[2]
                        })

reverseMap = {v: k for k, v in indexDict.items()}
# with open('distshapes.csv', 'w') as f:
#     csv_tools.writeCSV(f,
#                        ['abbrev', 'median', 'medianOverMean', 'quartileBalance'],
#                        tblRecs)
    
# for k,ind in indexDict.items():
#     print '%s: %s' % (k, valVec[ind])

features = whiten(valVec)
print 'finished whiten'
book, distortion = kmeans(features, 3)
print 'finished kmeans'
# print book
# print distortion
code, dist = vq(features, book)
# print code

trackers = ['EDNA', 'ELIZ', 'HSOU', 'NNRC', 'ORRH', 'PALM', 'SCRT']
clrs = ['red', 'blue', 'green']
fig, ax = plt.subplots()

xIndex = 0
yIndex = 2

labels = ['median', 'median/mean', 'q3-q1/median']

xVals = [v[xIndex] for v in valVec]
plt.xlabel(labels[xIndex])

yVals = [v[yIndex] for v in valVec]
plt.ylabel(labels[yIndex])

cVals = [clrs[code[i]] for i in xrange(len(valVec))]
ax.scatter(xVals, yVals, c=cVals)

for abbrev, offset in indexDict.items():
    xy = (xVals[offset], yVals[offset])
    xytext = (xVals[offset]+0.5, yVals[offset]+0.05)
    ax.annotate(abbrev, xy=xy, xytext=xytext)

xT = []
yT = []
cT = []
for i, x, y, c in zip(xrange(len(valVec)), xVals, yVals, cVals):
    if reverseMap[i] in trackers:
        xT.append(x)
        yT.append(y)
        cT.append(c)
ax.scatter(xT, yT, c=cT, marker='+', s=200)

# xCtr = book[:, xIndex]
# yCtr = book[:, yIndex]
# ax.scatter(xCtr, yCtr, marker='^', s=200)

ax.grid(True)
fig.tight_layout()
plt.show()


