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
import scipy.optimize as op
import matplotlib.pyplot as plt

from scipy.cluster.vq import whiten, kmeans, vq
from scipy.stats import lognorm


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


def lnLik(fitParms, xVec):
    """
    model is k*lognorm( sampVec, shape=sigma, scale=exp(mu), loc=0.0 ) + (1-k)*(1.0/365.0)
    """
    k = fitParms[0]
    mu = fitParms[1]
    sigma = fitParms[2]
    # print "k= %s, sigma= %s, mu= %s" % (k, sigma, mu)

    lpVec = np.log(k * lognorm.pdf(xVec, sigma, scale=math.exp(mu), loc=0.0)
                   + (1.0 - k) * (1.0 / 365.0))
    result = np.sum(lpVec)
#     print ("sum = %s for k = %s, mu = %s, sigma = %s on %d samples" %
#            (result, k, mu, sigma, len(xVec)))
    return result


def modelFit(sampVec):
    initialGuess = [1.0, math.log(33.843041745424408), 1.0689299528774443]

    def nll(*args):
        return -lnLik(*args)

    try:
        result = op.minimize(nll, initialGuess, args=sampVec,
                             bounds=[(0.005, 1.0), (0.05, None), (0.05, None)])
        print ("success %s: nll = %s for k = %s, mu = %s, sigma = %s on %d samples" %
               (result.success, result.fun, result.x[0], result.x[1], result.x[2], len(sampVec)))
        if not result.success:
            print result
        return result.x
    except Exception, e:
        print 'Exception: %s' % e
        sys.exit(e)


indexDict = {}
valVec = []
offset = 0
losListDict = importLOSTable('/home/welling/git/rhea-dante/test/nursing_home_CI_decolonization_2014/Length_of_Stay_2007_to_2009_OC_Nursing_Homes-12-18-11_SMB_with_abbrev_RHEA.csv')
tblRecs = []
for abbrev, losList in losListDict.items():
    if len(losList) >= 10:
        indexDict[abbrev] = offset
        fitVec = modelFit(losList)
        valVec.append(fitVec)
        offset += 1
        tblRecs.append({'abbrev': abbrev})

reverseMap = {v: k for k, v in indexDict.items()}

features = whiten(valVec)
book, distortion = kmeans(features, 3)
# print book
# print distortion
code, dist = vq(features, book)
# print code

trackers = ['EDNA', 'ELIZ', 'HSOU', 'NNRC', 'ORRH', 'PALM', 'SCRT']
clrs = ['red', 'blue', 'green']
fig1, axes = plt.subplots()
scatterAx = axes
fig2, histoAxes = plt.subplots(nrows=1, ncols=3)

xIndex = 0
yIndex = 2

labels = ['k', 'mu', 'sigma']

xVals = [v[xIndex] for v in valVec]

yVals = [v[yIndex] for v in valVec]

cVals = [clrs[code[i]] for i in xrange(len(valVec))]
scatterAx.scatter(xVals, yVals, c=cVals)

for abbrev, offset in indexDict.items():
    xy = (xVals[offset], yVals[offset])
    xytext = (xVals[offset]+0.0125, yVals[offset]+0.0125)
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
            samples.append(losListDict[abbrev])
    histoAxes[i].hist(samples, bins=100, range=(0.0, 400.0), stacked=True)
    histoAxes[i].set_title('locations in ' + clrs[i])

fig1.tight_layout()
fig1.canvas.set_window_title("Clustering")
fig2.tight_layout()
fig2.canvas.set_window_title("Cluster LOS Histograms")
plt.show()
