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

# Truncate the sample set at how many days during fitting?
truncLim = 300


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


def fullPDF(fitParms, xVec, truncLim):
    k = fitParms[0]
    mu = fitParms[1]
    sigma = fitParms[2]
    if truncLim is None:
        alpha = 1.0/float(max(xVec))
    else:
        alpha = 1.0/float(truncLim)
    # print "k= %s, sigma= %s, mu= %s, alpha= %s" % (k, sigma, mu, alpha)

    pVec = (k * lognorm.pdf(xVec, sigma, scale=math.exp(mu), loc=0.0)
            + (1.0 - k) * alpha)
    return pVec


def lnLik(fitParms, xVec, truncLim):
    """
    model is k*lognorm( sampVec, shape=sigma, scale=exp(mu), loc=0.0 ) + (1-k)*(alpha)
    """
    k = fitParms[0]
    mu = fitParms[1]
    sigma = fitParms[2]
    if truncLim is None:
        alpha = 1.0/float(max(xVec))
    else:
        alpha = 1.0/float(truncLim)
    # print "k= %s, sigma= %s, mu= %s, alpha= %s" % (k, sigma, mu, alpha)

    lpVec = np.log(k * lognorm.pdf(xVec, sigma, scale=math.exp(mu), loc=0.0)
                   + (1.0 - k) * alpha)
    result = np.sum(lpVec)
#     print ("sum = %s for k = %s, mu = %s, sigma = %s, alpha = %s on %d samples" %
#            (result, k, mu, sigma, alpha, len(xVec)))
    return result

def truncatedLnLik(fitParms, xVec, truncLim):
    k = fitParms[0]
    mu = fitParms[1]
    sigma = fitParms[2]
    #alpha = 1.0/float(max(xVec))
    alpha = 1.0/float(truncLim)
    cdfAtBound = (k*lognorm.cdf(1.0/alpha, sigma, scale=math.exp(mu), loc=0.0)
                  + (1-k))
    pVec = (k * lognorm.pdf(xVec, sigma, scale=math.exp(mu), loc=0.0)
            + (1.0 - k) * alpha)
    lpVec = np.log(pVec)
    result = np.sum(lpVec) - len(xVec) * np.log(cdfAtBound)
    #print '%s -> %s %s -> %s' % (fitParms, np.sum(lpVec),  len(xVec) * np.log(cdfAtBound), result)
    return result


def modelFit(sampVec, fun=lnLik, truncLim=None):
    initialGuess = [0.5, math.log(33.843041745424408), 1.0689299528774443]
    bounds = [(0.001, 0.999), (0.05, None), (0.05, None)]

    def nll(*args):
        return -fun(*args)

    try:
        result = op.minimize(nll, initialGuess, args=(sampVec, truncLim), bounds=bounds)
        print ("success %s: nll/samp = %s for %s on %d samples" %
               (result.success, result.fun/len(sampVec), result.x, len(sampVec)))
        if not result.success:
            print result
        return result.x
    except Exception, e:
        print 'Exception: %s' % e
        raise


def truncatedModelFit(sampVec, fun=truncatedLnLik):
    return modelFit([s for s in sampVec if s <= truncLim], fun=fun, truncLim=truncLim)


def plotCurve(axes, parmVec, scale, rng, nbins, pattern='r-'):
    curveX = np.linspace(rng[0], rng[1], nbins)
    curveY = (fullPDF(parmVec, curveX, truncLim) * scale * ((rng[1]-rng[0])/nbins))
    axes.plot(curveX, curveY, pattern, lw=2, alpha=0.6)


indexDict = {}
valVec = []
offset = 0
losListDict = importLOSTable('/home/welling/git/rhea-dante/test/'
                             'nursing_home_CI_decolonization_2014/'
                             'Length_of_Stay_2007_to_2009_OC_Nursing_Homes-12-18-11_SMB_with_abbrev_RHEA.csv')
tblRecs = []
aggregateLosList = []
lLD = losListDict.items()[:]
lLD.sort()
for abbrev, losList in lLD:
    if len(losList) >= 10:
        indexDict[abbrev] = offset
        print '%s:' % abbrev,
        #fitVec = modelFit(losList)
        fitVec = truncatedModelFit(losList)
        valVec.append(fitVec)
        offset += 1
        tblRecs.append({'abbrev': abbrev})
    else:
        print '%s: only %d samples' % (abbrev, len(losList))

for losList in losListDict.values():
    for v in losList:
        aggregateLosList.append(v)
print 'aggregated: ',
aggregateFitVec = truncatedModelFit(aggregateLosList)

reverseMap = {v: k for k, v in indexDict.items()}

features = whiten(valVec)
book, distortion = kmeans(features, 3)
# print book
# print distortion
code, dist = vq(features, book)
# print code

#trackers = ['EDNA', 'ELIZ', 'NNRC', 'ORRH', 'PALM', 'SCRT']
#trackers = ['WLNT', 'SNMR', 'FREE', 'COVI', 'GPCC', 'CAPO']
trackers = ['EXTW', 'ROYL', 'FREE', 'COVI', 'BGYL', 'STAN', 'NSUB']
clrs = ['red', 'blue', 'green', 'yellow']
fig1, axes = plt.subplots()
scatterAx = axes
fig2, histoAxes = plt.subplots(nrows=1, ncols=3)
fig3, locAxes = plt.subplots(nrows=1, ncols=len(trackers))
fig4, allAxes = plt.subplots()

xIndex = 0
yIndex = 2

labels = ['k', 'mu', 'sigma', 'alpha']

histoRange = (0.0, 400.0)

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
    if samples:
        histoAxes[i].hist(samples, bins=100, range=histoRange, stacked=True)
        allSamples = [v for vl in samples for v in vl]  # Thank you stackoverflow
        print '%s:' % clrs[i],
        fitVec = truncatedModelFit(allSamples)
        plotCurve(histoAxes[i], fitVec, len(allSamples), rng=histoRange, nbins=100)
    histoAxes[i].set_title('locations in ' + clrs[i])

allAxes.hist(aggregateLosList, bins=100, range=histoRange)
allAxes.set_title('All locations')
plotCurve(allAxes, aggregateFitVec, len(aggregateLosList),
          rng=histoRange, nbins=100)

# guess1 = [0.75, 3.0, 0.75]
# guess2 = [0.75, 3.25, 0.75]
# plotCurve(allAxes, guess1, len(aggregateLosList),
#           rng=histoRange, nbins=100, pattern='g-')
# plotCurve(allAxes, guess2, len(aggregateLosList),
#           rng=histoRange, nbins=100, pattern='k-')
# print 'guess1 %s gives %s' % (guess1, truncatedLnLik(guess1, aggregateLosList, truncLim))
# print 'guess2 %s gives %s' % (guess2, truncatedLnLik(guess2, aggregateLosList, truncLim))

nbins = 50
for i in xrange(len(trackers)):
    abbrev = trackers[i]
    locAxes[i].hist(losListDict[abbrev], bins=nbins, range=histoRange)
    locAxes[i].set_title(abbrev)
    fitParms = valVec[indexDict[abbrev]]
    plotCurve(locAxes[i], fitParms, len(losListDict[abbrev]),
              rng=histoRange, nbins=nbins)

fig1.tight_layout()
fig1.canvas.set_window_title("Clustering")
fig2.tight_layout()
fig2.canvas.set_window_title("Cluster LOS Histograms")
fig3.canvas.set_window_title("Tracked Locations")
plt.show()
