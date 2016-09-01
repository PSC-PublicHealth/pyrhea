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

"""
Estimate LOS model distribution parameters by optimizing the Kolomogorov-Smirnov statistic
"""

import math
import sys
import os.path
import phacsl.utils.formats.csv_tools as csv_tools

import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import matplotlib.path as path
import matplotlib.patches as patches

from scipy.cluster.vq import whiten, kmeans, vq
from scipy.stats import lognorm, expon

# Truncate the sample set at how many days during fitting?
truncLim = 300


def importLOSHistoTable(fname):
    """Parse the input csv file into histogram bands by facility ID"""
    with open(fname, 'r') as f:
        keys, losRecs = csv_tools.parseCSV(f)  # @UnusedVariable
    losHistoDict = {}
    for rec in losRecs:
        key = rec['UNIQUE_ID']
        band = (rec['lowbound'], rec['highbound'], rec['cumulative_total'])
        if key in losHistoDict:
            losHistoDict[key].append(band)
        else:
            losHistoDict[key] = [band]
    # Filter to fill gaps with zeros
    for key in losHistoDict:
        newLst = []
        lst = losHistoDict[key]
        lst.sort()
        low, high, ct = lst[0]
        if low != 1:
            newLst.append((1, low-1, 0))
        assert high >= low, ('LOS table %s bin %s has invalid width' %
                             (key, low))
        newLst.append((low, high, ct))
        for band in lst[1:]:
            assert band[0] >= high + 1, ('LOS table %s bin %s is not sequential' %
                                         (key, band[0]))
            if band[0] != high + 1:
                newLst.append((high+1, band[0] - 1, 0))
            assert band[1] >= band[0], ('LOS table %s bin %s has invalid width' %
                                        (key, low))
            newLst.append(band)
            low = band[0]
            high = band[1]
        losHistoDict[key] = newLst
    return losHistoDict

def fullPDF(fitParms, xVec):
    xArr = np.asarray(xVec, np.float64)
    mu = fitParms[0]
    sigma = fitParms[1]
    # print "fullPDF: sigma= %s, mu= %s" % (sigma, mu)
    pVec = lognorm.pdf(xArr, sigma, scale=math.exp(mu), loc=0.0)
    return pVec


def fullLogPDF(fitParms, xVec):
    xArr = np.asarray(xVec, np.float64)
    mu = fitParms[0]
    sigma = fitParms[1]
    # print "fullLogPDF: sigma= %s, mu= %s" % (sigma, mu)
    lnTerm = lognorm.logpdf(xArr, sigma, scale=math.exp(mu), loc=0.0)
    return lnTerm


def fullCDF(fitParms, limVec):
    limArr = np.asarray(limVec, np.float64)
    mu = fitParms[0]
    sigma = fitParms[1]
    cdfVec = lognorm.cdf(limArr, sigma, scale=math.exp(mu), loc=0.0)
    return cdfVec


def lnLik(fitParms, xVec, dummy):
    """
    model is lognorm( sampVec, shape=sigma, scale=exp(mu), loc=0.0 )
    """
    result = np.sum(fullLogPDF(fitParms, xVec))
    return result


def truncatedLnLik(fitParms, xVec, truncLim):
    cdfAtBound = fullCDF(fitParms, truncLim)
    lpVec = np.log(fullPDF(fitParms, xVec))
    result = np.sum(lpVec) - len(xVec) * np.log(cdfAtBound)
    # print '%s -> %s %s -> %s' % \
    #   (fitParms, np.sum(lpVec),  len(xVec) * np.log(cdfAtBound), result)
    return result

gblDebug = False

def kSScore(fitParms, losHistoList, truncLim=None):
    """
    This is actually the maximum difference between the predicted number
    of counts in a histo band and the actual number of counts.  It would
    become a Kolomagorov-Smirnov score if we knew there were no gaps
    between bands and the bands were suitably narrow.
    """
    try:
        ctVec = np.asarray([float(ct) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                           np.float64)
        tot = np.sum(ctVec)
        fracVec = ctVec / tot
        lowVec = np.asarray([float(llim) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                            np.float64) + 0.0001
        highVec = np.asarray([float(hlim) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                             np.float64) + 0.9999
        lowCDFs = fullCDF(fitParms, lowVec)
        highCDFs = fullCDF(fitParms, highVec)
        predictedVec = (highCDFs - lowCDFs)
        maxSep = np.max(np.fabs(predictedVec - fracVec))
        if gblDebug:
            print '------------------'
            print fitParms
            print fracVec[:8]
            print lowVec[:8]
            print highVec[:8]
            print lowCDFs[:8]
            print highCDFs[:8]
            print predictedVec[:8]
            print maxSep
            
    except Exception, e:
        print 'Exception: %s' % e
        print fitParms
        print ctVec
        print lowVec
        print highVec
        print lowCDFs
        print highCDFs
        print predictedVec
        print maxSep
        sys.exit('failed')
    return maxSep


def modelFit(losHistoList, fun=kSScore, truncLim=None):
    initialGuess = [3.246, 0.5]
    bounds = [(0.05, 400.0), (0.05, None)]
#     initialGuess = [math.log(33.843041745424408), 1.0689299528774443]
#     bounds = [(0.05, None), (0.05, None)]

    def minimizeMe(*args):
        """Construct the function to be minimized"""
        return fun(*args)

    try:
        result = op.minimize(minimizeMe, initialGuess,
                             args=(losHistoList, truncLim),
                             bounds=bounds)
        print ("success %s: min val = %s for %s on %d bands" %
               (result.success, result.fun, result.x, len(losHistoList)))
        if not result.success:
            print result
        return result.x, result.fun, result.success
    except Exception, e:
        print 'Exception: %s' % e
        raise


def truncatedModelFit(sampVec, fun=truncatedLnLik):
    return modelFit([s for s in sampVec if s <= truncLim], fun=fun, truncLim=truncLim)


def plotCurve(axes, parmVec, scale, rng, nbins, pattern='r-'):
    curveX = np.linspace(rng[0], rng[1], nbins)
#     curveY = (fullPDF(parmVec, curveX) * scale * ((rng[1]-rng[0])/nbins))
    curveY = fullPDF(parmVec, curveX) * scale
    axes.plot(curveX, curveY, pattern, lw=2, alpha=0.6)


def resampleBand(bandTuple, sampVec):
    lo, hi, ct = bandTuple
    if hi > sampVec.shape[0]:
        sampVec = np.pad(sampVec, (0, hi + 1 - sampVec.shape[0]),
                         'constant',
                         constant_values=(0.0, 0.0))
    sval = float(ct) / (hi + 1 - lo)
    for i in xrange(lo, hi+1):
        sampVec[i] += sval
    return sampVec


def plotAsHistogram(vec, axes):
    """
    vec is assumed to be an array representing sample counts
    in unit-width bins from 0.0 to float(vec.shape[0]).  axes
    is a pyplot axes structure.
    """
    # get the corners of the rectangles for the histogram
    bins = np.asarray(xrange(vec.shape[0] + 1), np.float64)
    left = np.array(bins[:-1])
    right = np.array(bins[1:])
    bottom = np.zeros(len(left))
    top = bottom + vec

    # we need a (numrects x numsides x 2) numpy array for the path helper
    # function to build a compound path
    XY = np.array([[left, left, right, right], [bottom, top, top, bottom]]).T

    # get the Path object
    barpath = path.Path.make_compound_path_from_polys(XY)

    # make a patch out of it
    patch = patches.PathPatch(barpath, facecolor='blue', edgecolor='gray', alpha=0.8)
    axes.add_patch(patch)

    # update the view limits
    axes.set_xlim(left[0], right[-1])
    axes.set_ylim(bottom.min(), top.max())
    

def main():

    modelDir = '/home/welling/git/pyrhea/models/ChicagoLand'
    losHistoPath = os.path.join(modelDir, 'HistogramTable_ActualLOS_PROTECT_082516.csv')
    losHistoDict = importLOSHistoTable(losHistoPath)

    tblRecs = []
    indexDict = {}
    valVecList = []
    offset = 0
    for abbrev, losHistoList in losHistoDict.items():
        if len(losHistoList) >= 2:
            indexDict[abbrev] = offset
            print '%s:' % abbrev,
            fitVec, accuracyMeasure, success = modelFit(losHistoList)
            valVecList.append(fitVec)
            offset += 1
            nSamples = sum([band[2] for band in losHistoList])
            nBins = len(losHistoList)
            note = "" if success else "fitting iteration did not converge"
            tblRecs.append({'abbrev': abbrev, 'KSScore': accuracyMeasure,
                            'mu': fitVec[0], 'sigma': fitVec[1],
                            'nsamples': nSamples, 'nbins': nBins,
                            'notes': note})
            print tblRecs[-1]
        else:
            print ('No fit for %s: %s histo bins is not enough' %
                   (abbrev, len(losHistoList)))

#     ofName = 'nh_los_model_fit_parms.csv'
    ofName = 'los_model_fit_parms_ks.csv'
    print 'writing summary file %s' % ofName
    with open(ofName, 'w') as f:
        csv_tools.writeCSV(f, ['abbrev', 'mu', 'sigma', 'KSScore',
                               'nsamples', 'nbins', 'notes'],
                           tblRecs)

#     for losList in losListDict.values():
#         for v in losList:
#             aggregateLosList.append(v)
#     print 'aggregated: ',
#     aggregateFitVec, aggregateNllPerSamp = modelFit(aggregateLosList)  # @UnusedVariable

    reverseMap = {val: key for key, val in indexDict.items()}

    features = whiten(valVecList)
    book, distortion = kmeans(features, 3)  # @UnusedVariable
    print 'Book: %s' % book
    print 'Distortion: %s' % distortion
    code, dist = vq(features, book)  # @UnusedVariable
    print 'code: %s' % code

    trackers = []
    pairs = [(rec['KSScore'], rec['abbrev']) for rec in tblRecs]
    pairs = sorted(pairs, reverse=True)
    print 'Top six pairs:'
    for val, key in pairs[:6]:
        print '  %s: %s' % (key, val)

    trackers = [key for val, key in pairs[:6]]
#     trackers = ['GGCH', 'STAN', 'TOWN', 'NEWO', 'ELIZ', 'LAKE', 'TERR']

    clrs = ['red', 'blue', 'green', 'yellow']
    fig1, scatterAx = plt.subplots()
    fig2, histoAxes = plt.subplots(nrows=1, ncols=3)
    fig3, locAxes = plt.subplots(nrows=1, ncols=len(trackers))

    xIndex = 0
    yIndex = 1
    labels = ['mu', 'sigma']

    histoRange = (0.0, 400.0)

    xVals = [v[xIndex] for v in valVecList]

    yVals = [v[yIndex] for v in valVecList]

    cVals = [clrs[code[i]] for i in xrange(len(valVecList))]
    scatterAx.scatter(xVals, yVals, c=cVals)

    for abbrev, offset in indexDict.items():
        xy = (xVals[offset], yVals[offset])
        xytext = (xVals[offset]+0.0125, yVals[offset]+0.0125)
        scatterAx.annotate(abbrev, xy=xy, xytext=xytext)

    xT = []
    yT = []
    cT = []
    for i, x, y, c in zip(xrange(len(valVecList)), xVals, yVals, cVals):
        if reverseMap[i] in trackers:
            xT.append(x)
            yT.append(y)
            cT.append(c)
    scatterAx.scatter(xT, yT, c=cT, marker='+', s=200)

    unwhitenedBook = book * np.std(valVecList, axis=0)[None, :]
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
                samples.extend(losHistoDict[abbrev])
        vec = np.zeros(400)
        if samples:
            for samp in samples:
                vec = resampleBand(samp, vec)
            plotAsHistogram(vec, histoAxes[i])
            allSampList = []
            for step in xrange(vec.shape[0]):
                allSampList.append((step, step, vec[step]))
            global gblDebug
            gblDebug = (i == 0)
            print '%s:' % clrs[i],
            fitVec, accuracyMeasure, success = modelFit(allSampList)  # @UnusedVariable
            plotCurve(histoAxes[i], fitVec, np.sum(vec), rng=histoRange,
                      nbins=vec.shape[0])
        histoAxes[i].set_title('locations in ' + clrs[i])

    nbins = 50
    for i in xrange(len(trackers)):
        abbrev = trackers[i]
        vec = np.zeros(400)
        for band in losHistoDict[abbrev]:
            vec = resampleBand(band, vec)
        plotAsHistogram(vec, locAxes[i])
        locAxes[i].set_title(abbrev)
        fitParms = valVecList[indexDict[abbrev]]
        plotCurve(locAxes[i], fitParms, np.sum(vec),
                  rng=histoRange, nbins=nbins)

    fig1.tight_layout()
    fig1.canvas.set_window_title("Clustering")
    fig2.tight_layout()
    fig2.canvas.set_window_title("Cluster LOS Histograms")
    fig3.canvas.set_window_title("Tracked Locations")
    plt.show()

############
# Main hook
############

if __name__ == "__main__":
    main()
