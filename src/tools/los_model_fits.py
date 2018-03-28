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
import traceback
import os.path
from collections import defaultdict
import phacsl.utils.formats.csv_tools as csv_tools
import tools_util as tu
from pyrheautils import pathTranslate

import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import matplotlib.path as path
import matplotlib.patches as patches

from scipy.cluster.vq import whiten, kmeans, vq
from scipy.stats import lognorm, expon, weibull_min

# Truncate the sample set at how many days during fitting?
truncLim = 365


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

def importLOSLineRecs(fname, key1=None, key2=None):
    """
    Parse the input csv file into a format matching the histogram bands produced
    by importLOSHistoTable above.  key1 must be the key for the facility abbreviation
    column; key2 must be the key for the stay duration column
    """
    if key1 is None:
        raise RuntimeError('key1 (column header 1) is required')
    if key2 is None:
        raise RuntimeError('key2 (column header 2) is required')
    with open(fname, 'r') as f:
        keys, losRecs = csv_tools.parseCSV(f)  # @UnusedVariable
    assert key1 in keys, 'column header key 1 is incorrect'
    assert key2 in keys, 'column header key 2 is incorrect'
    losHistoDict = defaultdict(list)
    for rec in losRecs:
        fac = rec[key1]
        band = (rec[key2], rec[key2] + 1, 1)
        losHistoDict[fac].append(band)
    return losHistoDict


class LOSModel(object):
    """
    Base class for a collection of possible statistical models for LOS distributions
    """
    def __init__(self, initialVals=None, optimizationBounds=None):
        if not initialVals:
            initialVals = []
        if not optimizationBounds:
            optimizationBounds = [(0.05, 400.0), (0.05, None)]
        self.fitParms = np.asarray(initialVals, np.float64)
        self.optBounds = optimizationBounds[:]
        self.nameToIndexMap = {}

    def fullPDF(self, xVec, fitParms=None):
        raise RuntimeError('base fullPDF method called')
    def fullLogPDF(self, xVec, fitParms=None):
        raise RuntimeError('base fullLogPDF method called')
    def fullCDF(self, xVec, fitParms=None):
        raise RuntimeError('base fullCDF method called')
    def strDesc(self):
        raise RuntimeError('base strDesc method called')
    def setFitParms(self, fitParms):
        self.fitParms = fitParms.copy()

    def modelFit(self, losHistoList, funToMinimize=None, funToMaximize=None,
                 truncLim=None, debug=False):
        if funToMinimize:
            if funToMaximize:
                raise RuntimeError('Only one of funToMinimize and'
                                   ' funToMaximize can be specified')
            else:
                minimizeMe = funToMinimize
        else:
            if funToMaximize:
                def minimizeMe(*args):
                    return -funToMaximize(*args)
            else:
                raise RuntimeError('Must specify either funToMinimize'
                                   ' or funToMaximize')
    
        try:
            initialGuess = self.fitParms.copy()
            bounds = self.optBounds[:]
            result = op.minimize(minimizeMe, initialGuess,
#                                  method='TNC',
#                                  options={'maxiter': 1000},
                                 args=(losHistoList, self, truncLim, debug),
                                 bounds=bounds)
            nSamp = sum([ct for lo, hi, ct in losHistoList])  # @UnusedVariable
            if not result.success:
                print result
                # On failure the min value seen is not set; fix that
                result.fun = minimizeMe(result.x, losHistoList, self, truncLim)
            self.setFitParms(result.x)
            if minimizeMe == funToMinimize:
                bestVal = result.fun
            else:
                bestVal = -result.fun
            print ("success %s: min val = %s per sample for %s on %d bands, %d samples" %
                   (result.success, bestVal/nSamp, result.x,
                    len(losHistoList), nSamp))
            return result.x, bestVal, result.success
        except Exception:
            traceback.print_exc(file=sys.stdout)
            sys.exit('Exception during minimization')


class WeibullLOSModel(LOSModel):
    """
    A simple exponential distribution
    """
    def __init__(self, initialVals=None, optimizationBounds=None):
        if not initialVals:
            initialVals = [0.01, 1.0]
        if not optimizationBounds:
            optimizationBounds = [(0.0001, None), (0.99, 1.01)]
        super(WeibullLOSModel, self).__init__(initialVals=initialVals,
                                              optimizationBounds=optimizationBounds)
        self.nameToIndexMap.update({'c': 0, 's': 1})

    def fullPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        c, s = fitParms[:2]
        pVec = weibull_min.pdf(xArr, c, scale=s, loc=0.0)
        return pVec

    def fullLogPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        c, s = fitParms[:2]
        lnTerm = weibull_min.logpdf(xArr, c, scale=s, loc=0.0)
        return lnTerm

    def fullCDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        c, s = fitParms[:2]
        cdfVec = weibull_min.cdf(xArr, c, scale=s, loc=0.0)
        return cdfVec

    def intervalCDF(self, lowVec, highVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        c, s = fitParms[:2]
        xMean = weibull_min.mean(c, scale=s, loc=0.0)
        lowHighFlags = lowVec >= xMean
        highLowFlags = highVec < xMean
        bothLowFlags = highLowFlags  # since low must be lower
        bothHighFlags = lowHighFlags  # ditto
        lowCDFVec = weibull_min.cdf(lowVec, c, scale=s, loc=0.0)
        highCDFVec = weibull_min.cdf(highVec, c, scale=s, loc=0.0)
        bothLowCDFVec = highCDFVec - lowCDFVec
        lowSFVec = weibull_min.sf(lowVec, c, scale=s, loc=0.0)
        highSFVec = weibull_min.sf(highVec, c, scale=s, loc=0.0)
        bothHighCDFVec = lowSFVec - highSFVec
        mixedCDFVec = 1.0 - (highSFVec + lowCDFVec)
        dCDFVec = np.select([bothLowFlags, bothHighFlags],
                            [bothLowCDFVec, bothHighCDFVec],
                            default=mixedCDFVec)
        return dCDFVec

    def strDesc(self):
        return "weibell_min(lmda=$0)"


class ExponLOSModel(LOSModel):
    """
    A simple exponential distribution
    """
    def __init__(self, initialVals=None, optimizationBounds=None):
        if not initialVals:
            initialVals = [0.01]
        if not optimizationBounds:
            optimizationBounds = [(0.0001, None)]
        super(ExponLOSModel, self).__init__(initialVals=initialVals,
                                            optimizationBounds=optimizationBounds)
        self.nameToIndexMap.update({'lmda': 0})

    def fullPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        lmda = fitParms[0]
        # print "fullPDF: lmda= %s" % lmda
        pVec = expon.pdf(xArr, scale=1.0/lmda, loc=0.0)
        return pVec

    def fullLogPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        lmda = fitParms[0]
        # print "fullLogPDF: lmda= %s" % lmda
        lnTerm = expon.logpdf(xArr, scale=1.0/lmda, loc=0.0)
        return lnTerm

    def fullCDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        lmda = fitParms[0]
        # print "fullCDF: lmda= %s" % lmda
        cdfVec = expon.cdf(xArr, scale=1.0/lmda, loc=0.0)
        return cdfVec

    def intervalCDF(self, lowVec, highVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        lmda = fitParms[0]
        invL = 1.0/lmda
        xMean = expon.mean(scale=invL, loc=0.0)
        lowHighFlags = lowVec >= xMean
        highLowFlags = highVec < xMean
        bothLowFlags = highLowFlags  # since low must be lower
        bothHighFlags = lowHighFlags  # ditto
        lowCDFVec = expon.cdf(lowVec, scale=invL, loc=0.0)
        highCDFVec = expon.cdf(highVec, scale=invL, loc=0.0)
        bothLowCDFVec = highCDFVec - lowCDFVec
        lowSFVec = expon.sf(lowVec, scale=invL, loc=0.0)
        highSFVec = expon.sf(highVec, scale=invL, loc=0.0)
        bothHighCDFVec = lowSFVec - highSFVec
        mixedCDFVec = 1.0 - (highSFVec + lowCDFVec)
        dCDFVec = np.select([bothLowFlags, bothHighFlags],
                            [bothLowCDFVec, bothHighCDFVec],
                            default=mixedCDFVec)
        return dCDFVec

    def strDesc(self):
        return "expon(lmda=$0)"


class LogNormLOSModel(LOSModel):
    """
    A simple log-normal distribution
    """
    def __init__(self, initialVals=None, optimizationBounds=None):
        if not initialVals:
            initialVals = [3.246, 0.5]
        if not optimizationBounds:
            optimizationBounds = [(0.05, 400.0), (0.05, None)]
        super(LogNormLOSModel, self).__init__(initialVals=initialVals,
                                              optimizationBounds=optimizationBounds)
        self.nameToIndexMap.update({'mu': 0, 'sigma': 1})

    def fullPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        mu = fitParms[0]
        sigma = fitParms[1]
        # print "fullPDF: sigma= %s, mu= %s" % (sigma, mu)
        pVec = lognorm.pdf(xArr, sigma, scale=math.exp(mu), loc=0.0)
        return pVec

    def fullLogPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        mu = fitParms[0]
        sigma = fitParms[1]
        # print "fullLogPDF: sigma= %s, mu= %s" % (sigma, mu)
        lnTerm = lognorm.logpdf(xArr, sigma, scale=math.exp(mu), loc=0.0)
        return lnTerm

    def fullCDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        mu = fitParms[0]
        sigma = fitParms[1]
        cdfVec = lognorm.cdf(xArr, sigma, scale=math.exp(mu), loc=0.0)
        return cdfVec

    def intervalCDF(self, lowVec, highVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        mu = fitParms[0]
        sigma = fitParms[1]
        xMean = lognorm.mean(sigma, scale=math.exp(mu), loc=0.0)
        lowHighFlags = lowVec >= xMean
        highLowFlags = highVec < xMean
        bothLowFlags = highLowFlags  # since low must be lower
        bothHighFlags = lowHighFlags  # ditto
        lowCDFVec = lognorm.cdf(lowVec, sigma, scale=math.exp(mu), loc=0.0)
        highCDFVec = lognorm.cdf(highVec, sigma, scale=math.exp(mu), loc=0.0)
        bothLowCDFVec = highCDFVec - lowCDFVec
        lowSFVec = lognorm.sf(lowVec, sigma, scale=math.exp(mu), loc=0.0)
        highSFVec = lognorm.sf(highVec, sigma, scale=math.exp(mu), loc=0.0)
        bothHighCDFVec = lowSFVec - highSFVec
        mixedCDFVec = 1.0 - (highSFVec + lowCDFVec)
        dCDFVec = np.select([bothLowFlags, bothHighFlags],
                            [bothLowCDFVec, bothHighCDFVec],
                            default=mixedCDFVec)
        return dCDFVec

    def strDesc(self):
        return "lognorm(mu=$0,sigma=$1)"


class TwoPopLOSModel(LOSModel):
    """
    The population is divided between a group following a log normal LOS
    distribution and a group with a slow exponential LOS distribution
    """
    def __init__(self, initialVals=None, optimizationBounds=None):
        if not initialVals:
            initialVals = [0.8, 3.246, 0.5, 0.0015]
        if not optimizationBounds:
            optimizationBounds = [(0.001, 0.999), (0.05, 400.0), (0.05, None), (0.00001, 0.1)]
        self.innerLNModel = LogNormLOSModel(initialVals=initialVals[1:3],
                                            optimizationBounds=optimizationBounds[1:3])
        super(TwoPopLOSModel, self).__init__(initialVals=initialVals,
                                              optimizationBounds=optimizationBounds)
        self.nameToIndexMap.update({'k': 0, 'mu': 1, 'sigma': 2, 'lmda': 3})        

    def fullPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        k = fitParms[0]
        lmda = fitParms[3]
        pVec = ((k * self.innerLNModel.fullPDF(xVec, fitParms[1:3]))
                + ((1.0 - k) * expon.pdf(xArr, scale=1.0/lmda)))
        return pVec

    def fullLogPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        k = fitParms[0]
        lmda = fitParms[3]
        # print "fullLogPDF: k= %s, sigma= %s, mu= %s, lmda= %s" % (k, sigma, mu, lmda)
        lnTerm = k * self.innerLNModel.fullPDF(xArr, fitParms[1:3])
        expTerm = (1.0 - k) * expon.pdf(xArr, scale=1.0/lmda)
        oldErrInfo = np.seterr(all='ignore')
        lnVersion = (math.log(k) + self.innerLNModel.fullLogPDF(xArr, fitParms[1:3])
                     + np.log(1.0 + (expTerm / np.where(lnTerm == 0.0, 1.0, lnTerm))))
        expVersion = (math.log(1.0-k) + expon.logpdf(xArr, scale=1.0/lmda)
                      + np.log(1.0 + (lnTerm / np.where(expTerm == 0.0, 1.0, expTerm))))
        np.seterr(**oldErrInfo)
        lpVec = np.where(lnTerm > expTerm, lnVersion, expVersion)
        return lpVec

    def fullCDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        k = fitParms[0]
        lmda = fitParms[3]
        lnTerm = k * self.innerLNModel.fullCDF(xArr, fitParms[1:3])
        expTerm = (1.0 - k) * expon.cdf(xArr, scale=1.0/lmda)
        return lnTerm + expTerm

    def intervalCDF(self, lowVec, highVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        k = fitParms[0]
        lmda = fitParms[3]
        lnTerm = self.innerLNModel.intervalCDF(lowVec, highVec, fitParms[1:3])
        
        xMean = 1.0 / lmda
        lowHighFlags = lowVec >= xMean
        highLowFlags = highVec < xMean
        bothLowFlags = highLowFlags  # since low must be lower
        bothHighFlags = lowHighFlags  # ditto
        
        lowCDFVec = expon.cdf(lowVec, scale=1.0/lmda)
        highCDFVec = expon.cdf(highVec, scale=1.0/lmda)
        bothLowCDFVec = highCDFVec - lowCDFVec

        lowSFVec = expon.sf(lowVec, scale=1.0/lmda)
        highSFVec = expon.sf(highVec, scale=1.0/lmda)
        bothHighCDFVec = lowSFVec - highSFVec
        mixedCDFVec = 1.0 - (highSFVec + lowCDFVec)
        dCDFVec = np.select([bothLowFlags, bothHighFlags],
                            [bothLowCDFVec, bothHighCDFVec],
                            default=mixedCDFVec)
        return (k * lnTerm) + ((1.0 - k) * dCDFVec)

    def modelFit(self, losHistoList, funToMinimize=None, funToMaximize=None,
                 truncLim=None, debug=False):
        resultX, bestVal, resultSuccess = \
            super(TwoPopLOSModel, self).modelFit(losHistoList,
                                                 funToMinimize=funToMinimize,
                                                 funToMaximize=funToMaximize,
                                                 truncLim=truncLim,
                                                 debug=debug)
        self.innerLNModel.setFitParms(resultX[1:3])
        return resultX, bestVal, resultSuccess

    def strDesc(self):
        return "$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)"


class TwoLogNormLOSModel(LOSModel):
    """
    The population is divided between two groups, each with its own log normal distribution.
    """
    def __init__(self, initialVals=None, optimizationBounds=None):
        if not initialVals:
            initialVals = [0.5, 2.2, 0.5, 3.3, 0.5]
        if not optimizationBounds:
            optimizationBounds = [(0.001, 0.999), (0.05, 400.0), (0.05, None),
                                  (0.05, 400.0), (0.05, None)]
        self.leftInnerLN = LogNormLOSModel(initialVals=initialVals[1:3],
                                           optimizationBounds=optimizationBounds[1:3])
        self.rightInnerLN = LogNormLOSModel(initialVals=initialVals[3:],
                                            optimizationBounds=optimizationBounds[3:])
        super(TwoLogNormLOSModel, self).__init__(initialVals=initialVals,
                                                 optimizationBounds=optimizationBounds)
        self.nameToIndexMap.update({'k': 0, 'mu': 1, 'sigma': 2, 'mu2':3, 'sigma2': 4})

    def fullPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        k = fitParms[0]
        pVec = ((k * self.leftInnerLN.fullPDF(xArr, fitParms[1:3]))
                + ((1.0 - k) * self.rightInnerLN.fullPDF(xArr, fitParms[3:])))
        return pVec

    def fullLogPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        k = fitParms[0]
        leftTerm = k * self.leftInnerLN.fullPDF(xArr, fitParms[1:3])
        rightTerm = (1.0 - k) * self.rightInnerLN.fullPDF(xArr, fitParms[3:])
        oldErrInfo = np.seterr(all='ignore')

        leftVersion = (math.log(k) + self.leftInnerLN.fullLogPDF(xArr, fitParms[1:3])
                       + np.log(1.0 + (rightTerm / np.where(leftTerm == 0.0, 1.0, leftTerm))))
        rightVersion = (math.log(1.0 - k) + self.rightInnerLN.fullLogPDF(xArr, fitParms[3:])
                        + np.log(1.0 + (leftTerm / np.where(rightTerm == 0.0, 1.0, rightTerm))))
        np.seterr(**oldErrInfo)
        lpVec = np.where(leftTerm > rightTerm, leftVersion, rightVersion)
        return lpVec

    def fullCDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        k = fitParms[0]
        leftTerm = k * self.leftInnerLN.fullCDF(xArr, fitParms[1:3])
        rightTerm = (1.0 - k) * self.rightInnerLN.fullCDF(xArr, fitParms[3:])
        return leftTerm + rightTerm

    def intervalCDF(self, lowVec, highVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        k = fitParms[0]
        leftTerm = k * self.leftInnerLN.intervalCDF(lowVec, highVec, fitParms[1:3])
        rightTerm = (1.0 - k) * self.rightInnerLN.intervalCDF(lowVec, highVec, fitParms[3:])
        
        return leftTerm + rightTerm

    def modelFit(self, losHistoList, funToMinimize=None, funToMaximize=None,
                 truncLim=None, debug=False):
        resultX, bestVal, resultSuccess = \
            super(TwoLogNormLOSModel, self).modelFit(losHistoList,
                                                     funToMinimize=funToMinimize,
                                                     funToMaximize=funToMaximize,
                                                     truncLim=truncLim,
                                                     debug=debug)
        self.leftInnerLN.setFitParms(resultX[1:3])
        self.rightInnerLN.setFitParms(resultX[3:])
        return resultX, bestVal, resultSuccess

    def strDesc(self):
        return "$0*lognorm(mu=$1,sigma=$2)+(1-$0)*lognorm(mu=$3, sigma=$4)"

class TwoExponLOSModel(LOSModel):
    """
    The population is divided between two groups, each with its own exponential distribution.
    """
    def __init__(self, initialVals=None, optimizationBounds=None):
        if not initialVals:
            initialVals = [0.5, 0.1, 0.001]
        if not optimizationBounds:
            optimizationBounds = [(0.0, 1.0), (0.01, None), (0.0001, None)]
        self.leftInnerEx = ExponLOSModel(initialVals=initialVals[1:2],
                                           optimizationBounds=optimizationBounds[1:2])
        self.rightInnerEx = ExponLOSModel(initialVals=initialVals[2:],
                                            optimizationBounds=optimizationBounds[2:])
        super(TwoExponLOSModel, self).__init__(initialVals=initialVals,
                                               optimizationBounds=optimizationBounds)
        self.nameToIndexMap.update({'k': 0, 'lmda': 1, 'lmda2': 2})

    def fullPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        k = fitParms[0]
        pVec = ((k * self.leftInnerEx.fullPDF(xArr, fitParms[1:2]))
                + ((1.0 - k) * self.rightInnerEx.fullPDF(xArr, fitParms[2:])))
        return pVec

    def fullLogPDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        k = fitParms[0]
        leftTerm = k * self.leftInnerEx.fullPDF(xArr, fitParms[1:2])
        rightTerm = (1.0 - k) * self.rightInnerEx.fullPDF(xArr, fitParms[2:])
        oldErrInfo = np.seterr(all='ignore')
        leftVersion = (math.log(k) + self.leftInnerEx.fullLogPDF(xArr, fitParms[1:2])
                       + np.log(1.0 + (rightTerm / np.where(leftTerm == 0.0, 1.0, leftTerm))))
        rightVersion = (math.log(1.0 - k) + self.rightInnerEx.fullLogPDF(xArr, fitParms[2:])
                        + np.log(1.0 + (leftTerm / np.where(rightTerm == 0.0, 1.0, rightTerm))))
        np.seterr(**oldErrInfo)
        lpVec = np.where(leftTerm > rightTerm, leftVersion, rightVersion)
        return lpVec

    def fullCDF(self, xVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        xArr = np.asarray(xVec, np.float64)
        k = fitParms[0]
        leftTerm = k * self.leftInnerEx.fullCDF(xArr, fitParms[1:2])
        rightTerm = (1.0 - k) * self.rightInnerEx.fullCDF(xArr, fitParms[2:])
        return leftTerm + rightTerm

    def intervalCDF(self, lowVec, highVec, fitParms=None):
        if fitParms is None:
            fitParms = self.fitParms
        k = fitParms[0]
        leftTerm = k * self.leftInnerEx.intervalCDF(lowVec, highVec, fitParms[1:2])
        rightTerm = (1.0 - k) * self.rightInnerEx.intervalCDF(lowVec, highVec, fitParms[2:])
        return leftTerm + rightTerm

    def modelFit(self, losHistoList, funToMinimize=None, funToMaximize=None,
                 truncLim=None, debug=False):
        resultX, bestVal, resultSuccess = \
            super(TwoExponLOSModel, self).modelFit(losHistoList,
                                                   funToMinimize=funToMinimize,
                                                   funToMaximize=funToMaximize,
                                                   truncLim=truncLim,
                                                   debug=debug)
        self.leftInnerEx.setFitParms(resultX[1:2])
        self.rightInnerEx.setFitParms(resultX[2:])
        return resultX, bestVal, resultSuccess

    def strDesc(self):
        return "$0*expon(lmda=$1)+(1-$0)*expon(lmda2=$2)"


def lnLik(fitParms, losHistoList, losModel, truncLim = None, debug=False):
    ctVec = lowVec = highVec = None
    if truncLim:
        losHistoList = [(llim, hlim, ct) for llim, hlim, ct in losHistoList
                        if hlim <= truncLim]
    try:
        ctVec = np.asarray([float(ct) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                           np.float64)
        lowVec = (np.asarray([float(llim) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                            np.float64)
                  + 0.0001)
        highVec = (np.asarray([float(hlim) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                             np.float64)
                   + 0.9999)
        lnLikVal = np.sum(ctVec * np.log(losModel.intervalCDF(lowVec, highVec, fitParms)))
        if truncLim is not None:
            cdfAtBound = losModel.fullCDF(truncLim, fitParms)
            lnLikVal -= np.sum(ctVec) * np.log(cdfAtBound)
        if debug:
            print '%s -> %s' % (fitParms, lnLikVal/np.sum(ctVec))
#             print ctVec
#             print lowVec
#             print highVec
#             print losModel.intervalCDF(lowVec, highVec, fitParms)
#             print ctVec * np.log(losModel.intervalCDF(lowVec, highVec, fitParms))
#             print 'log likelihood %s' % lnLikVal
        return lnLikVal
    except Exception:
        traceback.print_exc(file=sys.stdout)
        print fitParms
        print ctVec
        print lowVec
        print highVec
        print losModel.intervalCDF(lowVec, highVec, fitParms)
        sys.exit('lnLik failed')

def boundedLnLik(fitParms, losHistoList, losModel, truncLim = None, debug=False):
    ctVec = lowVec = highVec = None
    if truncLim:
        losHistoList = [(llim, hlim, ct) for llim, hlim, ct in losHistoList
                        if hlim <= truncLim]
    try:
        ctVec = np.asarray([float(ct) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                           np.float64)
        lowVec = (np.asarray([float(llim) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                            np.float64)
                  + 0.0001)
        highVec = (np.asarray([float(hlim) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                             np.float64)
                   + 0.9999)
        lowBnd = min([llim for llim, hlim, ct in losHistoList])
        highBnd = max([hlim for llim, hlim, ct in losHistoList])
        cdfScale = losModel.intervalCDF([lowBnd], [highBnd], fitParms)
        lnLikVal = (np.sum(ctVec * np.log(losModel.intervalCDF(lowVec, highVec, fitParms)))
                    - np.sum(ctVec) * np.log(cdfScale))
        if truncLim is not None:
            cdfAtBound = losModel.fullCDF(truncLim, fitParms)
            lnLikVal -= np.sum(ctVec) * np.log(cdfAtBound)
        if debug:
            print '%s -> %s' % (fitParms, lnLikVal/np.sum(ctVec))
#             print ctVec
#             print lowVec
#             print highVec
#             print losModel.intervalCDF(lowVec, highVec, fitParms)
#             print ctVec * np.log(losModel.intervalCDF(lowVec, highVec, fitParms))
#             print 'log likelihood %s' % lnLikVal
        return lnLikVal
    except Exception:
        traceback.print_exc(file=sys.stdout)
        print fitParms
        print ctVec
        print lowVec
        print highVec
        print losModel.intervalCDF(lowVec, highVec, fitParms)
        sys.exit('lnLik failed')


def kSScore(fitParms, losHistoList, losModel, truncLim=None):
    """
    This is actually the maximum difference between the predicted number
    of counts in a histo band and the actual number of counts.  It would
    become a Kolomagorov-Smirnov score if we knew there were no gaps
    between bands and the bands were suitably narrow.
    """
    if truncLim is not None:
        raise RuntimeError('kSScore: truncLim is not implemented')
    ctVec = fracVec = lowVec = highVec = lowCDFs = highCDFs = predictedVec = maxSep = None
    try:
        ctVec = np.asarray([float(ct) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                           np.float64)
        tot = np.sum(ctVec)
        fracVec = ctVec / tot
        lowVec = np.asarray([float(llim) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                            np.float64) + 0.0001
        highVec = np.asarray([float(hlim) for llim, hlim, ct in losHistoList],  # @UnusedVariable
                             np.float64) + 0.9999
        lowCDFs = losModel.fullCDF(lowVec, fitParms)
        highCDFs = losModel.fullCDF(highVec, fitParms)
        predictedVec = (highCDFs - lowCDFs)
        maxSep = np.max(np.fabs(predictedVec - fracVec))
    except Exception, e:
        traceback.print_exc(file=sys.stdout)
        print fitParms
        print ctVec
        print lowVec
        print highVec
        print lowCDFs
        print highCDFs
        print predictedVec
        print maxSep
        sys.exit('kSScore failed')
    return maxSep


def plotCurve(axes, parmVec, losModel, scale, rng, nbins, pattern='r-'):
    curveX = np.linspace(rng[0], rng[1], nbins)
#     curveY = (fullPDF(parmVec, curveX) * scale * ((rng[1]-rng[0])/nbins))
    curveY = losModel.fullPDF(curveX, parmVec) * scale
    axes.plot(curveX, curveY, pattern, lw=2, alpha=0.6)


def resampleBand(bandTuple, sampVec, range=None):
    """
    If given, range should be a tuple (low, high).  If range is none,
    (0, sampVec.shape[0]) is used.
    """
    if range is None:
        rLow = 0.0
        rHigh = float(sampVec.shape[0])
    else:
        rLow, rHigh = range
    sampW = sampVec.shape[0]
    lo, hi, ct = bandTuple
    if hi > rHigh or lo < rLow:
        print 'sample band (%s, %s) out of range (%s, %s)' % (lo, hi, rLow, rHigh)
    else:
        cLow = int(sampW*((lo - rLow)/(rHigh-rLow)) + 0.0001)
        cHigh = int(sampW*((hi - rLow)/(rHigh-rLow)) - 0.0001)
        sval = float(ct) / (cHigh + 1 - cLow)
        #print lo, hi, ct, sampW, cLow, cHigh, sval
        for i in xrange(cLow, cHigh+1):
            sampVec[i] += sval
    return sampVec


def plotAsHistogram(vec, axes, range=None):
    """
    vec is assumed to be an array representing sample counts.
    range gives the width of the vector; if range is none the
    vector is assumed to span from 0.0 to float(vec.shape[0]).
    axes is a pyplot axes structure.
    """
    # get the corners of the rectangles for the histogram
    if range is None:
        rLow = 0.0
        rHigh = float(vec.shape[0])
    else:
        rLow, rHigh = range
    delta = (rHigh - rLow)/vec.shape[0]
    bins = np.asarray([rLow + delta * i for i in xrange(vec.shape[0] + 1)],
                      np.float64)
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


def performClustering(valVecList, nClusters=3):
    """
    Given a list of fit parameter values, separate into nClusters groups
    by k-means clustering.  Returns a vector of cluster codes and a vector
    of cluster center information
    """
    features = whiten(valVecList)
    book, distortion = kmeans(features, nClusters)  # @UnusedVariable
#     print 'Book: %s' % book
#     print 'Distortion: %s' % distortion
    code, dist = vq(features, book)  # @UnusedVariable
#     print 'code: %s' % code
    unwhitenedBook = book * np.std(valVecList, axis=0)[None, :]
    return code, unwhitenedBook


def codeByCategory(reverseMap, facDict):
    """
    For those facilities listed in reverseMap, make up an integer code
    representing its facility type.  Returns a vector of those codes and
    a dictionary of category names by code.
    """
    codeDict = {}
    codes = []
    curCode = 0
    for i in xrange(len(reverseMap)):
        abbrev = reverseMap[i]
        category = facDict[abbrev]['category']
        if category not in codeDict:
            codeDict[category] = curCode
            curCode += 1
        codes.append(codeDict[category])
    return codes, codeDict


def makeScatterPlot(ax, markerTupleList):
    for tpl in markerTupleList:
        if len(tpl) == 5:
            xVals, yVals, cVals, mrk, sz = tpl
            ax.scatter(xVals, yVals, marker=mrk, c=cVals, s=sz)
        elif len(tpl) == 6:
            xVals, yVals, cVals, mrk, sz, label = tpl
            ax.scatter(xVals, yVals, marker=mrk, c=cVals, s=sz, label=label)
        else:
            raise RuntimeError('Cannot parse a tuple of length %d' % len(tpl))

def main():

    inputDict = tu.readModelInputs('../sim/twoyear_run_OC.yaml')
    facDict = tu.getFacDict(inputDict)
    
    #losHistoPath = os.path.join(pathTranslate('$(MODELDIR)', 'Histogram_09292016.csv')
    losHistoPath = os.path.join(pathTranslate('$(MODELDIR)'),
                                "OC_Nursing_Home_LOS_Line-Lists_for_RHEA_2.0_-_2011-2015"
                                "_-_Adult_Only_-_09-29-2017_FINAL_NH_LOS_Line_List.csv"
                                )
    losHistoDict = importLOSLineRecs(losHistoPath,
                                     key1='CODE',
                                     key2='LOS')


#    LOSModelType = LogNormLOSModel
#    LOSModelType = TwoPopLOSModel
#    LOSModelType = TwoLogNormLOSModel
#    LOSModelType = ExponLOSModel
#    LOSModelType = TwoExponLOSModel
    LOSModelType = WeibullLOSModel
    tblRecs = []
    indexDict = {}
    valVecList = []
    offset = 0
    for abbrev, losHistoList in losHistoDict.items():
        if abbrev not in facDict:
            continue  # some facilities have been dropped from consideration
        nSamples = sum([band[2] for band in losHistoList])
        nBins = len(losHistoList)
        if len([tpl for tpl in losHistoList if tpl[2] != 0]) >= 2:
            indexDict[abbrev] = offset
            print '%s:' % abbrev,
            losModel = LOSModelType()
            fitVec, accuracyMeasure, success = losModel.modelFit(losHistoList,
                                                                 funToMaximize=boundedLnLik)
            valVecList.append(fitVec)
            offset += 1
            note = "" if success else "fitting iteration did not converge"
            newRec = {'abbrev': abbrev,
                            'lnLikPerSample': accuracyMeasure/nSamples,
                            'nsamples': nSamples, 'nbins': nBins,
                            'notes': note, 'pdf': losModel.strDesc()}
            newRec.update({k: fitVec[v] for k, v in losModel.nameToIndexMap.items()})
            newRec.update({'_fitVec': fitVec})
            tblRecs.append(newRec)
            print tblRecs[-1]
        else:
            print ('No fit for %s: %s histo bins is not enough' %
                   (abbrev, len(losHistoList)))
            note = "Not enough histogram bins to estimate LOS"
            tblRecs.append({'abbrev': abbrev,
                            'lnLikPerSample': None,
                            'nsamples': nSamples, 'nbins': nBins,
                            'notes': note, 'pdf': 'NA'})

    ofName = 'los_model_fit_parms.csv'
    print 'writing summary file %s' % ofName
    with open(ofName, 'w') as f:
        csv_tools.writeCSV(f, ['abbrev', 'k', 'mu', 'sigma', 'mu2', 'sigma2',
                               'lmda', 'lmda2', 'lnLikPerSample', 'nsamples', 'nbins', 'notes'],
                           tblRecs)

    reverseMap = {val: key for key, val in indexDict.items()}

    nTrackers = 6
    nClusters = 3
    xLabel = 'c'
#    yLabel = 'lmda2'
    yLabel = 's'
    annotateScatterPts = True
    prop_cycle = plt.rcParams['axes.prop_cycle']
    clrs = prop_cycle.by_key()['color']
#     histoRange = (0.0, 400.0)
#     histoBins = 400
    histoRange = (0.0, 1000.0)
    histoBins = 1000

    tupleL = [(-rec['lnLikPerSample'], rec['abbrev'], rec['_fitVec']) for rec in tblRecs
             if rec['lnLikPerSample'] is not None]
    tupleL = sorted(tupleL, reverse=True)
    print 'Top pairs:'
    for val, key, pdfS in tupleL[:nTrackers]:
        print '  %s: %s fit by %s' % (key, -val, pdfS)
    tupleL.reverse()
    print 'Bottom pairs:'
    for val, key, pdfS in tupleL[:nTrackers]:
        print '  %s: %s fit by %s' % (key, -val, pdfS)


    trackers = ([key for val, key, pdfS in tupleL[:nTrackers/2]]  # @UnusedVariable
                + [key for val, key, pdfS in tupleL[-(nTrackers - nTrackers/2):]])  # @UnusedVariable
#     trackers = ['EMBA_555_S', 'SH_700_S', 'CLAR_700_S', 'ADVO_17800_H', 'MORR_150_H']
#     trackers = ['HEAR_3100_S', 'MISE_6300_S', 'PARK_6125_V', 'NORT_800_H', 'SOME_5009_S',
#                 'RIVE_350_H']
#    trackers = [
#                 'ADVO_3435_L', 'PRES_100_L', 'RML_5601_L', 'THC_225_L', 'THC_2544_L',
#                 'THC_365_L', 'THC_4058_L', 'THC_6130_L', 'VIBR_9509_L',
#                'PETE_520_S', 'ADVO_17800_H', 'MORR_150_H', 'THC_2544_L', 'THC_225_L']

    for lbl in [xLabel, yLabel]:
        assert lbl in losModel.nameToIndexMap, '%s is not a parameter name' % lbl
    xIndex = losModel.nameToIndexMap[xLabel]
    yIndex = losModel.nameToIndexMap[yLabel]

    xVals = [val[xIndex] for val in valVecList]

    yVals = [val[yIndex] for val in valVecList]

    code, codeDict = codeByCategory(reverseMap, facDict)
#     print codeDict
    assert len(clrs) >= len(codeDict), 'Define more colors!'
    cVals = [clrs[code[i]] for i in xrange(len(valVecList))]
    trackedVals = [(x, y, clr)
                   for i, x, y, clr in zip(xrange(len(valVecList)), xVals, yVals, cVals)
                   if reverseMap[i] in trackers]
    trXVals, trYVals, trCVals = [list(x) for x in zip(*trackedVals)]
    fig1, scatterAx = plt.subplots()
    scatterAx.grid(True)
    scatterAx.set_title("LOS distribution characteristics by facility\n"
                        "colored by category")
    scatterAx.set_xlabel(xLabel)
    scatterAx.set_ylabel(yLabel)
    scatterAx.grid(True)
    scatterSets = []
    for category, codeInt in codeDict.items():
        scatterSets.append(([x for i, x in zip(code, xVals) if i == codeInt],
                            [y for i, y in zip(code, yVals) if i == codeInt],
                            [clrs[codeInt] for i in code if i == codeInt],
                            'o', 20, category))
    scatterSets.append((trXVals, trYVals, trCVals, '+', 200))
    makeScatterPlot(scatterAx, scatterSets)
    scatterAx.legend(loc='upper right', shadow=True)
    fig1.tight_layout()
    fig1.canvas.set_window_title("By Categories")

    clCode, unwhitenedBook = performClustering(valVecList, nClusters=nClusters)
    cVals = [clrs[clCode[i]] for i in xrange(len(valVecList))]
    xCtr = unwhitenedBook[:, xIndex]
    yCtr = unwhitenedBook[:, yIndex]
    cCtr = [clrs[i] for i in xrange(xCtr.shape[0])]
    trackedVals = [(x, y, c)
                   for i, x, y, c in zip(xrange(len(valVecList)), xVals, yVals, cVals)
                   if reverseMap[i] in trackers]
    trXVals, trYVals, trCVals = [list(x) for x in zip(*trackedVals)]
    fig2, clusterAx = plt.subplots()
    clusterAx.grid(True)
    clusterAx.set_title("LOS distribution characteristics by facility\n"
                        "colored by cluster")
    clusterAx.set_xlabel(xLabel)
    clusterAx.set_ylabel(yLabel)
    clusterAx.grid(True)
    makeScatterPlot(clusterAx, [(xVals, yVals, cVals, 'o', 20),
                                (xCtr, yCtr, cCtr, '^', 200),
                                (trXVals, trYVals, trCVals, '+', 200)
                                ])
    fig2.tight_layout()
    fig2.canvas.set_window_title("By Cluster")

    if annotateScatterPts:
        for abbrev, offset in indexDict.items():
            xy = (xVals[offset], yVals[offset])
            xytext = (xVals[offset]+0.00125, yVals[offset]+0.00125)
            clusterAx.annotate(abbrev, xy=xy, xytext=xytext)
            scatterAx.annotate(abbrev, xy=xy, xytext=xytext)

    fig3, histoAxes = plt.subplots(nrows=1, ncols=len(codeDict))
    if len(codeDict) == 1:
        histoAxes = [histoAxes]
    for category, i in codeDict.items():
        samples = []
        for abbrev, offset in indexDict.items():
            if code[offset] == i:
                samples.extend(losHistoDict[abbrev])
        vec = np.zeros(histoBins)
        if samples:
            for samp in samples:
                vec = resampleBand(samp, vec, range=histoRange)
        plotAsHistogram(vec, histoAxes[i], range=histoRange)
        allSampList = []
        rLow, rHigh = histoRange
        delta = (rHigh - rLow)/vec.shape[0]
        for step in xrange(vec.shape[0]):
            lo = rLow + step * delta
            hi = lo + delta
            allSampList.append((lo, hi, vec[step]/delta))
        allSampList.sort()
        losModel = LOSModelType()
        print ('%s: '%category),
        fitVec, accuracyMeasure, success = losModel.modelFit(allSampList,  # @UnusedVariable
                                                             funToMaximize=boundedLnLik)
        lowBandEdge = min([lo for lo, hi, ct in allSampList])  # @UnusedVariable
        highBandEdge = max([hi for lo, hi, ct in allSampList])  # @UnusedVariable
        curveScale = ((np.sum(vec) * (histoRange[1]-histoRange[0]))
                      / (histoBins*losModel.intervalCDF(lowBandEdge, highBandEdge, fitVec)))
        plotCurve(histoAxes[i], fitVec, losModel,
                  np.sum(vec)*((histoRange[1]-histoRange[0])/histoBins),
                  rng=histoRange, nbins=vec.shape[0])
        histoAxes[i].set_title(category)
        if category == 'HOSPITAL':
            print 'at src: %s %s' % (sum([ct for lo, hi, ct in allSampList]),curveScale)
            print 'First and last of %d bins: %s %s' % (len(allSampList), allSampList[0],
                                                        allSampList[-1])
            testSampL = allSampList[:]
            testVec = vec.copy()
    fig3.tight_layout()
    fig3.canvas.set_window_title("Category LOS Histograms")

    fig4, locAxes = plt.subplots(nrows=1, ncols=len(trackers))
    for i in xrange(len(trackers)):
        abbrev = trackers[i]
        vec = np.zeros(histoBins)
        for band in losHistoDict[abbrev]:
            vec = resampleBand(band, vec, range=histoRange)
        plotAsHistogram(vec, locAxes[i], range=histoRange)
        locAxes[i].set_title(abbrev)
        if abbrev in indexDict:
            fitParms = valVecList[indexDict[abbrev]]
            losModel = LOSModelType()

            lowBandEdge = min([a for a, b, c in losHistoDict[abbrev]])
            highBandEdge = max([b for a, b, c in losHistoDict[abbrev]])
            curveScale = ((np.sum(vec) * (histoRange[1] - histoRange[0]))
                      / (histoBins*losModel.intervalCDF(lowBandEdge, highBandEdge, fitVec)))

            plotCurve(locAxes[i], fitParms, losModel, curveScale,
                      #np.sum(vec)*((histoRange[1]-histoRange[0])/histoBins),
                      rng=histoRange, nbins=vec.shape[0])
        else:
            print 'No fit curve was calculated for %s' % abbrev
    fig4.canvas.set_window_title("Tracked Locations")

#     testFig, testAxes = plt.subplots()  # @UnusedVariable
# #    bandL = []
# #     for i in xrange(100):
# #         band = (float(i), float(i+1), 100.0*math.exp(-(float(i)+0.5)/100.0))
# #         bandL.append(band)
#     bandL = testSampL
#     losModel = TwoExponLOSModel()
#     fitVec, accuracyMeasure, success = losModel.modelFit(bandL,  # @UnusedVariable
#                                                          funToMaximize=boundedLnLik)
#     vec = testVec
# #     vec = np.zeros(histoBins)
# #     for band in bandL:
# #         vec = resampleBand(band, vec, range=histoRange)
#     plotAsHistogram(vec, testAxes, range=histoRange)
#     print 'sum is %s' % np.sum(vec)
#     print 'alternate sum is %s' % np.sum([ct for lo, hi, ct in bandL])
#     lowBandEdge = min([lo for lo, hi, ct in bandL])
#     highBandEdge = max([hi for lo, hi, ct in bandL])
#     curveScale = ((np.sum(vec) * (histoRange[1] - histoRange[0]))
#                   / (histoBins*losModel.intervalCDF(lowBandEdge, highBandEdge, fitVec)))
#     print 'BandEdges: %s %s' % (lowBandEdge, highBandEdge)
#     print 'scaling: %s %s %s' % (np.sum(vec), losModel.intervalCDF(lowBandEdge, highBandEdge, fitVec), curveScale)
#     print 'fitVec: %s' % fitVec
#     print 'First and last of %d bins: %s %s' % (len(bandL), bandL[0], bandL[-1])
#     plotCurve(testAxes, fitVec, losModel, curveScale,
#               rng=histoRange, nbins=histoBins)
#     testAxes.set_title('Test')

    plt.show()

############
# Main hook
############

if __name__ == "__main__":
    main()
