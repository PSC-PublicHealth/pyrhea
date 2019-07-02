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

import logging
import types
import sys
import random
from functools import wraps
from math import fabs, log, exp
import unittest
import cStringIO

import numpy as np
from scipy.stats import kstest
from scipy.stats.distributions import rv_continuous, lognorm, expon, weibull_min
from scipy.optimize import minimize_scalar


logger = logging.getLogger(__name__)


class LogNormPlusExp(rv_continuous):
    """
    Defined to yield k*lognorm(sigma=s, scale=math.exp(mu)) + (1-k)*expon(scale=1/lmda)
    """
    def _pdf(self, x, s, mu, k, lmda):
        """s is the non-keyword parameter of lognorm, which corresponds to sigma.  """
        return ((lognorm.pdf(x, s, scale=np.exp(mu)) * k)
                + (expon.pdf(x, scale=(1.0 / lmda)) * (1.0 - k)))

    def _cdf(self, x, s, mu, k, lmda):
        return ((lognorm.cdf(x, s, scale=np.exp(mu)) * k)
                + (expon.cdf(x, scale=(1.0 / lmda)) * (1.0 - k)))

    def _argcheck(self, s, mu, k, lmda):
        return np.all(s > 0.0) and np.all(k >= 0.0) and np.all(k <= 1.0) and np.all(lmda > 0.0)

lognormplusexp = LogNormPlusExp(name='lognormplusexp', a=0.0)


class DoubleWeibull(rv_continuous):
    """
    defined to yield k*weibull(c=shape1, lmda=scale1) + (1-k)*weibull(c=shape2, lmda=scale2)
    """
    def _pdf(self, x, k, shape1, scale1, shape2, scale2):
        return ((weibull_min.pdf(x, shape1, scale=scale1) * k)
                + (weibull_min.pdf(x, shape2, scale=scale2) * (1.0-k)))

    def _cdf(self, x, k, shape1, scale1, shape2, scale2):
        return ((weibull_min.cdf(x, shape1, scale=scale1) * k)
                + (weibull_min.cdf(x, shape2, scale=scale2) * (1.0-k)))

    def _argcheck(self, k, shape1, scale1, shape2, scale2):
        return (np.all(k >= 0.0) and np.all(k<= 1.0) and np.all(shape1 > 0.0)
                and np.all(scale1 > 0.0) and np.all(shape2 > 0.0) and np.all(scale2 > 0.0))

doubleweibull = DoubleWeibull(name='doubleweibull', a=0.0)

class DoubleExpon(rv_continuous):
    """
    defined to yield k*expon(scale=1.0/lmda1) + (1-k)*expon(scale=1.0/lmda2)
    """
    def _argcheck(self, k, lmda1, lmda2):
        return (np.all(k >= 0.0) and np.all(k<= 1.0) and np.all(lmda1 > 0.0)
                and np.all(lmda2 > 0.0))
    def _pdf(self, x, k, lmda1, lmda2):
        submodels = [expon(scale=1.0/lmda1), expon(scale=1.0/lmda2)]
        return (k * submodels[0].pdf(x)) + ((1.0 - k) * submodels[1].pdf(x))

    def _cdf(self, x, k, lmda1, lmda2):
        submodels = [expon(scale=1.0/lmda1), expon(scale=1.0/lmda2)]
        return (k * submodels[0].cdf(x)) + ((1.0 - k) * submodels[1].cdf(x))

    def _logpdf(self, x, k, lmda1, lmda2):
        # pdf = k lmda1 exp(-lmda1 x) + (1-k)lmda2 exp(-lmda2 x)
        #     = k lmda1 exp(-lmda1 x)(1 + ((1-k)/k)(lmda2/lmda1)exp(lmda1-lmda2 x))
        return np.choose(k > 0.5,
                         [(np.log(1.0-k) + np.log(lmda2) - (lmda2 * x)
                           + np.log(1.0 + (((k/(1.0-k))*(lmda2/lmda2)*np.exp((lmda2 - lmda2) * x))))),
                          (np.log(k) + np.log(lmda1) - (lmda1 * x)
                           + np.log(1.0 + (((1.0-k)/k)*(lmda2/lmda1)*np.exp((lmda1 - lmda2) * x))))])

    def _rvs(self, k, lmda1, lmda2):
        size = self._size
        submodels = [expon(scale=1.0/lmda1), expon(scale=1.0/lmda2)]
        submodel_choices = np.random.choice(2, size, p=[k, 1.0-k])
        submodel_samples = [submodel.rvs(size=size)
                            for submodel in submodels]
        return np.choose(submodel_choices, submodel_samples)

doubleexpon = DoubleExpon(name='doubleexpon', a=0.0)


def _paramGrabberWrapper(mthd, keyword, update_mthd):
    @wraps(mthd)
    def wrapped_mthd(self, *args, **kwds):
        if keyword in kwds:
#             print('metaclass grabbed %s = %s from call to %s' % (keyword, kwds[keyword],
#                                                                  mthd))
            getattr(self, update_mthd)(kwds[keyword])
            del kwds[keyword]
        else:
            pass
        return mthd(self, *args, **kwds)
    return wrapped_mthd


class ParamGrabberMetaClass(type):
    def __new__(cls, name, parents, classDct):
        keyword = classDct['_meta_keyword_']
        update_mthd = classDct['_meta_update_mthd_']
        newClassDct = {}
        base_cls = parents[0]
        for entry in dir(base_cls):
            if entry[0] != '_':
                attr = getattr(base_cls, entry)
                if callable(attr):
                    attr = _paramGrabberWrapper(attr, keyword, update_mthd)
                    newClassDct[entry] = attr
        for entry, attr in classDct.items():
            if callable(attr) and entry[0] != '_':
                attr = _paramGrabberWrapper(attr, keyword, update_mthd)
            newClassDct[entry] = attr
        return type.__new__(cls, name, parents, newClassDct)


class Empirical(rv_continuous):
    """
    The real empirical distribution is discrete, but we define this one to be
    continuous, with the samples for a given day presumed to occur uniformly over
    the course of the day.  Thus the CDF linearly interpolates across the course
    of each day, and the PDF is a step function and constant within a day.  Samples
    drawn from this crv are presumed to correspond to real-world samples associated
    with the day which is the floor() of the sample time.
    """
    __metaclass__ = ParamGrabberMetaClass
    _meta_keyword_ = 'histogram_dict'
    _meta_update_mthd_ = '_histogram_update'

    def freeze(self, *args, **kwargs):
        rslt = super(Empirical, self).freeze(*args, **kwargs)
        rslt.dist._histogram_update(getattr(self, Empirical._meta_keyword_))
#         print('%s %s now has %s %s' % (type(rslt.dist), rslt.dist, Empirical._meta_keyword_,
#                                        getattr(rslt.dist, Empirical._meta_keyword_)))
        return rslt

    def _histogram_update(self, histo_dct):
        # print('%s %s just got update for %s' % (type(self), self, histo_dct))
        if histo_dct is None or not histo_dct:
            raise ValueError('Sample histogram is missing or empty')
        self.histogram_dict = histo_dct
        self.hist_min = min(histo_dct.keys())
        self.hist_max = max(histo_dct.keys())
        valL, ctL = zip(*histo_dct.items())
        self.hist_sampV = np.asarray(valL)
        self.hist_sampP = np.asarray(ctL, dtype=np.float)
        samp_count = self.hist_sampP.sum()
        if samp_count <= 0.0:
            raise ValueError('Sample histogram must have a positive number of samples')
        self.hist_sampP /= samp_count
        self.hist_pdfV = np.zeros(self.hist_max + 1 - self.hist_min)
        for k, v in self.histogram_dict.items():
            if v < 0:
                raise ValueError('Sample histogram entry has a negative number of counts')
            self.hist_pdfV[k - self.hist_min] = float(v)/samp_count
        self.hist_cdfV = np.cumsum(self.hist_pdfV)
        self.hist_cdfV[-1] = 1.0  # avoid rounding error
#         print 'hist_min: ', self.hist_min
#         print 'hist_max: ', self.hist_max
#         print 'hist_sampV: ', self.hist_sampV
#         print 'hist_sampP: ', self.hist_sampP
#         print 'hist_pdfV: ', self.hist_pdfV
#         print 'hist_cdfV: ', self.hist_cdfV

    def __init__(self, *args, **kwargs):
        super(Empirical, self).__init__(*args, **kwargs)
        self.histogram_dict = None
        self.hist_sampV = None
        self.hist_sampP = None
        self.hist_pdfV = None
        self.hist_cdfV = None
        self.hist_min = None  # also serves as base offset for the arrays
        self.hist_max = None
        self.hist_totcts = None

    def _argcheck(self, *args, **kwargs):
        return self.histogram_dict is not None

    def _pdf(self, x):
        idxV = np.floor(x).astype(np.int) - self.hist_min
        maskV = np.logical_and((idxV >= 0),(idxV <= self.hist_max - self.hist_min))
        idxV = np.maximum(idxV, 0)
        idxV = np.minimum(idxV, self.hist_max - self.hist_min)
        return np.choose(maskV, [0.0, self.hist_pdfV[idxV]])

    def _cdf(self, x):
        return np.interp(x,
                         np.linspace(self.hist_min, self.hist_max,
                                     self.hist_max + 1 - self.hist_min),
                         self.hist_cdfV)
 
    def _rvs(self, *args, **kwargs):
        return self._random_state.choice(self.hist_sampV, size=self._size,
                                         replace=True, p=self.hist_sampP)


empirical = Empirical(name='empirical', a=0.0)


class AlwaysZeroCRV(rv_continuous):
    """
    The pdf and cdf of this distribution are 0.0 for all finite x.  Useful when we want to
    turn off some feature which is driven by a pdf.
    """
    def _argcheck(self):
        return True

    def _pdf(self, x):
        return np.zeros_like(x)

    def _cdf(self, x):
        return np.zeros_like(x)

    def _logpdf(self, x):
        return np.full(x.shape, np.NINF)

    def _ppf(self, x):
        rslt = np.choose(np.greater(x, 0.0),
                         [np.full(x.shape, np.inf),
                          np.full(x.shape, 0.0)])
        return rslt

    def _rvs(self, k, lmda1, lmda2):
        size = self._size
        if size > 0:
            raise RuntimeError('Tried to take samples from the AlwaysZero distribution!')

alwayszerocrv = AlwaysZeroCRV(name='alwayszerocrv', a=0.0)


class CachedCDFGenerator:
    """
    This supports a common operation performed by Facilities- given a treatment interval (relative
    to a start date of zero), return a likelihood.  For example, this may be the likelihood that
    the patient will be discharged.  Since the interval bounds are integers, caching of the
    generated values is very effective.
    """
    def __init__(self, frozenCRV):
        self.frozenCRV = frozenCRV
        self.mdn = frozenCRV.median()
        self.cache = {}

    def intervalProb(self, start, end):
        key = (start, end)
        if key in self.cache:
            return self.cache[key]
        else:
            if end <= self.mdn:
                # both must be to the left of median
                sv = self.frozenCRV.cdf(start)
                cP = ((self.frozenCRV.cdf(end) - sv) / (1.0 - sv))

            elif start <= self.mdn:
                # spans median
                sv = self.frozenCRV.cdf(start)
                esf = self.frozenCRV.sf(end)
                cP = (1.0 - (sv + esf)) / (1.0 - sv)
            else:
                # both must be to the right of median
                ssf = self.frozenCRV.sf(start)
                esf = self.frozenCRV.sf(end)
                cP = (ssf - esf) / ssf
            self.cache[key] = cP
            return cP


class BayesTree(object):
    def __init__(self, v1, v2=None, p=1.0, tag=None):
        """
        Caution: neither v1 nor v2 should be a sequence type!
        """
        if isinstance(v1, BayesTree):
            v1TT = v1.tagTree
            v1 = v1.tree
        else:
            v1TT = None
        if isinstance(v2, BayesTree):
            v2TT = v2.tagTree
            v2 = v2.tree
        else:
            v2TT = None
        if p == 1.0:
            if tag is None and v2TT is None:
                # Trim this trivial level
                self.tree = v1
                self.tagTree = v1TT
            else:
                # Cannot trim the level without losing the tag
                self.tree = (p, v1, v2)
                self.tagTree = (tag, v1TT, v2TT)
        else:
            assert v2 is not None, 'BayesTree defined with transition to None'
            self.tree = (p, v1, v2)
            self.tagTree = (tag, v1TT, v2TT)

    @staticmethod
    def fromTuple(tpl, tag=None):
        """
        This convenience function expects a tuple of the form (p, v1, v2) where
        the components match the inputs to __init__().
        """
        if isinstance(tpl, types.TupleType):
            p, v1, v2 = tpl
            return BayesTree(v1, v2, p, tag=tag)
        else:
            return BayesTree(tpl, tag=tag)

    @staticmethod
    def fromLinearCDF(linearCDF, tol=0.00001, tag=None):
        """
        Given a linear CDF in list form, return an equivalent CDF in BayesTree form.

        The input CDF is of the form:

           [(p1, v1), (p2, v2), ...]

        such that the sum of p1, p2, etc is 1.0 to within tolerance tol.
        """
        lCDF = len(linearCDF)
        if lCDF == 1:
            p, r = linearCDF[0]
            assert fabs(p - 1.0) <= tol, 'CDF terms do not sum to 1 (%s)' % linearCDF
            return BayesTree(r, tag=tag)
        elif lCDF == 2:
            p1, r1 = linearCDF[0]
            p2, r2 = linearCDF[1]
            assert fabs(p1 + p2 - 1.0) <= tol, 'CDF terms do not sum to 1 (%s)' % linearCDF
            if p2 >= p1:
                return BayesTree(r1, r2, p1, tag=tag)
            else:
                return BayesTree(r2, r1, p2, tag=tag)
        else:
            linearCDF.sort()
            part1 = linearCDF[:lCDF/2]
            part2 = linearCDF[lCDF/2:]
            pivot = part1[-1][0]
            if pivot == 0.0:
                return BayesTree.fromLinearCDF(part2, tag=tag)
            else:
                w1 = sum([p for p, r in part1])
                w2 = sum([p for p, r in part2])
                assert fabs(w1 + w2 - 1.0) <= tol, 'CDF terms do not sum to 1 (%s)' % linearCDF
                return BayesTree(BayesTree.fromLinearCDF([(p / w1, r) for p, r in part1]),
                                 BayesTree.fromLinearCDF([(p / w2, r) for p, r in part2]),
                                 w1, tag=tag)

    @staticmethod
    def _innerDump(tree, tagTree, indent=0, ofile=sys.stdout):
        if isinstance(tree, types.TupleType):
            prob, path1, path2 = tree
            if tagTree:
                topTag, tagPath1, tagPath2 = tagTree
                if topTag:
                    ofile.write('%s(%s tag=%s\n' % (' '*indent, prob, topTag))
                else:
                    ofile.write('%s(%s\n' % (' '*indent, prob))
                BayesTree._innerDump(path1, tagPath1, indent+4, ofile=ofile)
                BayesTree._innerDump(path2, tagPath2, indent+4, ofile=ofile)
            else:
                ofile.write('%s(%s\n' % (' '*indent, prob))
                BayesTree._innerDump(path1, None, indent+4, ofile=ofile)
                BayesTree._innerDump(path2, None, indent+4, ofile=ofile)
            ofile.write('%s)\n' % (' '*indent))
        else:
            ofile.write('%s%s\n' % (' '*indent, tree))

    def dump(self, ofile=sys.stdout):
        self._innerDump(self.tree, self.tagTree, ofile=ofile)

    @staticmethod
    def _innerTraverse(tree, rng):
        # PatientAgent._printBayesTree(tree)
        while True:
            if isinstance(tree, types.TupleType):
                if isinstance(tree[0], (types.FunctionType, types.MethodType)):
#                     import pdb
#                     pdb.set_trace()
                    raise RuntimeError('Cannot traverse BayesTree; encountered unevaluated %s'
                                       % repr(tree[0]))
                elif rng() <= tree[0]:
                    tree = tree[1]
                else:
                    tree = tree[2]
            else:
                # print 'walk returning %s' % str(tree)
                return tree

    def traverse(self, randomNumberGenerator=None):
        """Traverse the tree, returning the value from one randomly chosen leaf.

        If randomNumberGenerator is provided, randomNumberGenerator.random() will be used
        to generate the p-values used for the traversal.
        """
        # self.dump()
        if randomNumberGenerator is None:
            rng = random.random
        else:
            rng = randomNumberGenerator.random
        return self._innerTraverse(self.tree, rng)

    @staticmethod
    def _innerFindTag(tree, tagTree, target):
        if isinstance(tree, types.TupleType):
            assert isinstance(tagTree, types.TupleType), ('Ran out of tree looking for tag %s'
                                                          % target)
            prob, lTree, rTree = tree  # @UnusedVariable
            topTag, lTag, rTag = tagTree
            if topTag == target:
                return BayesTree.fromTuple(tree, tag=tagTree)
            else:
                lRslt = BayesTree._innerFindTag(lTree, lTag, target)
                if lRslt:
                    return lRslt
                else:
                    return BayesTree._innerFindTag(rTree, rTag, target)

        else:
            assert not isinstance(tagTree, types.TupleType), ('Ran out of tags looking for tag %s'
                                                              % target)
            if tagTree == target:
                return BayesTree.fromTuple(tree, tag=tagTree)
            else:
                return None

    def findTag(self, target):
        return self._innerFindTag(self.tree, self.tagTree, target)

    def replaceSubtree(self, target, replTree):
        if isinstance(self.tree, types.TupleType):
            prob, lTree, rTree = self.tree
            topTag, lTag, rTag = self.tagTree
            if topTag == target:
                return replTree
            else:
                lBT = BayesTree.fromTuple(lTree, tag=lTag).replaceSubtree(target, replTree)
                if lBT.tree is not lTree:
                    return BayesTree(lBT,
                                     BayesTree.fromTuple(rTree, tag=rTag),
                                     prob)
                else:
                    return BayesTree(BayesTree.fromTuple(lTree, tag=lTag),
                                     BayesTree.fromTuple(rTree, tag=rTag)).replaceSubtree(target,
                                                                                          replTree)
        else:
            if self.tagTree == target:
                return replTree
            else:
                return BayesTree.fromTuple(self.tree, tag=self.tagTree)

    def getParts(self):
        """
        This method returns a tuple (leftSubTree, rightSubTree, p, topTag).
        Note that if this tree is already a leaf, it will return an essentially identical tree
        in place of its left subtree and 'None' for the right subtree.
        """
        if isinstance(self.tree, types.TupleType):
            prob, lTree, rTree = self.tree
            topTag, lTag, rTag = self.tagTree
            return (BayesTree.fromTuple(lTree, tag=lTag),
                    BayesTree.fromTuple(rTree, tag=rTag),
                    prob, topTag)
        else:
            return (BayesTree.fromTuple(self.tree, tag=self.tagTree),
                    None,
                    1.0, None)

    def copy(self):
        """
        Returns a shallow copy of the BayesTree
        """
        lTree, rTree, prob, topTag = self.getParts()
        if rTree is None:
            return lTree
        else:
            return BayesTree(lTree.copy(), rTree.copy(), prob, tag=topTag)

def fullLogNormCRVFromMean(mean, sigma):
    """
    Returns a scipy.stats.distributions frozen lognorm continuous random variable from a mean
    and sigma.
    """
    mu = log(mean) - (0.5 * sigma * sigma)
    return lognorm(sigma, scale=np.exp(mu), loc=0.0)


def fullCRVFromPDFModel(pdfModel):
    """
    Returns a scipy.stats.distributions frozen continuous random variable from the
    data in a PDF model.  The PDF model has the form:
        pdfModel = {'pdf': 'some descriptive text string'
                    'parms': iterable of floats }
    The algorithm is a simple case statement with matches for known pdf strings.
    """
    modelStr = pdfModel['pdf'].strip('"\'')
    if modelStr == 'lognorm(mu=$0,sigma=$1)':
        mu, sigma = pdfModel['parms']
        return lognorm(sigma, scale=exp(mu), loc=0.0)
    elif modelStr == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)':
        k, mu, sigma, lmda = pdfModel['parms']
        return lognormplusexp(s=sigma, mu=mu, k=k, lmda=lmda)
    elif modelStr == 'expon(lambda=$0)':
        lmda = pdfModel['parms'][0]
        return expon(scale=1.0/lmda)
    elif modelStr == '$0*weibull(k=$1, lmda=$2)+(1-$0)*weibull(k=$3, lmda=$4)':
        k, shape1, scale1, shape2, scale2 = pdfModel['parms']
        return doubleweibull(k, shape1, scale1, shape2, scale2)
    elif modelStr == '$0*expon(lambda=$1)+(1-$0)*expon(lambda=$2)':
        k, lmda1, lmda2 = pdfModel['parms']
        return doubleexpon(k, lmda1, lmda2)
    elif modelStr == 'alwayszero()':
        return alwayszerocrv()
    elif modelStr == 'empirical(loc=$0,scale=$1)':
        assert 'sampleHistogram' in pdfModel, 'empirical pdfModel requires a sampleHistogram entry'
        loc, scale = pdfModel['parms']
        return empirical(loc=loc, scale=scale, histogram_dict=pdfModel['histo_dict'])
    else:
        raise RuntimeError('Unknown LOS model %s' % pdfModel['pdf'])

def jointSamples(nSamps, aCRV, bCRV, frac):
    """
    Generate nSamps samples from the joint distribution frac*(aCRV + bCRV) + (1.0-frac)*aCRV

    This is useful in decomposing mixed populations like hospital patients
    """
    aSamps = aCRV.rvs(nSamps)
    bSamps = bCRV.rvs(nSamps)
    choices = np.random.choice(2, nSamps, p=[frac, 1.0 - frac])
    return np.choose(choices, [aSamps + bSamps, aSamps])

def correctLogNormCRVForComponent(meanLOSHosp, lowBound, highBound,
                                  hospCRV, fracICU, icuCRV,
                                  nSamples=1000000):
    def func(val, nSamps, meanLOSHosp, hospCRV, fracICU, icuCRV):
        guessCRV = fullLogNormCRVFromMean(meanLOSHosp, val)
        samps = jointSamples(nSamps, guessCRV, icuCRV, fracICU)
        rslt = kstest(samps, hospCRV.cdf).statistic
        #print '%s -> %s' % (val, rslt)
        return rslt
    rslt = minimize_scalar(func, [lowBound, highBound],
                           args=(nSamples, meanLOSHosp, hospCRV, fracICU, icuCRV))
    assert rslt.success, 'Minimization failed; result was %s' % rslt
    return rslt.x

def createMapFun(parms):
    s = parms['s']
    mu = parms['mu']
    k = parms['k']
    lmda = parms['lmda']
    print 'testing %s' % parms

    def testfun(xVec):
        return ((lognorm.pdf(xVec, s, scale=np.exp(mu)) * k)
                + expon.pdf(xVec, scale=(1.0 / lmda)) * (1.0 - k))

    return testfun


_lognormplusexp_test_parmsets = [{'s': 2.0, 'mu': 2.0, 'k': 0.0, 'lmda': 3.0},
                                 {'s': 1.5, 'mu': 1.0, 'k': 0.25, 'lmda': 2.0},
                                 {'s': 1.0, 'mu': 2.0, 'k': 0.5, 'lmda': 1.0},
                                 {'s': 0.5, 'mu': 1.0, 'k': 0.75, 'lmda': 0.5},
                                 {'s': 0.3, 'mu': 2.0, 'k': 1.0, 'lmda': 0.5},
                                 ]

_bayestree_test_dump_string = """\
(0.7
    no change
    (0.3
        (0.333333333333
            line 1
            line 2
        )
        (0.285714285714
            line 4
            line 3
        )
    )
)

"""

_bayestree_test_traversal_dict = {'line 1': 280, 'line 2': 565, 'line 3': 1518,
                                  'line 4': 612, 'no change': 7025}


_empirical_test_dict = {1:5, 2:3, 4:3, 5:6, 18:1}

_empirical_test_dict2 = {12:7}

def main():
    "This is a simple test routine which takes kvp files as arguments"

    import matplotlib.pyplot as plt

    # Test empirical CRV
    print 'testing empirical CRV'
    TestStats('test_stats_empirical').debug()

    # Test CachedCDFGenerator
    print 'testing CachedCDFGenerator'
    TestStats('test_stats_cachedcdfgenerator').debug()

    # Test BayesTree
    print 'testing BayesTree constructors and dump'
    TestStats('test_stats_bayestree_dump').debug()

    print 'testing BayesTree traversal'
    TestStats('test_stats_bayestree_traversal').debug()

    # Test the lognormplusexp distribution
    print 'Plotting implementation of the lognormplusexp distribution'
    print '(this takes a minute to generate samples)'
    fig, ax = plt.subplots(len(_lognormplusexp_test_parmsets), 1)
    for ind, parms in enumerate(_lognormplusexp_test_parmsets):
        x = np.linspace(lognormplusexp.ppf(0.01, **parms),
                        lognormplusexp.ppf(0.99, **parms), 100)
        ax[ind].plot(x, lognormplusexp.pdf(x, **parms),
                     'r-', lw=5, alpha=0.6, label='lognormplusexp pdf')
        rv = lognormplusexp(**parms)
        ax[ind].plot(x, rv.pdf(x), 'b-', lw=2, label='frozen pdf')
        testFun = createMapFun(parms)
        ax[ind].plot(x, testFun(x), 'k+', lw=3, label='analytic')
        r = lognormplusexp.rvs(size=10000, **parms)
        ax[ind].hist(r, 100, density=True, histtype='stepfilled', alpha=0.4)
        ax[ind].legend(loc='best', frameon=False)
        ax[ind].set_title(str(parms))

        samplePts = [0.001, 0.5, 0.999]
        vals = lognormplusexp.ppf(samplePts, **parms)
        assert np.allclose(samplePts, lognormplusexp.cdf(vals, **parms)), \
            'allclose test failed for parms %s' % str(parms)

    fig.tight_layout()
    plt.show()


class TestStats(unittest.TestCase):
    def test_stats_empirical(self):
        # test missing or bad parameters
        rslt = empirical.pdf([12.3, 4.7, 0.2])
        # The best implementation so far only detects this at _argcheck, returning NaNs.
        self.assertTrue(np.all(np.isnan(rslt)))
        with self.assertRaises(ValueError):
            ignored = empirical.pdf([12.3, 4.7, 0.2], histogram_dict={})  # @UnusedVariable
        with self.assertRaises(ValueError):
            ignored = empirical.pdf([12.3, 4.7, 0.2],    # @UnusedVariable
                                    histogram_dict={1:1, 2:-1, 3:1})
        # unscaled pdf
        rslt = empirical.pdf(np.asarray([12.3, 4.7, 0.2]), histogram_dict=_empirical_test_dict)
        self.assertTrue(np.allclose(rslt, np.array([0., 0.16666667, 0.0])))
        # scaled pdf
        rslt = empirical.pdf(np.asarray([12.3, 4.7, 0.2]), histogram_dict=_empirical_test_dict,
                             loc=7., scale=3.)
        self.assertTrue(np.allclose(rslt, np.array([0.09259259, 0.0, 0.0])))
        # cdf
        rslt = empirical.cdf(np.asarray([12.3, 4.7, 0.2]), histogram_dict=_empirical_test_dict)
        self.assertTrue(np.allclose(rslt, np.array([0.94444444, 0.84444444, 0.27777778])))
        # create and test frozen crv
        crv = empirical(histogram_dict=_empirical_test_dict2, loc=7., scale=3.)
        rslt = crv.pdf(np.asarray([12.3, 4.7, 0.2]))
        self.assertTrue(np.allclose(rslt, np.array([0.0, 0.0, 0.0])))
        crv = empirical(histogram_dict=_empirical_test_dict)
        rslt = crv.pdf(np.asarray([12.3, 4.7, 0.2]))
        self.assertTrue(np.allclose(rslt, np.array([0.0,  0.16666667, 0.0])))
        # generate some samples
        rS = np.random.RandomState(1234)
        drawV = crv.rvs(size=10, random_state=rS)
        self.assertTrue(np.all(drawV == np.array([1, 4, 2, 5, 5, 1, 1, 5, 5, 5], dtype=np.int)))
        drawV = empirical(histogram_dict=_empirical_test_dict).rvs(size=12, random_state=rS)
        self.assertTrue(np.all(drawV == np.array([[2, 4, 5, 5, 2, 4, 4, 1, 5, 5, 2, 4]],
                                                 dtype=np.int)))

    def test_stats_lognormplusexp(self):
        for parms in _lognormplusexp_test_parmsets:
            samplePts = [0.001, 0.5, 0.999]
            vals = lognormplusexp.ppf(samplePts, **parms)
            self.assertTrue(np.allclose(samplePts, lognormplusexp.cdf(vals, **parms)),
                            msg='allclose test failed for parms %s' % str(parms))

    def test_stats_cachedcdfgenerator(self):
        cCG = CachedCDFGenerator(expon(0.7))
        v1 = cCG.intervalProb(1.2, 3.4)
        v2 = cCG.intervalProb(0.1, 2.3)
        v3 = cCG.intervalProb(1.2, 3.4)

        self.assertTrue(v1 == v3)
        self.assertTrue(v1 is v3, msg='value was not cached')
        self.assertTrue(v2 == (expon.cdf(2.3, 0.7)
                               - expon.cdf(0.1, 0.7))/(1.0 - expon.cdf(0.1, 0.7)))

    def test_stats_bayestree_dump(self):
        changeProb = 0.7
        with self.assertRaises(AssertionError):
            bT = BayesTree('no change',  # @UnusedVariable
                           BayesTree.fromLinearCDF([(0.1, 'line 1'),
                                                    (0.2, 'line 2'),
                                                    (0.5, 'line 3'),
                                                    (0.3, 'line 4')]),
                           changeProb)

        with self.assertRaises(AssertionError):
            bT = BayesTree('no change',  # @UnusedVariable
                           BayesTree.fromLinearCDF([(0.1, 'line 1'),
                                                    (0.2, 'line 2'),
                                                    (0.5, 'line 3'),
                                                    (0.1, 'line 4'),
                                                    (0.2, 'line 5')]),
                           changeProb)

        bT = BayesTree('no change',
                       BayesTree.fromLinearCDF([(0.1, 'line 1'),
                                                (0.2, 'line 2'),
                                                (0.5, 'line 3'),
                                                (0.2, 'line 4')]),
                       changeProb)
        fakeF = cStringIO.StringIO()
        bT.dump(fakeF)
        echoF = cStringIO.StringIO(fakeF.getvalue())
        fakeTest = cStringIO.StringIO(_bayestree_test_dump_string)
        for r1, r2 in zip(echoF.readlines(), fakeTest.readlines()):
            self.assertTrue(r1 == r2, 'dump readback: <%s> vs. <%s>' % (r1, r2))

    def test_stats_bayestree_traversal(self):
        bT = BayesTree('no change',
                       BayesTree.fromLinearCDF([(0.1, 'line 1'),
                                                (0.2, 'line 2'),
                                                (0.5, 'line 3'),
                                                (0.2, 'line 4')]),
                       0.7)

        rNG = random.Random()
        rNG.seed(1234)
        samps = {}
        for i in xrange(10000):  # @UnusedVariable
            val = bT.traverse(randomNumberGenerator=rNG)
            if val in samps:
                samps[val] += 1
            else:
                samps[val] = 1
        self.assertTrue(samps == _bayestree_test_traversal_dict, 'Traversal samples do not match')

############
# Main hook
############

if __name__ == "__main__":
    main()
