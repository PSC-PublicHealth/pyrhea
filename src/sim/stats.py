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

import types
import sys
import random
from scipy.stats.distributions import rv_continuous, lognorm, expon
from math import fabs, log, exp
import logging

import unittest
import cStringIO
import numpy as np
import matplotlib.pyplot as plt

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


class CachedCDFGenerator:
    """
    This supports a common operation performed by Facilities- given a treatment interval (relative
    to a start date of zero), return a likelihood.  For example, this may be the likelihood that
    the patient will be discharged.  Since the interval bounds are integers, caching of the
    generated values is very effective.
    """
    def __init__(self, frozenCRV):
        self.frozenCRV = frozenCRV
        self.cache = {}

    def intervalProb(self, start, end):
        key = (start, end)
        if key in self.cache:
            return self.cache[key]
        else:
            sv = self.frozenCRV.cdf(start)
            cP = ((self.frozenCRV.cdf(end) - sv)
                  / (1.0 - sv))
            self.cache[key] = cP
            return cP


class BayesTree(object):
    def __init__(self, v1, v2=None, p=1.0):
        """
        Caution: neither v1 nor v2 should be a sequence type!
        """
        if isinstance(v1, BayesTree):
            v1 = v1.tree
        if isinstance(v2, BayesTree):
            v2 = v2.tree
        if p == 1.0:
            self.tree = v1
        else:
            assert v2 is not None, 'BayesTree defined with transition to None'
            self.tree = (p, v1, v2)

    @staticmethod
    def fromTuple(tpl):
        """
        This convenience function expects a tuple of the form (p, v1, v2) where
        the components match the inputs to __init__().
        """
        p, v1, v2 = tpl
        return BayesTree(v1, v2, p)

    @staticmethod
    def fromLinearCDF(linearCDF, tol=0.00001):
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
            return BayesTree(r)
        elif lCDF == 2:
            p1, r1 = linearCDF[0]
            p2, r2 = linearCDF[1]
            assert fabs(p1 + p2 - 1.0) <= tol, 'CDF terms do not sum to 1 (%s)' % linearCDF
            if p2 >= p1:
                return BayesTree(r1, r2, p1)
            else:
                return BayesTree(r2, r1, p2)
        else:
            linearCDF.sort()
            part1 = linearCDF[:lCDF/2]
            part2 = linearCDF[lCDF/2:]
            pivot = part1[-1][0]
            if pivot == 0.0:
                return BayesTree.fromLinearCDF(part2)
            else:
                w1 = sum([p for p, r in part1])
                w2 = sum([p for p, r in part2])
                assert fabs(w1 + w2 - 1.0) <= tol, 'CDF terms do not sum to 1 (%s)' % linearCDF
                return BayesTree(BayesTree.fromLinearCDF([(p / w1, r) for p, r in part1]),
                                 BayesTree.fromLinearCDF([(p / w2, r) for p, r in part2]),
                                 w1)

    @staticmethod
    def _innerDump(tree, indent=0, ofile=sys.stdout):
        if isinstance(tree, types.TupleType):
            prob, path1, path2 = tree
            ofile.write('%s(%f\n' % (' '*indent, prob))
            BayesTree._innerDump(path1, indent+4, ofile=ofile)
            BayesTree._innerDump(path2, indent+4, ofile=ofile)
            ofile.write('%s)\n' % (' '*indent))
        else:
            ofile.write('%s%s\n' % (' '*indent, tree))

    def dump(self, ofile=sys.stdout):
        self._innerDump(self.tree, ofile=ofile)

    @staticmethod
    def _innerTraverse(tree, rng):
        # PatientAgent._printBayesTree(tree)
        while True:
            if isinstance(tree, types.TupleType):
                if rng() <= tree[0]:
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


def fullLogNormCRVFromMean(mean, sigma):
    """
    Returns a scipy.stats.distributions frozen lognorm continuous random variable from a mean
    and sigma.
    """
    mu = log(mean) - (0.5 * sigma * sigma)
    return lognorm(sigma, scale=exp(mu), loc=0.0)


def fullCRVFromPDFModel(pdfModel):
    """
    Returns a scipy.stats.distributions frozen continuous random variable from the
    data in a PDF model.  The PDF model has the form:
        pdfModel = {'pdf': 'some descriptive text string'
                    'parms': iterable of floats }
    The algorithm is a simple case statement with matches for known pdf strings.
    """
    if pdfModel['pdf'] == 'lognorm(mu=$0,sigma=$1)':
        mu, sigma = pdfModel['parms']
        return lognorm(sigma, scale=exp(mu), loc=0.0)
    elif pdfModel['pdf'] == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)':
        k, mu, sigma, lmda = pdfModel['parms']
        return lognormplusexp(s=sigma, mu=mu, k=k, lmda=lmda)
    elif pdfModel['pdf'] == 'expon(lambda=$0)':
        lmda = pdfModel['parms'][0]
        return expon(scale=1.0/lmda)
    else:
        raise RuntimeError('Unknown LOS model %s' % pdfModel['pdf'])


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
(0.700000
    no change
    (0.300000
        (0.333333
            line 1
            line 2
        )
        (0.285714
            line 4
            line 3
        )
    )
)

"""

_bayestree_test_traversal_dict = {'line 1': 280, 'line 2': 565, 'line 3': 1518,
                                  'line 4': 612, 'no change': 7025}


def main():
    "This is a simple test routine which takes kvp files as arguments"

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
        ax[ind].hist(r, 100, normed=True, histtype='stepfilled', alpha=0.4)
        ax[ind].legend(loc='best', frameon=False)
        ax[ind].set_title(str(parms))

        samplePts = [0.001, 0.5, 0.999]
        vals = lognormplusexp.ppf(samplePts, **parms)
        assert np.allclose(samplePts, lognormplusexp.cdf(vals, **parms)), \
            'allclose test failed for parms %s' % str(parms)

    fig.tight_layout()
    plt.show()


class TestStats(unittest.TestCase):
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
