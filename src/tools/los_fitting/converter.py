#!/usr/bin/env python
# encoding: utf-8
'''
los_fitting.converter -- convert the output of nh-fitting-2weibull.R to 'standard' csv form

los_fitting.converter is a command line program which expects to find the outputs of a particular
R script in the current working directory.  It converts those outputs to CSV form, matching the
format produced by the los_model_fits family of programs.  The input filenames are fixed to
match the outputs of the R script.

It defines a main.  That's about it.

@author:     Joel Wellign

@copyright:  2018 PSC. All rights reserved.

@license:    license

@contact:    welling@psc.edu
@deffield    updated: Updated
'''

from __future__ import print_function
import six
import sys
import os
import six

from optparse import OptionParser
import pandas as pd
import numpy as np
from scipy.stats import weibull_min

import phacsl.utils.formats.csv_tools as csv_tools

__all__ = []
__version__ = 0.1
__date__ = '2018-09-01'
__updated__ = '2018-09-01'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

def parseDistParams(fn):
    """Despite the name, these aren't really CSVs"""
    rslt = {}
    idx = None
    try:
        with open(fn, 'rU') as f:
            for idx, line in enumerate(f):
                words = line.split()
                if idx == 0:
                    assert words == ['summary', 'of', 'weibullRMM_SEM', 'object:']
                elif idx == 1:
                    assert words == ['comp', '1', 'comp', '2']
                elif idx in [2, 3, 4]:
                    assert len(words) == 3
                    key, val1, val2 = words
                    val1 = float(val1)
                    val2 = float(val2)
                    rslt[key] = (val1, val2)
                elif idx == 5:
                    assert words[:3] == ['loglik', 'at', 'estimate:']
                    rslt['loglik'] = float(words[3])
                elif idx == 6:
                    assert words[1:] == ['%', 'of', 'the', 'data', 'right', 'censored']
                    if float(words[0]) >= 2.0:
                        sys.exit('%s was generated with too low a LOS censorship threshold'
                                 % fn)
                elif idx == 7:
                    assert words[0] == 'nobs:'
                    rslt['nobs'] = int(words[1])
                else:
                    sys.exit('%s has an unexpected format' % fn)
    except AssertionError:
        sys.exit('%s line %s has an unexpected format' % (fn, idx))
    assert sum(rslt['lambda']) == 1.0, 'lambda values in %s do not sum to 1.0' % fn

    # By convention we want the low-mean component first, and lambda should be the weight
    # of the low-mean component
    orderL = []
    for idx, (shape, scale) in enumerate(zip(rslt['shape'], rslt['scale'])):
        orderL.append((weibull_min(shape, scale=scale).mean(), idx))
    orderL.sort()
    orderL = [b for a, b in orderL]  # @UnusedVariable
    orderedRslt = rslt.copy()
    for key in ['shape', 'scale']:
        orderedRslt[key] = [rslt[key][idx] for idx in orderL]
    orderedRslt['lambda'] = rslt['lambda'][orderL[0]]
    print(orderedRslt)
    return orderedRslt

def main(argv=None):
    '''Command line options.'''

    program_name = os.path.basename(sys.argv[0])
    program_version = "v0.1"
    program_build_date = "%s" % __updated__

    program_version_string = '%%prog %s (%s)' % (program_version, program_build_date)
    #program_usage = '''usage: spam two eggs''' # optional - will be autogenerated by optparse
    program_longdesc = '''''' # optional - give further explanation about what the program does
    program_license = "Copyright 2018 (CMU)                                            \
                Licensed under the Apache License 2.0\nhttp://www.apache.org/licenses/LICENSE-2.0"

    if argv is None:
        argv = sys.argv[1:]
    try:
        # setup option parser
        parser = OptionParser(version=program_version_string, epilog=program_longdesc,
                              description=program_license)
        parser.add_option("-o", "--out", dest="outfile", help="set output path [default: %default]",
                          metavar="FILE")
        parser.add_option("-v", "--verbose", dest="verbose", action="count",
                          help="set verbosity level [default: %default]")

        # set defaults
        parser.set_defaults(outfile="./los_model_fits_2weibull.csv")

        # process options
        (opts, args) = parser.parse_args(argv)
        if args:
            parser.error('No arguments, only options')

        if opts.verbose > 0:
            print("verbosity level = %d" % opts.verbose)
        if opts.outfile:
            print("outfile = %s" % opts.outfile)

        # MAIN BODY #
        pFKDF = pd.read_csv('per-facility-k.csv')
        print(pFKDF.columns)

        dpDct = parseDistParams('dist-params.csv')
        extwDct = parseDistParams('extw-dist-params.csv')

        df = pd.DataFrame(columns=['abbrev', 'k', 'shape1', 'scale1', 'shape2', 'scale2',
                                   'lnLikPerSample', 'nsamples'])
        df = df.append({'abbrev': 'EXTW', 'k': extwDct['lambda'],
                        'shape1': extwDct['shape'][0], 'scale1':extwDct['scale'][0],
                        'shape2': extwDct['shape'][1], 'scale2':extwDct['scale'][0],
                        'lnLikPerSample': extwDct['loglik']/extwDct['nobs'],
                        'nsamples': extwDct['nobs'], 'notes': 'fit separately'
                        }, ignore_index=True)
        newDF = pFKDF[['code', 'k', 'n.obs']].rename(index=str, columns={'code':'abbrev',
                                                                         'n.obs':'nsamples'})
        newDF['shape1'] = dpDct['shape'][0]
        newDF['scale1'] = dpDct['scale'][0]
        newDF['shape2'] = dpDct['shape'][1]
        newDF['scale2'] = dpDct['scale'][1]
        newDF['notes'] = ''
        newDF['k'] = 1.0 - newDF['k']
        df = df.append(newDF, ignore_index=True)
        df = df.sort_values(by=['abbrev']).reset_index(drop=True)
        print(df)
        df.to_csv(opts.outfile, index=False)

    except TypeError, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-v")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'los_fitting.converter_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())