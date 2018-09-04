#!/usr/bin/env python
# encoding: utf-8
'''
los_model_tester -- Given a set of LOS samples and a LOS model, generate test statistics and graphs

los_model_tester is a command line program to test the LOS models for a set of facilities against
the LOS samples for those facilities

It's mainly just a command line utility.

@author:     Joel Welling

@copyright:  2018 CMU. All rights reserved.

@license:    license

@contact:    welling@psc.edu
@deffield    updated: Updated
'''

from __future__ import print_function

import six
import sys
import os
import re

from optparse import OptionParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import tools_util as tu
import stats
from notes_plotter import LOSPlotter, loadFacilityTypeConstants, findFacImplCategory
import pyrheautils as pu

__all__ = []
__version__ = 0.1
__date__ = '2018-09-04'
__updated__ = '2018-09-04'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

IMPL_TIER_DICT = {'HOSPITAL': 'HOSP', 'NURSINGHOME': 'NURSING',
                  'LTAC': 'LTAC', 'VSNF': 'VENT'}


def collectBarSamples(df):
    quantum = 1.0
    barWidth = 1.6 * quantum
    counts = df.groupby(['CODE', 'LOS']).size().reset_index(name='Count')
    return ((counts['LOS'].values - 0.5*quantum), counts['Count'].values, barWidth)


def main(argv=None):
    '''Command line options.'''

    program_name = os.path.basename(sys.argv[0])
    program_version = "v0.1"
    program_build_date = "%s" % __updated__

    program_version_string = '%%prog %s (%s)' % (program_version, program_build_date)
    program_usage = '''usage: %%prog [-f abbrev [-f abbrev [...]]] -s samples.csv rundesc.yaml'''
    program_longdesc = ''''''
    program_license = "Copyright 2018 user_name (organization_name)                                            \
                Licensed under the Apache License 2.0\nhttp://www.apache.org/licenses/LICENSE-2.0"

    if argv is None:
        argv = sys.argv[1:]
    if True:
#     try:
        # setup option parser
        parser = OptionParser(version=program_version_string, epilog=program_longdesc,
                              description=program_license)
        parser.add_option("-f", "--fac",
                          help="add facility to analyze (default is all tracked facilities)",
                          action='append')
        parser.add_option("-s", "--samples", help="csv file providing LOS samples (required)")
        parser.add_option("-v", "--verbose", dest="verbose", action="count",
                          help="set verbosity level [default: %default]")

        # process options
        (opts, args) = parser.parse_args(argv)
        if len(args) != 1:
            parser.error('A run description file is required\n')

        if not opts.samples:
            parser.error('You must use --fac to provide a csv file of LOS samples\n')

        if opts.verbose > 0:
            print("verbosity level = %d" % opts.verbose)

        # MAIN BODY #

        inputDict = tu.readModelInputs('../sim/twoyear_run_OC.yaml')
        pu.prepPathTranslations(inputDict)
        facDict = tu.getFacDict(inputDict)

        catNames = list(set([rec['category'] for rec in facDict.values()]))

        if 'facilitySelectors' in inputDict:
            facImplRules = [(re.compile(rule['category']), rule['implementation'])
                           for rule in inputDict['facilitySelectors']]
        else:
            facImplRules = [(re.compile(cat), cat)
                            for cat in catNames]  # an identity map
        implDir = pu.pathTranslate('$(IMPLDIR)')
        constDir = pu.pathTranslate('$(CONSTANTS)')
        catToImplDict = {cat: findFacImplCategory(implDir, facImplRules, cat)
                         for cat in catNames}

        testThese = inputDict['trackedFacilities'][:] if opts.fac is None else opts.fac[:]
        print('targets: %s' % testThese)

        sampDF = pd.read_csv(opts.samples)
        for fac in testThese:
            facSampDF = sampDF[sampDF['CODE'] == fac]
            nSamples = facSampDF['CODE'].count()
            if not nSamples:
                print('%s: no LOS samples in table' % fac)
            else:
                facCRV = stats.fullCRVFromPDFModel(facDict[fac]['losModel'])
                cat = facDict[fac]['category']
                print('%s: %d samples, mean %s median %s' % (fac, nSamples,
                                                             facCRV.mean(), facCRV.median()))
                constants = loadFacilityTypeConstants(catToImplDict[cat].lower(), implDir)
                fig, axes = plt.subplots(1, 1)
                lP = LOSPlotter(facDict[fac], constants, catToImplDict)
                bins, counts, barWidth = collectBarSamples(facSampDF)
                rects = axes.bar(bins, counts, width=barWidth,  # @UnusedVariable
                                        color='b')
                lP.plot(IMPL_TIER_DICT[catToImplDict[cat]], axes, int(max(bins)-min(bins)),
                        min(bins) - barWidth, max(bins) + barWidth, sum(counts))

                axes.set_title('%s (%d samples)' % (fac, nSamples))
        plt.show()

#     except Exception, e:
#         indent = len(program_name) * " "
#         sys.stderr.write(program_name + ": " + repr(e) + "\n")
#         sys.stderr.write(indent + "  for help use --help")
#         return 2


if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-v")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'los_model_tester_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())