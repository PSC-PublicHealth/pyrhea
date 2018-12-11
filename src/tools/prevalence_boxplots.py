#!/usr/bin/env python
# encoding: utf-8
'''
prevalence_boxplots -- This module generates boxplots of prevalence values from collections of test runs.

prevalence_boxplots is a tool that generates boxplots of prevalence values from collections of test runs.
A typical usage command would be something like:

  ./prevalence_boxplots.py -s $workdir/cre_prev_run_\*.mpk --glob -t $workdir/expected.pkl -y $workdir/taumod_config.yaml

where $workdir is the directory containing run inputs and results.  taumod_config is used to supply the
sampling days to be used.

Output is a set of files prevalence_TIER_nn.png (where TIER is for example HOSP and nn is an integer)
and prevalence_pooled_over_tiers.png .

It defines facBoxPlots, tierBoxPlot

@author:     welling

@copyright:  2018 Pittsubrgh Supercomputing Center (PSC). All rights reserved.

@license:    AGPL-3.0

@contact:    welling@psc.edu
@deffield    updated: Updated
'''

###################################################################################
# Copyright   2018, Pittsburgh Supercomputing Center (PSC).  All Rights Reserved. #
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


import sys
import os
import glob
import yaml
from collections import defaultdict

import cPickle as pickle
from optparse import OptionParser
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Hack for no display w/ pyplot. Change pyplot backend *before* loading it.
import matplotlib.pyplot as plt # use plt.savefig(<path>) rather than plt.show()

__all__ = []
__version__ = 0.1
__date__ = '2018-05-02'
__updated__ = '2018-12-11'

DEBUG = False
VERBOSE = False
TESTRUN = 0
PROFILE = 0

BOXES_PER_FIG = 10

SHOW_OPTIONS_DICT = {'prevalence': ('prev_sample', 'prevalence', 'prevalence'),
                     'colonizations': ('newColonized_sample', 'colonizations',
                                       'total colonizations per day'),
                     'population': ('TOTAL_sample', 'population', 'population'),
                     'arrivecolonized': ('arrive_colonized_sample', 'arrivecolonized',
                                         'fraction arriving colonized'),
                     'incidence': ('incidence_sample', 'incidence',
                                   'per-patient incidence'),
                     'arrivals': ('arrivals_sample', 'arrivals', 'arrivals')
                     }

def expandGlobbedList(pathList):
    """
    Given a list of path strings, expand with globbing
    """
    newPathList = []
    for fName in pathList:
        #print '%s yields %s' % (fName, glob.glob(fName))
        newPathList += glob.glob(fName)
    return newPathList


def buildLbl(tier, figlbl, longlbl):
    if longlbl is None:
        if figlbl is None:
            return '%s' % tier
        else:
            return '%s %s' % (figlbl, tier)
    else:
        if figlbl is None:
            return '%s %s' % (tier, longlbl)
        else:
            return '%s %s %s' % (figlbl, tier, longlbl)


def facBoxPlots(sampDF, targetD, tier, logy=False, label='prevalence', colKey='prev_sample',
                figlbl=None, longlbl=None):
    tierDF = sampDF[sampDF['tier']==tier]
    tierGps = tierDF.groupby('abbrev')
    nBlocks = (len(tierGps) / BOXES_PER_FIG) + 1
    tupleL = [tpl for tpl in tierGps]
    fig, axes = plt.subplots(1, 1)  # @UnusedVariable
    if VERBOSE:
        print 'starting tier %s' % tier
    for block in range(nBlocks):
        blockNum = block + 1
        tupleBlockL = tupleL[block * BOXES_PER_FIG : (block + 1) * BOXES_PER_FIG]
        if not tupleBlockL:
            break  # annoying edge case
        labelL = []
        sampL = []
        markerXL = []
        markerYL = []
        for idx, (key, subDF) in enumerate(tupleBlockL):
            if DEBUG:
                print 'fac %s tier %s' % (key, tier)
            labelL.append(key)
            sampL.append(subDF[colKey].dropna())
            if targetD is not None and (key, tier) in targetD:
                markerXL.append(1.0 + idx)
                markerYL.append(targetD[(key, tier)])
        if DEBUG:
            print 'plotting %s' % str(labelL)
        axes.tick_params(axis='x', labelsize='small', labelrotation=45.0)
        axes.boxplot(sampL, labels=labelL)
        if markerYL:
            axes.plot(markerXL, markerYL, 'D')

        if nBlocks > 1:
            axes.set_title('%s %d' % (buildLbl(tier, figlbl, longlbl), blockNum))
        else:
            axes.set_title(buildLbl(tier, figlbl, longlbl))
        if logy:
            axes.set_yscale('log')
        plt.tight_layout()
        plt.savefig('%s_%s_%02d.png' % (label, tier, block))
        plt.cla()  # to save memory


def tierBoxPlot(sampDF, targetD, dayL, logy=False, label='prevalence',
                colKey='prev_sample', figlbl=None, longlbl=None):

    tierTargetD = defaultdict(set)
    if targetD is not None:
        for (fac, tier), val in targetD.items():  # @UnusedVariable
            tierTargetD[tier].add(val)

    labelL = []
    sampL = []
    markerXL = []
    markerYL = []
    for col, tier in enumerate(sampDF['tier'].unique()):
        labelL.append(tier)
        tierDF = sampDF[sampDF['tier']==tier]
        sampL.append(tierDF[colKey].dropna())
        if targetD is not None:
            for val in tierTargetD[tier]:
                markerXL.append(1.0 + col)
                markerYL.append(val)

    fig, axes = plt.subplots(1, 1)  # @UnusedVariable
    if VERBOSE:
        print 'plotting tier collective boxplots'
    axes.boxplot(sampL, labels=labelL)
    axes.set_title(buildLbl('All Tiers', figlbl, longlbl))
    if markerYL:
        axes.plot(markerXL, markerYL, 'D')
    if logy:
        axes.set_yscale('log')
    plt.tight_layout()
    plt.savefig('%s_pooled_over_tiers.png' % label)


def deriveCols(df, dayL):
    df['prev_sample'] = (df['COLONIZED_sum'].astype(float)
                             / df['TOTAL_sum'].astype(float))
    df['newColonized_sample'] = df['newColonized_sum'].astype(float) / len(dayL)
    df['TOTAL_sample'] = df['TOTAL_sum'].astype(float)/len(dayL)
    df['arrive_colonized_sample'] = (df['creArrivals_sum'].astype(float)
                                         / df['arrivals_sum'].astype(float))
    df['incidence_sample'] = (df['newColonized_sum'].astype(float)
                              / df['localtierdepartures_sum'].astype(float))
    df['arrivals_sample'] = (df['arrivals_sum'].astype(float) / len(dayL))
    return df


def main(argv=None):
    '''Command line options.'''

    global VERBOSE
    global DEBUG

    program_name = os.path.basename(sys.argv[0])
    program_version = "v0.1"
    program_build_date = "%s" % __updated__

    program_version_string = '%%prog %s (%s)' % (program_version, program_build_date)
    #program_usage = '''usage: spam two eggs''' # optional - will be autogenerated by optparse
    program_longdesc = '''''' # optional - give further explanation about what the program does
    program_license = "Copyright 2018 Joel Welling (PSC)                                            \
                Licensed under the GNU Affero General Public License\nhttp://www.gnu.org/licenses/agpl-3.0"

    if argv is None:
        argv = sys.argv[1:]
    try:
        # setup option parser
        parser = OptionParser(version=program_version_string, epilog=program_longdesc,
                              description=program_license)
        parser.add_option("-s", "--sampfile", dest="sampfile", action="append",
                          help="msgpack file containing a pandas dataframe of samples",
                          metavar="FILE")
        parser.add_option("--glob", dest="glob", action="store_true",
                          help="apply filename globbing to sampfile values", metavar="FLAG")
        parser.add_option("-t", "--target", dest="target", action="store",
                          help="pkl file containing target value info [default: %default]",
                          metavar="FILE")
        parser.add_option("-y", "--tauopts", dest="tauopts", action="store",
                          help=("optional yaml file containing taumod options [default: %default]"
                                " If omitted, the last 90 days of data will be used"),
                          metavar="FILE")
        parser.add_option("-v", "--verbose", dest="verbose", action="count",
                          help="set verbosity level [default: %default]")
        parser.add_option("--log", dest="log", action="store_true",
                          help="use log scale on the Y axis", metavar="FLAG")
        parser.add_option("--show", dest="show", action="store",
                          help=("What to plot.  One of {0} [default: %default]"
                                .format(', '.join(SHOW_OPTIONS_DICT.keys()))))
        parser.add_option("--label", dest="figlbl", action="store",
                          help="A label for the figures [default: %default]")

        # set defaults
        parser.set_defaults(tauopts=None, target=None, log=False,
                            show='prevalence', figlbl=None)

        # process options
        (opts, args) = parser.parse_args(argv)

        if args:
            parser.error('This program takes no arguments, just options')

        if opts.verbose > 0:
            VERBOSE = True
            if opts.verbose > 1:
                DEBUG = True

        if not opts.sampfile:
            parser.error('At least one sample file is required')

        if opts.show not in SHOW_OPTIONS_DICT:
            parser.error('Invalid option to --show')

        # MAIN BODY #
        if opts.tauopts is None:
            dayL = range(90)
        else:
            with open(opts.tauopts, 'rU') as f:
                tauOpts = yaml.load(f)
            dayL = tauOpts['DayList']

        if opts.target is None:
            targetD = None
        else:
            with open(opts.target, 'rU') as f:
                targetD = pickle.load(f)

        sampFileL = opts.sampfile
        if opts.glob:
            sampFileL = expandGlobbedList(sampFileL)
        sampDFL = []
        maxDay = None
        for idx, sampFile in enumerate(sampFileL):
            print 'parsing ', sampFile
            df = pd.read_msgpack(sampFile)
            df['run'] = idx
            mD = df['day'].max()
            if maxDay is None:
                maxDay = mD
            else:
                if mD != maxDay:
                    raise RuntimeError('File %s max day %s does not match %s from earlier sample files'
                                       % (sampFile, mD, maxDay))
            days = [maxDay - day for day in dayL]
            sampDFL.append(df[df.day.isin(days)])
        sampDF = pd.concat(sampDFL)
        sampDF['TOTAL'] = sampDF[['COLONIZED', 'CLEAR', 'INFECTED', 'CHRONIC',
                                  'RECOVERED', 'UNDETCOLONIZED']].sum(axis=1)
        sampDF = sampDF.groupby(['tier', 'abbrev', 'run']).sum()  # Sum over days and wards
        sampDF = sampDF.add_suffix('_sum').reset_index()
        poolDF = deriveCols(sampDF.groupby(['tier', 'run']).sum().reset_index(), dayL)
        sampDF = deriveCols(sampDF, dayL)


        colKey, label, lnglbl = SHOW_OPTIONS_DICT[opts.show]
        for tier in sampDF['tier'].unique():
            if tier != 'HOME':
                facBoxPlots(sampDF, targetD, tier, logy=opts.log,
                            label=label, figlbl=opts.figlbl, longlbl=lnglbl,
                            colKey=colKey)
        tierBoxPlot(poolDF, targetD, dayL, logy=opts.log,
                    label=label, figlbl=opts.figlbl, longlbl=lnglbl,
                    colKey=colKey)

    except Exception, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help\n")
        raise
        return 2


if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'prevalence_boxplots_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
