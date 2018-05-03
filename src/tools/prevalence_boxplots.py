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

It defines prevalenceBoxPlots, tierPrevalenceBoxPlots

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
from matplotlib import rcParams

__all__ = []
__version__ = 0.1
__date__ = '2018-05-02'
__updated__ = '2018-05-02'

DEBUG = False
VERBOSE = False
TESTRUN = 0
PROFILE = 0

BOXES_PER_FIG = 10

def expandGlobbedList(pathList):
    """
    Given a list of path strings, expand with globbing
    """
    newPathList = []
    for fName in pathList:
        #print '%s yields %s' % (fName, glob.glob(fName))
        newPathList += glob.glob(fName)
    return newPathList


def prevalenceBoxPlots(sampDF, targetD, tier):
    tierDF = sampDF[sampDF['tier']==tier]
    tierGps = tierDF.groupby('fac')
    nBlocks = (len(tierGps) / BOXES_PER_FIG) + 1
    tupleL = [tpl for tpl in tierGps]
    fig, axes = plt.subplots(1, 1)  # @UnusedVariable
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
            #print 'fac %s tier %s' % (key, tier)
            #print subDF
            labelL.append(key)
            sampL.append(subDF['prev_sample'].dropna())
            markerXL.append(1.0 + idx)
            markerYL.append(targetD[(key, tier)])
        print 'plotting %s' % str(labelL)
        axes.tick_params(axis='x', labelsize='small', labelrotation=45.0)
        axes.boxplot(sampL, labels=labelL)
        axes.plot(markerXL, markerYL, 'D')
        if nBlocks > 1:
            axes.set_title('%s %d' % (tier, blockNum))
        else:
            axes.set_title(tier)
        axes.set_yscale('log')
        plt.tight_layout()
        plt.savefig('prevalence_%s_%02d.png' % (tier, block))
        plt.cla()  # to save memory


def tierPrevalenceBoxPlot(sampDF, targetD):
    # Need a different summing pattern for this copy of sampDF
    sampDF = sampDF.groupby(['tier', 'day', 'run']).sum()  # Sum over wards within a sample
    sampDF = sampDF.add_suffix('_sum').reset_index()
    sampDF['prev_sample'] = sampDF['COLONIZED_sum'].astype(float)/sampDF['TOTAL_sum'].astype(float)

    tierTargetD = defaultdict(set)
    for (fac, tier), val in targetD.items():  # @UnusedVariable
        tierTargetD[tier].add(val)

    labelL = []
    sampL = []
    markerXL = []
    markerYL = []
    for col, tier in enumerate(sampDF['tier'].unique()):
        labelL.append(tier)
        tierDF = sampDF[sampDF['tier']==tier]
        sampL.append(tierDF['prev_sample'].dropna())
        for val in tierTargetD[tier]:
            markerXL.append(1.0 + col)
            markerYL.append(val)

    fig, axes = plt.subplots(1, 1)  # @UnusedVariable
    print 'plotting tier collective boxplots'
    axes.boxplot(sampL, labels=labelL)
    axes.plot(markerXL, markerYL, 'D')
    axes.set_yscale('log')
    plt.savefig('prevalence_pooled_over_tiers.png')


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
        parser = OptionParser(version=program_version_string, epilog=program_longdesc, description=program_license)
        parser.add_option("-s", "--sampfile", dest="sampfile", action="append",
                          help="msgpack file containing a pandas dataframe of samples", metavar="FILE")
        parser.add_option("--glob", dest="glob", action="store_true",
                          help="apply filename globbing to sampfile values", metavar="FLAG")
        parser.add_option("-t", "--target", dest="target", action="store",
                          help="pkl file containing target prevalence info [default: %default]", metavar="FILE")
        parser.add_option("-y", "--tauopts", dest="tauopts", action="store",
                          help="yaml file containing taumod options [default: %default]", metavar="FILE")
        parser.add_option("-v", "--verbose", dest="verbose", action="count",
                          help="set verbosity level [default: %default]")

        # set defaults
        parser.set_defaults(tauopts="./taumod_config.yaml", target="expected.pkl")

        # process options
        (opts, args) = parser.parse_args(argv)

        if opts.verbose > 0:
            VERBOSE = True
            if opts.verbose > 1:
                DEBUG = True

        if not opts.sampfile:
            parser.error('At least one sample file is required')

        # MAIN BODY #
        with open(opts.tauopts, 'rU') as f:
            tauOpts = yaml.load(f)
        dayL = tauOpts['DayList']

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
        savSampDF = sampDF.copy()
        sampDF = sampDF.groupby(['tier', 'fac', 'day', 'run']).sum()  # Sum over wards within a sample
        sampDF = sampDF.drop(columns=['ward'])
        sampDF = sampDF.add_suffix('_sum').reset_index()
        sampDF['prev_sample'] = sampDF['COLONIZED_sum'].astype(float)/sampDF['TOTAL_sum'].astype(float)
        for tier in sampDF['tier'].unique():
            prevalenceBoxPlots(sampDF, targetD, tier)
        tierPrevalenceBoxPlot(savSampDF, targetD)

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