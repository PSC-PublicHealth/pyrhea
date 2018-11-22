#!/usr/bin/env python
# encoding: utf-8
'''
plot_fuzzy_time_series -- plot a set of statistical time series curves

plot_fuzzy_time_series is a command line program to plot multiple statistical time series.

@author:     welling

@copyright:  2018 Pittsubrgh Supercomputing Center (PSC). All rights reserved.

@license:    AGPL-3.0

@contact:    welling@psc.edu
@deffield    updated: Updated
'''

import sys
import os

from optparse import OptionParser
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Hack for no display w/ pyplot. Change pyplot backend *before* loading it.
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple

from prevalence_boxplots import expandGlobbedList
from plotting_utilities import pltCurvesWithBounds
from pathogenbase import PthStatus

__all__ = []
__version__ = 0.1
__date__ = '2018-11-21'
__updated__ = '2018-11-21'

VERBOSE = False
DEBUG = False
TESTRUN = 0
PROFILE = 0

def parseSampleFiles(sampFileL):
    """
    Given a list of msgpack files, parse them and return the collected dataset.
    
    'run' and 'TOTAL' columns are added.
    """
    sampDFL = []
    sampDFLWithRun = []
    maxDay = None
    minDay = None
    for sampFile in sampFileL:
        if VERBOSE:
            print 'parsing ', sampFile
        df = pd.read_msgpack(sampFile)
        if 'run' in df:
            sampDFLWithRun.append(df)
        else:
            sampDFL.append(df)
        mxD = df['day'].max()
        mnD = df['day'].min()
        if maxDay is None:
            maxDay = mxD
        else:
            if mxD != maxDay:
                raise RuntimeError('File %s max day %s does not match %s from earlier sample files'
                                   % (sampFile, mxD, maxDay))
        if minDay is None:
            minDay = mnD
        else:
            if mnD != minDay:
                raise RuntimeError('File %s min day %s does not match %s from earlier sample files'
                                   % (sampFile, mnD, minDay))

        #sampDFL.append(df[df.day.isin(days)])
    runBase = 0
    for sampDF in sampDFLWithRun:
        sampDF['run'] += runBase
        runBase = sampDF['run'].max() + 1
    for sampDF in sampDFL:
        sampDF['run'] = runBase
        runBase += 1
    sampDF = None
    if VERBOSE:
        print 'Merging datasets...'
    for idx, df in enumerate(sampDFLWithRun + sampDFL):
        if 'ward' in df:
            # Sum over wards within a sample
            df = df.groupby(['tier', 'abbrev', 'run', 'day']).sum().reset_index()
            df = df.drop(columns=['ward'])
        sampDF = pd.concat([sampDF, df])
        if DEBUG:
            print 'merging dataset %s' % idx
    if VERBOSE:
        print 'Merge compete'
    pthNameL = PthStatus.names.values()[:]
    sampDF['TOTAL'] = sampDF[pthNameL].sum(axis=1)
    if DEBUG:
        print 'Added TOTAL'

    return sampDF

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
    program_license = "Copyright 2018 Joel Welling, Pittsburgh Supercomputing Center, CMU\
                Licensed under the GNU AGPL 3.0\nhttp://www.gnu.org/licenses/agpl-3.0.html"

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
        parser.add_option("-v", "--verbose", dest="verbose", action="count",
                          help="set verbosity level [default: %default]")
        parser.add_option("--log", dest="log", action="store_true",
                          help="use log scale on the Y axis", metavar="FLAG")
        parser.add_option("--show", dest="show", action="store",
                          help="What to plot.  One of 'prevalence', 'incidence' [default: %default]")
        parser.add_option("--label", dest="figlbl", action="store",
                          help="A label for the figure [default: %default]")
        parser.add_option("-o", "--out", dest="out", action="store",
                          help="Name for output image file [default: %default]")

        # set defaults
        parser.set_defaults(tauopts="./taumod_config.yaml", target="expected.pkl", log=False,
                            show='prevalence', figlbl=None, out="timeseries.png")

        # process options
        (opts, args) = parser.parse_args(argv)

        if opts.verbose > 0:
            VERBOSE = True
            if opts.verbose > 1:
                DEBUG = True

        if not opts.sampfile:
            parser.error('At least one sample file is required')
        if opts.show not in ['prevalence', 'incidence']:
            parser.error('Invalid option to --show')

        # MAIN BODY #

        sampFileL = opts.sampfile
        if opts.glob:
            sampFileL = expandGlobbedList(sampFileL)
        sampDF = parseSampleFiles(sampFileL)

        if VERBOSE:
            print 'Merging individual facilities'
            sampDF = sampDF.groupby(['day', 'tier', 'run']).sum().reset_index()

        if DEBUG:
            print 'Calculating prevalence'
        sampDF['prevalence'] = sampDF['COLONIZED'].astype(float)/sampDF['TOTAL'].astype(float)

        if DEBUG:
            print 'Dropping irrelevant columns'
        for col in sampDF.columns:
            if col not in ['day', 'tier', 'prevalence', 'run', 'newColonized']:
                sampDF = sampDF.drop(columns=[col])

        fig, axes = plt.subplots(1, 1)
        allArtists = []
        allLabels = []
        key, lbl = {'prevalence': ('prevalence', 'prevalence'),
                    'incidence': ('newColonized', 'incidence')}[opts.show]
        if VERBOSE:
            print 'Starting pltCurvesWithBounds'
        tierL = sampDF['tier'].unique()[:]
        tierL.sort()
        for tier in tierL:
            if DEBUG:
                print 'plotting %s' % tier
            artistL, labelL = pltCurvesWithBounds(sampDF, axes, key, 'day', 'tier', [tier],
                                   ['baseline %s' % tier])
            allArtists += artistL
            allLabels += labelL
        if DEBUG:
            print 'Finished pltCurvesWithBounds'
        axes.legend(allArtists, allLabels, handler_map={tuple: HandlerTuple()})
        titleStr = '{0} {1} ({2} realizations)'.format('' if opts.figlbl is None else opts.figlbl,
                                                    lbl, sampDF['run'].nunique())
        axes.set_title(titleStr)
        if opts.log:
            axes.set_yscale('log')
        fig.set_size_inches(8.0, 6.0)
        plt.savefig(opts.out)

    except Exception, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
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
        profile_filename = 'plot_fuzzy_time_series_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())