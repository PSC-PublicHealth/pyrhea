#!/usr/bin/env python
# encoding: utf-8
'''
parse_taumod_out -- parse the output stream produced by taumod.py to produce a pandas dataframe

parse_taumod_out is a utility

It defines parseTaumodOut

@author:     welling

@copyright:  2018 PSC. All rights reserved.

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
from collections import defaultdict
from optparse import OptionParser

import pandas as pd

__all__ = []
__version__ = 0.1
__date__ = '2018-04-30'
__updated__ = '2018-04-30'

TESTRUN = 0
PROFILE = 0
VERBOSITY = 0
DEBUG = 0


def _stripLoggerInfo(someFile):
    recs = someFile.readlines()
    if VERBOSITY:
        print 'read complete; %d lines' % len(recs)
    outRecs = []
    currentRec = None
    dbg = DEBUG
    for rec in recs:
        if 'INFO' in rec:
            dbg = DEBUG
            offset = rec.find('INFO')
            if currentRec is None:
                currentRec = rec[:offset]
            else:
                currentRec = currentRec + rec[:offset]
            if dbg:
                print 'rec: %s' % rec
                print 'currentRec: %s' % currentRec
        else:
            if currentRec is None:
                outRecs.append(rec)
                dbg = False
            else:
                currentRec += rec
                outRecs.append(currentRec)
                currentRec = None
            if dbg:
                print 'rec: %s' % rec
                print 'currentRec: %s' % currentRec
                print 'appended'
    return outRecs


def parseTaumodOut(fname):
    """
    fname points to a file containing the output of taumod.  Returns a Pandas Dataframe.
    """
    with open(fname, 'rU') as f:
        outRecs = _stripLoggerInfo(f)
    parsedRecs = []
    counts = defaultdict(int)
    for rec in outRecs:
        sepL1 = rec.split(',')
        if len(sepL1) == 13:
            sepL10 = sepL1[0].split(':')
            sepL10 = [word.strip() for word in sepL10]
            assert sepL10[1] == 'total', 'bad first sep %s'%sepL1
            fac, tier = sepL10[0].split()
            ent = {'fac': fac, 'tier': tier, 'total': int(sepL10[-1])}
            for phrase in sepL1[1:]:
                words = phrase.split()
                ent[words[0].strip()] = float(words[1].strip())
            ent['samp'] = counts[(ent['fac'], ent['tier'])]
            counts[(ent['fac'], ent['tier'])] += 1
            parsedRecs.append(ent)
    if VERBOSITY:
        print '%d useful records' % len(parsedRecs)
    if DEBUG:
        print parsedRecs[:5]
        print parsedRecs[-5:]
    return pd.DataFrame(parsedRecs)


def main(argv=None):
    '''Command line options.'''

    program_name = os.path.basename(sys.argv[0])
    program_version = "v0.1"
    program_build_date = "%s" % __updated__

    program_version_string = '%%prog %s (%s)' % (program_version, program_build_date)
    #program_usage = '''usage: spam two eggs''' # optional - will be autogenerated by optparse
    program_longdesc = '''''' # optional - give further explanation about what the program does
    program_license = "Copyright 2018 welling (PSC)                                            \
                Licensed under the AGPL-3.0 http://www.gnu.org/licenses/agpl-3.0.html"

    if argv is None:
        argv = sys.argv[1:]
    try:
        # setup option parser
        parser = OptionParser(version=program_version_string, epilog=program_longdesc,
                              description=program_license)
        parser.add_option("-i", "--in", dest="infile", help="set input path [default: %default]",
                          metavar="FILE")
        parser.add_option("-o", "--out", dest="outfile", help="set output path [default: %default]",
                          metavar="FILE")
        parser.add_option("-v", "--verbose", dest="verbose", action="count",
                          help="set verbosity level [default: %default]")

        # set defaults
        parser.set_defaults(outfile="./taumod.mpk", infile="./taumod.out")

        # process options
        (opts, args) = parser.parse_args(argv)
        if args:
            raise RuntimeError('This program takes no arguments, only options')

        if opts.verbose > 0:
            print("verbosity level = %d" % opts.verbose)
            global VERBOSITY, DEBUG
            VERBOSITY = opts.verbose
            DEBUG = (VERBOSITY > 1)
        if opts.infile:
            print("infile = %s" % opts.infile)
        if opts.outfile:
            print("outfile = %s" % opts.outfile)

        df = parseTaumodOut(opts.infile)
        #print df
        df.to_msgpack(opts.outfile)

    except Exception, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
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
        profile_filename = 'parse_taumod_out_profile.cprf'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())