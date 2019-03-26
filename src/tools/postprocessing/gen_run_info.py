#!/usr/bin/env python
# encoding: utf-8
'''
postprocessing.gen_run_info -- This generates some run-related info in a convenient form for bash scripts

postprocessing.gen_run_info is a tool to generate some run-related info in a convenient form for bash scripts

It defines nothing of use elsewhere.

@author:     welling

@copyright:  2019 PSC, CMU

@license:    AGPL-3

@contact:    welling@psc.edu
@deffield    updated: Updated
'''

import sys
import os

from optparse import OptionParser
import os.path
import sys
import yaml
import tools_util as tu
import pyrheautils as pu

__all__ = []
__version__ = 0.1
__date__ = '2019-03-26'
__updated__ = '2019-03-26'

DEBUG = 0
TESTRUN = 0
PROFILE = 0

def main(argv=None):
    '''Command line options.'''

    program_name = os.path.basename(sys.argv[0])
    program_version = "v0.1"
    program_build_date = "%s" % __updated__

    program_version_string = '%%prog %s (%s)' % (program_version, program_build_date)
    #program_usage = '''usage: spam two eggs''' # optional - will be autogenerated by optparse
    program_longdesc = '''''' 
    program_license = "Copyright 2019 welling (PSC, CMU) \
                Licensed under the AGPL 3.0"

    print('#! /usr/bin/bash')

    if argv is None:
        argv = sys.argv[1:]
    try:
        # setup option parser
        parser = OptionParser(version=program_version_string, epilog=program_longdesc,
                              description=program_license)
        parser.add_option("-v", "--verbose", dest="verbose", action="count",
                          help="set verbosity level [default: %default]")

        # process options
        (opts, args) = parser.parse_args(argv)

        if opts.verbose > 0:
            print("#verbosity level = %d" % opts.verbose)

        runYaml = args[0]
        topDir = args[1]
        bundleDir = args[2]
        pyrheaDir = args[3]
        workDir = args[4]
        scenario = args[5]

        fn = os.path.join(bundleDir, runYaml)
        inputDict = tu.readModelInputs(fn)
        pu.prepPathTranslations(inputDict)

        print('burnindays=%s' % inputDict['burnInDays'])
        print('scenariowaitdays=%s' % inputDict['scenarioWaitDays'])
        print('rundurationdays=%s' % inputDict['runDurationDays'])
        print('totalrundays=%s' % (inputDict['burnInDays'] + inputDict['runDurationDays']))


    except Exception, e:
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
        profile_filename = 'postprocessing.gen_run_info_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())