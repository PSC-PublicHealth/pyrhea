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

_rhea_svn_id_="$Id$"

import sys
import os.path
import optparse

import schemautils


def main(myargv=None):
    "Provides a few test routines"

    if myargv is None:
        myargv = sys.argv

    parser = optparse.OptionParser(usage="""
    %prog [-R] schemaFile checkFile

    The schema and the file(s) to be checked can be JSON or YAML.
    """)

    parser.add_option("-R", action="store_true",
                      help="If checkFile is a directory, check files recursively")

    opts, args = parser.parse_args()

    if len(args) == 2:
        validator = schemautils.getValidator(args[0])

        if opts.R and os.path.isdir(args[1]):
            for root, dirs, files in os.walk(args[1]):  # @UnusedVariable
                for fName in files:
                    if fName.lower().endswith(('.yaml', '.yml', '.json', '.jsn')):
                        path = os.path.join(root, fName)
                        try:
                            nErrors = sum([1 for e in
                                           validator.iter_errors(schemautils.fileToJSON(path))])
                            if nErrors:
                                print 'FAIL: %s: %d errors' % (path, nErrors)
                            else:
                                print 'PASS: %s' % path
                        except Exception, e:  # @UnusedVariable
                            print ('FAIL: %s %s' % (path, e))
        else:
            try:
                nErrors = 0
                for error in validator.iter_errors(schemautils.fileToJSON(args[1])):
                    print '%s: %s' % (' '.join([str(word) for word in error.path]), error.message)
                    nErrors += 1
                if nErrors:
                    print 'FAIL: %s: %d errors' % (args[1], nErrors)
                else:
                    print 'PASS: %s' % args[1]
            except Exception, e:
                print 'FAIL: %s: %s' % (args[1], e)
    else:
        parser.print_help()


############
# Main hook
############

if __name__ == "__main__":
    main()
