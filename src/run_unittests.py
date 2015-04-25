#!/usr/bin/env python

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

""" run_unittests.py
This module collects and runs unit tests from other modules, based on the python
module 'unittest'.
"""

_rhea_svn_id_ = "$Id$"

import sys
import unittest

moduleNames = ['pyrheautils']


def main():
    modulesToTest = []
    if len(sys.argv) > 1:
        for a in sys.argv[1:]:
            if a in moduleNames:
                modulesToTest.append(__import__(a))
            else:
                print 'Unknown module %s' % a
    else:
        modulesToTest = [__import__(a) for a in moduleNames]

    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for m in modulesToTest:
        suite.addTests(loader.loadTestsFromModule(m))
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)

############
# Main hook
############

if __name__ == "__main__":
    main()
