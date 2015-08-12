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

_hermes_svn_id_="$Id: full_import_paths.py 2262 2015-02-09 14:38:25Z stbrown $"

import sys, os, os.path

cwd = os.path.dirname(__file__)

if cwd not in sys.path:
    try:
        sys.path.append(cwd)
    except:
        pass

paths = ['../../src/',
         '../../src/sim/',
         '../../src/tools/',
         ]

try:
    fullpaths = [os.path.join(cwd, path) for path in paths]
    fullpaths = [p for p in fullpaths if p not in sys.path]
    sys.path = fullpaths + sys.path
except:
    pass

