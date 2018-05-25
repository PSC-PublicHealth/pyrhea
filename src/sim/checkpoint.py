#! /usr/bin/env python

###################################################################################
# Copyright   2017, Pittsburgh Supercomputing Center (PSC).  All Rights Reserved. #
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

from __future__ import print_function
import random
import sys

class checkpoint(object):
    def __init__(self, ticks):
        import dmtcp
        self.ticks = ticks
        self.dmtcp = dmtcp

        if dmtcp.isEnabled is False:
            raise RuntimeError("dmtcp is not enabled.  Checkpointing can't work without it")
        print("checkpointing turned on")

    def checkpoint(self, ticks):
        if ticks != self.ticks:
            return
        
        dmtcp = self.dmtcp
        
        dmtcp.checkpoint()

        if dmtcp.isResume():
            print("The checkpoint filename is: %s"%dmtcp.checkpointFilename())
            sys.exit()

        # at this point we should be a restarted image.  Let's invoke a new random seed.
        random.seed()

        # is there anything that really needs done?
    
