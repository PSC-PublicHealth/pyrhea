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

_rhea_svn_id_ = "$Id$"

###############################
# A good command line for the resulting .dot file is:
#
#   neato -Gspline=true -Nfontsize=5 -Nheight=0.1 -Efontsize=5 -n2 -Tsvg -ograph.svg graph.dot
#
###############################

import os.path

from map_transfer_matrix import parseFacilityData, writeDotGraph


def main():
    #facDict = parseFacilityData('/home/welling/workspace/pyRHEA/models/OrangeCounty/'
    #                            'facilityfactsCurrent')
    modelDir = '/home/welling/workspace/pyRHEA/models/OrangeCounty2013/'
    facDict = parseFacilityData([
                                 os.path.join(modelDir, 'facilityfactsCurrent2013'),
                                 os.path.join(modelDir, 'synthCommunities')
                                 ])

    title = 'All Facilities'
    transferDict = {}

    transInDict, transOutDict = writeDotGraph('graph.dot', title,  # @UnusedVariables
                                              facDict, transferDict)


if __name__ == "__main__":
    main()
