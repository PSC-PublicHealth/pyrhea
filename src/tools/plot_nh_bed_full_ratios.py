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

import yaml_tools
import numpy as np
import matplotlib.pyplot as plt


def main():
    vertBars = False
    keys, recs = yaml_tools.parse_all_simplified('/home/welling/workspace/pyRHEA/models/OrangeCounty/facilityfactsCurrent')
    recs = [r for r in recs if r['category'] == 'NURSINGHOME']
    ind = np.arange(len(recs))
    pairs = [(float(r['meanPop'])/float(r['nBeds']) if 'meanPop' in r else 0.0,
              r['abbrev']) for r in recs]
    heights = [a for a, b in pairs]  # @UnusedVariable
    labels = [b for a, b in pairs]
    fig, ax = plt.subplots()  # @UnusedVariable
    if vertBars:
        ax.bar(ind, heights, 1.0)
        ax.set_xticklabels(labels)
        ax.set_xticks(ind+0.5)
        ax.set_title("meanPop/nBeds for Nursing Homes")
        ax.plot([0.0, len(heights)+1.0], [1.0, 1.0], 'r-')
        locs, lbls = plt.xticks()  # @UnusedVariable
        plt.setp(lbls, rotation=90)
    else:
        ax.barh(ind, heights, 1.0)
        ax.set_yticklabels(labels)
        ax.set_yticks(ind+0.5)
        ax.set_title("meanPop/nBeds for Nursing Homes")
#         locs, lbls = plt.xticks()  # @UnusedVariable
#        plt.setp(lbls, rotation=90)
        ax.plot([1.0, 1.0], [0.0, len(heights)+1.0], 'r-')
    plt.show()

############
# Main hook
############

if __name__ == "__main__":
    main()
