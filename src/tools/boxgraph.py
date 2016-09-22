#! /usr/bin/env python
from difflib import Match

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

"""
Estimate LOS model distribution parameters by optimizing the Kolomogorov-Smirnov statistic
"""

import math
import sys
import traceback
import os.path
import phacsl.utils.formats.csv_tools as csv_tools
from map_transfer_matrix import parseFacilityData

import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import matplotlib.path as path
import matplotlib.patches as patches


def fixNoneEntries(recL):
    """If the lnLikPerSample comes in as the string 'None', replace with None"""
    newL = []
    for rec in recL:
        if rec['lnLikPerSample'] == 'None':
            rec['lnLikPerSample'] = None
        newL. append(rec)
    return newL


def main():

    modelDir = '/home/welling/git/pyrhea/models/ChicagoLand'
    with open(os.path.join(modelDir, 'los_model_fit_parms_lognorm.csv')) as f:
        keys, lnRecs = csv_tools.parseCSV(f)
    lnRecs = fixNoneEntries(lnRecs)
    l = []
    lnDict = {r['abbrev']: r for r in lnRecs}
#     with open(os.path.join(modelDir, 'los_model_fit_parms_twopop.csv')) as f:
    with open(os.path.join(modelDir, 'los_model_fit_parms_bimodal.csv')) as f:
        keys, tpRecs = csv_tools.parseCSV(f)
    tpRecs = fixNoneEntries(tpRecs)
    tpDict = {r['abbrev']: r for r in tpRecs}
    facDict = parseFacilityData(os.path.join(modelDir, 'facilityfacts'))
    assert set(lnDict.keys()) == set(tpDict.keys()), 'abbrevs do not Match'
    
    ratioListDict = {}
    for abbrev in lnDict.keys():
        if abbrev not in facDict:
            continue
        assert abbrev in tpDict, ('%s has a model fit for lognormal but not for two-pop?' %
                                  abbrev)
        category = facDict[abbrev]['category']
        if category not in ratioListDict:
            ratioListDict[category] = []
        if lnDict[abbrev]['lnLikPerSample'] is None:
            if tpDict[abbrev]['lnLikPerSample'] is None:
                continue
            else:
                print ('%s: %s vs %s: both should be None' %
                       (abbrev, lnDict[abbrev]['lnLikPerSample'],
                        tpDict[abbrev]['lnLikPerSample']))
        else:
            if tpDict[abbrev]['lnLikPerSample'] is None:
                print ('%s: %s vs %s: both should be None' %
                       (abbrev, lnDict[abbrev]['lnLikPerSample'],
                        tpDict[abbrev]['lnLikPerSample']))
            else:
                ratio = lnDict[abbrev]['lnLikPerSample']/tpDict[abbrev]['lnLikPerSample']
                if ratio < 0.999:
                    print '!!!! ',
                if ratio < 10.0:
                    print '%s: %s vs %s -> %s' % (abbrev, lnDict[abbrev]['lnLikPerSample'],
                                                  tpDict[abbrev]['lnLikPerSample'], ratio)
        ratioListDict[category].append(ratio)
          
    fig, axes = plt.subplots()
    labels = []
    data = []
    for category, ratioList in ratioListDict.items():
        labels.append('%s (%d)' % (category, len(ratioList)))
        data.append(ratioList)
    axes.boxplot(data, labels=labels, showmeans=True)
#     axes.set_yscale('log')
    axes.grid(True)
    axes.set_title('Ln Likelihood Ratio (higher means 2pop is better)', fontsize=10)
    plt.show()

############
# Main hook
############

if __name__ == "__main__":
    main()
