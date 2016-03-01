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

import sys
import phacsl.utils.formats.yaml_tools as yaml_tools
import math

import numpy as np
import matplotlib.pyplot as plt

from scipy.cluster.vq import whiten, kmeans, vq
from scipy.stats import lognorm, expon

# Truncate the sample set at how many days during fitting?
truncLim = 300


# def plotCurve(axes, parmVec, scale, rng, nbins, pattern='r-'):
#     curveX = np.linspace(rng[0], rng[1], nbins)
#     curveY = (fullPDF(parmVec, curveX) * scale * ((rng[1]-rng[0])/nbins))
#     axes.plot(curveX, curveY, pattern, lw=2, alpha=0.6)


def main():

    gp1Keys, gp1Recs = yaml_tools.parse_all('/home/welling/workspace/pyRHEA/models/OrangeCounty/facilityfactsCurrent')
    losDict1 = {r['abbrev']: r['losModel'] for r in gp1Recs
                if 'losModel' in r and r['losModel']['pdf'] != 'undefined'
                and r['losModel']['pdf'].startswith('$0*lognorm')}
    gp2Keys, gp2Recs = yaml_tools.parse_all('/home/welling/workspace/pyRHEA/models/OrangeCounty2013/facilityfacts3')
    losDict2 = {r['abbrev']: r['losModel'] for r in gp2Recs
                if 'losModel' in r and r['losModel']['pdf'] != 'undefined'
                and r['losModel']['pdf'].startswith('$0*lognorm')}
    goodKeySet = set(losDict1.keys()).intersection(set(losDict2.keys()))
    goodKeyList = list(goodKeySet)
    goodKeyList.sort()

    clrs = ['red', 'blue', 'green', 'yellow']

    xIndex = 1
    yIndex = 2
    labels = ['k', 'mu', 'sigma', 'lmda']

#     xVals = ([losDict1[k]['parms'][xIndex] for k in goodKeyList]
#              + [losDict2[k]['parms'][xIndex] for k in goodKeyList])
#     yVals = ([losDict1[k]['parms'][yIndex] for k in goodKeyList]
#              + [losDict2[k]['parms'][yIndex] for k in goodKeyList])
#     xMin = min(xVals)
#     xMax = max(xVals)
#     yMin = min(yVals)
#     yMax = max(yVals)
    scatterAx = plt.subplot(121)
    polarAx = plt.subplot(122, polar=True)
#     scatterAx.set_xlim([xMin, xMax])
#     scatterAx.set_ylim([yMin, yMax])
    scatterAx.set_xlabel(labels[xIndex])
    scatterAx.set_ylabel(labels[yIndex])
#     for k in goodKeyList:
#         scatterAx.set_title(k)
#         for idx, src in enumerate([losDict1, losDict2]):
#             parms = src[k]['parms']
#             ll = src[k]['negLogLikPerSample']
#             circ = plt.Circle((parms[xIndex], parms[yIndex]), 0.01 * math.sqrt(ll),
#                               color=clrs[idx])
#             fig.gca().add_artist(circ)
    tplVec = []
    for k in goodKeyList:
        print k
        print losDict1[k]
        print losDict2[k]
        tplVec.append((losDict1[k]['parms'][xIndex],
                       losDict1[k]['parms'][yIndex],
                       losDict2[k]['parms'][xIndex],
                       losDict2[k]['parms'][yIndex],
                       k))
    scatterAx.scatter([t[0] for t in tplVec], [t[1] for t in tplVec],
                      c=clrs[0], marker='o', s=100)
    scatterAx.scatter([t[2] for t in tplVec], [t[3] for t in tplVec],
                      c=clrs[1], marker='o', s=100)
    for tpl in tplVec:
        x0, y0, x1, y1, abbrev = tpl
        scatterAx.arrow(x0, y0, x1-x0, y1-y0, length_includes_head=True)
        scatterAx.annotate(abbrev, xy=(x0, y0), xytext=(x0 + 0.0125, y0 + 0.0125))
    scatterAx.grid(True)

    polarVec = []
    for x0, y0, x1, y1, abbrev in tplVec:
        len = math.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0))
        theta = math.atan2((y1 - y0), (x1 - x0))
        polarVec.append((len, theta, abbrev))
    polarAx.scatter([t[1] for t in polarVec], [t[0] for t in polarVec],
                    marker='o', s=100)
    for r, theta, abbrev in polarVec:
        polarAx.annotate(abbrev, xy=(theta, r), xytext=(theta + 0.05, r))


    plt.show()


############
# Main hook
############

if __name__ == "__main__":
    main()
