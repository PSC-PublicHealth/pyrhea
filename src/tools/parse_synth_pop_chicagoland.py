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

import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools
import math
import os

import numpy as np
import matplotlib.pyplot as plt

#countyFIPS = ['06059']  # orange county
countyFIPS = [18073, 18089, 18091, 18111, 18127,
              17031, 17037, 17043, 17063, 17089, 17091, 17093, 17097, 17111, 17197,
              55059]
countyFIPS = [str(fips) for fips in countyFIPS]
countyFIPS = set(countyFIPS)
thisDir = os.path.dirname(__file__)
baseDir = '/home/jim/mnt/olympus_root/mnt/beegfs1/data/shared_group_data/syneco/input/west/north_america/united_states/'
srcDirs = [os.path.join(baseDir, '17'), os.path.join(baseDir, '18'), os.path.join(baseDir, '55')]
outDir = 'synthCommunities'

#srcDir = os.path.join(thisDir, '../../models/OrangeCounty2013/epi_synth_pop_2')
#outDir = os.path.join(thisDir, '../../models/OrangeCounty2013/synthCommunities')

           
def recurseFiles(dirs):
    dirStack = dirs

    while True:
        if not len(dirStack):
            break
        dir = dirStack.pop()

        try:  # deal with permission errors and the like.
            for fn in os.listdir(dir):
                ffn = os.path.join(dir, fn)
                if os.path.isfile(ffn):
                    yield(ffn)
                elif os.path.isdir(ffn):
                    dirStack.append(ffn)
                else:
                    pass
        except:
            pass


           
def importLOSTable(fname):
    with open(fname, 'r') as f:
        keys, losRecs = csv_tools.parseCSV(f)  # @UnusedVariable

    losListDict = {}
    for r in losRecs:
        k = r['# Abbreviation']
        v = r['LOS']
        if k not in losListDict:
            losListDict[k] = []
        losListDict[k].append(v)
#     for k,v in losListDict.items():
#         print '%s: %d' % (k,len(v))
    return losListDict


def fiveNumberSummary(vec):
    return [len(vec),
            np.mean(vec),
            np.std(vec),
            np.percentile(vec, 25.0),
            np.median(vec),
            np.percentile(vec, 75.0),
            ]


def betterSummary(vec):
    fNS = fiveNumberSummary(vec)
    N = fNS[0]
    mean = fNS[1]
    stdv = fNS[2]
    q1 = fNS[3]
    median = fNS[4]
    q3 = fNS[5]
    #return [median/mean, ((q3-q1)/stdv), q1/median, q3/median]
    #return [median, median/mean, (q3-median)/(median-q1)]
    return [median, median/mean, (q3-q1)/median]
    #return [float(N), median/mean, (q3-q1)/median]

def leftovers():
    indexDict = {}
    valVec = []
    offset = 0
    losListDict = importLOSTable('/home/welling/git/rhea-dante/test/nursing_home_CI_decolonization_2014/Length_of_Stay_2007_to_2009_OC_Nursing_Homes-12-18-11_SMB_with_abbrev_RHEA.csv')
    #losListDict = importLOSTable('/home/welling/Dropbox/RHEA Inputs/OC_County_Data_2013/OC_County_Data_2013-clean-FOR-JOEL/'
    #                             + 'length_of_stay_2011_2012_OC_Nursing_Homes_08-12-2014_begin0000-end2012+Fix-SILOS+SCLE.csv')
    #losListDict = importLOSTable('/home/welling/Dropbox/RHEA Inputs/OC_County_Data_2013/OC_County_Data_2013-clean-FOR-JOEL/'
    #                             + 'length_of_stay_2011_2012_OC_Nursing_Homes_08-12-2014_begin0000-end2012+SILOS+SCLE.csv')
    #del losListDict['CMLH']
    tblRecs = []
    allLosSamples = []
    for abbrev, losList in losListDict.items():
        allLosSamples.extend(losList)
        if len(losList) >= 10:
            indexDict[abbrev] = offset
            bSVec = betterSummary(losList)
            valVec.append(bSVec)
            offset += 1
            tblRecs.append({'abbrev': abbrev,
                            'median': bSVec[0],
                            'medianOverMean': bSVec[1],
                            'quartileBalance': bSVec[2]
                            })
    print 'Statistics over all LOS samples:'
    allLossSampleSummary = fiveNumberSummary(allLosSamples)
    for lbl, val in zip(['N', 'mean', 'stdv', 'Q1', 'median', 'Q3'], allLossSampleSummary):
        print '   %s: %s' % (lbl, val)

    reverseMap = {v: k for k, v in indexDict.items()}


def main():
    tractDict = {}
    comCount = 0
    communities = []
           
    for ffn in recurseFiles(srcDirs):
        fn = os.path.basename(ffn)

        if not fn.startswith('people_'):
            continue
        if not fn.endswith('.csv'):
            continue
        parts = os.path.splitext(fn)[0].split('_')
        fileFIPS = parts[1][:5]
        if fileFIPS not in countyFIPS:
            print "%s is from the wrong county, skipped"%fn
            continue

        print "reading %s"%fn
        tract = parts[1][5:]
        if tract in tractDict:
            raise RuntimeError('Redundant tract %s in %s' % (parts[3], fN))

        with open(ffn, 'r') as f:
            keys, recs = csv_tools.parseCSV(f)  # @UnusedVariable
            n = 1
            meanLon = recs[0]['longitude']
            meanLat = recs[0]['latitude']
            for r in recs[1:]:
                n += 1
                meanLon += r['longitude']
                meanLat += r['latitude']
            meanLon /= n
            meanLat /= n
            print 'Tract %s: n = %d, lon = %s, lat = %s' % (tract, n, meanLon, meanLat)
            communities.append({'name': 'SynthPop for FIPS %s tract %s' % (fileFIPS, tract),
                                'abbrev': 'C%03d' % comCount,
                                'category': 'COMMUNITY',
                                'meanPop': {'value': n,
                                            'prov': 'Synthetic population size for this tract'
                                            },
                                'longitude': meanLon,
                                'latitude': meanLat,
                                'lat_lon_prov': ('Means of lat and lon for all agents '
                                                 'in this tract'),
                                'FIPS': fileFIPS,
                                'censusTract': tract
                                })
            comCount += 1

    yaml_tools.save_all(outDir, communities)

        


############
# Main hook
############

if __name__ == "__main__":
    main()
