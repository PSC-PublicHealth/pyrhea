#! /usr/bin/env python

"""
Generate 'community' records based on no data
"""
_rhea_svn_id_ = "$Id$"

import yaml_tools


outDir = '/home/welling/workspace/pyRHEA/models/OrangeCounty/communitiesCurrent'
totPop = 3000000
nCommunities = 100

blockSz = totPop / nCommunities
newRecs = []
for i in xrange(nCommunities):
    abbrev = 'CM%02d' % i
    if i == nCommunities - 1:
        pop = totPop
    else:
        pop = blockSz
    totPop -= blockSz
    rec = {'abbrev': abbrev,
           'name': 'Synthetic community %d' % i,
           'category': 'COMMUNITY',
           'meanPop': pop,
           'meanPop_prov': 'completely synthetic from %s' % _rhea_svn_id_
           }
    newRecs.append(rec)

print '%d records created' % len(newRecs)
yaml_tools.save_all(outDir, newRecs)
