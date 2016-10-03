#! /usr/bin/env python

"""
This recreates a set of facilityfacts YAML files with population data added or updated.
The sources of the data are csv files; different sources are used by facility type.
"""

import sys
import os.path
import types
import yaml
import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools

def loadCSVByAbbrev(modelDir, fname, key=None):
    if key is None:
        key = 'abbrev'
    fullName = os.path.join(modelDir, fname)
    with open(fullName) as fl:
        keys, recs = csv_tools.parseCSV(fl)
    assert key in keys, ('%s has no "%s" field' % (fullName, key))
    return {rec[key]: rec for rec in recs}

def typeCheck(val):
    return (type(val) in [types.FloatType, types.IntType]) 

modelDir = '/home/welling/git/pyrhea/models/ChicagoLand'

allKeySet, recs = yaml_tools.parse_all(os.path.join(modelDir, 'facilityfacts'))
facDict = {r['abbrev']:r for r in recs}

bedLosFName = 'Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_v082116_BedLOS.csv'
bedLosInfo = loadCSVByAbbrev(modelDir, bedLosFName, key='UNIQUE_ID')

matLosFName = 'Matrices_LOS_09292016_Facilities_LOS.csv'
matLosInfo = loadCSVByAbbrev(modelDir, matLosFName, key='UNIQUE_ID')

blMeanLOSProvStr = '%s rev 415091df Avg_Daily_Census' % bedLosFName
mlMeanLOSProvStr = '%s rev 415091df Mean' % matLosFName

newRecs = []
for abbrev, rec in facDict.items():
    nR = rec.copy()
    changed = False
    if ('meanLOS' in nR and 'value' in nR['meanLOS'] and
            typeCheck(nR['meanLOS']['value'])):
        if (abbrev in bedLosInfo and typeCheck(bedLosInfo[abbrev]['LOS'])):
            assert bedLosInfo[abbrev]['LOS'] == nR['meanLOS']['value'], ('mismatch on %s for %s' %
                                                                         (bedLosFName, abbrev))
        elif (abbrev in matLosInfo and typeCheck(matLosInfo[abbrev]['Mean'])):
            assert matLosInfo[abbrev]['Mean'] == nR['meanLOS']['value'], ('mismatch on %s for %s' %
                                                                          (matLosFName, abbrev))            
    elif (abbrev in bedLosInfo and typeCheck(bedLosInfo[abbrev]['LOS'])):
        nR['meanLOS'] = {'value': bedLosInfo[abbrev]['LOS'],
                         'prov': blMeanLOSProvStr}
        changed = True
    elif (abbrev in matLosInfo and typeCheck(matLosInfo[abbrev]['Mean'])):
        nR['meanLOS'] = {'value': matLosInfo[abbrev]['Mean'],
                         'prov': mlMeanLOSProvStr}
        changed = True
    else:
        print '%s: no coverage'

    if changed:
        newRecs.append(nR)

print '%d of %d records modified' % (len(newRecs), len(recs))
yaml_tools.save_all(os.path.join(modelDir, 'facilityfactsUpdated'), newRecs)

