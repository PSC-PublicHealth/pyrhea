#! /usr/bin/env python

"""
This recreates a set of facilityfacts YAML files with LOS models added.  The
sources of the data are csv files; different sources are used by facility type.
"""

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

snfICUProvStr = 'SNFs do not have ICUs'
ltachICUProvStr = 'LTACHs do not have ICUs'

hospFracICUProvStr = '%s 415091df5d Avg_Daily_Census_ICU/Avg_Daily_Census' % bedLosFName

newRecs = []
for abbrev, rec in facDict.items():
    if rec['category'] in ['SNF', 'VSNF', 'NURSINGHOME']:
        nR = rec.copy()
        nR['fracAdultPatientDaysICU'] = {'value': 0.0,
                                         'prov': snfICUProvStr}
        nR['meanLOSICU'] = {'value': 0.0, 'prov': snfICUProvStr}
        newRecs.append(nR)
    elif rec['category'] in ['LTAC', 'LTACH']:
        nR = rec.copy()
        nR['fracAdultPatientDaysICU'] = {'value': 0.0,
                                         'prov': ltachICUProvStr}
        nR['meanLOSICU'] = {'value': 0.0, 'prov': ltachICUProvStr}
        newRecs.append(nR)
    elif rec['category'] in ['HOSPITAL']:
        if abbrev in bedLosInfo:
            v1 = bedLosInfo[abbrev]['Avg_Daily_Census_ICU']
            v2 = bedLosInfo[abbrev]['Avg_Daily_Census']
            if typeCheck(v1) and typeCheck(v2):
                nR = rec.copy()
                nR['fracAdultPatientDaysICU'] = {'value': float(v1)/float(v2),
                                                 'prov': hospFracICUProvStr}
                newRecs.append(nR)
                
print '%d of %d records modified' % (len(newRecs), len(recs))
yaml_tools.save_all(os.path.join(modelDir, 'facilityfactsUpdated'), newRecs)

    
    

