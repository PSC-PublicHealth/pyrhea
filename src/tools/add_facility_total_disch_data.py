#! /usr/bin/env python

"""
This recreates a set of facilityfacts YAML files with totalDischarges added.  The
sources of the data are csv files.
"""

import os.path
import yaml
import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools

def loadCSVByAbbrev(modelDir, fname, key='abbrev'):
    fullName = os.path.join(modelDir, fname)
    with open(fullName) as fl:
        keys, recs = csv_tools.parseCSV(fl)
    assert key in keys, ('%s has no "%s" field' % (fullName, key))
    return {rec[key]: rec for rec in recs}

modelDir = '/home/welling/git/pyrhea/models/ChicagoLand'

allKeySet, recs = yaml_tools.parse_all(os.path.join(modelDir, 'facilityfacts'))
facDict = {r['abbrev']:r for r in recs}

infoFName = 'Matrices_LOS_09292016_Facilities_LOS.csv'
infoDict = loadCSVByAbbrev(modelDir, infoFName, key='UNIQUE_ID')

totDischProvStr = """
Matrices_LOS_08092016_Facilities_LOS ccda79e0 (which matches
Matrices_LOS_08092016:Facilities_LOS), column 'N Obs'
"""

newRecs = []
for abbrev, rec in facDict.items():
    nR = rec.copy()
    if 'totalDischarges' in nR:
        del nR['totalDischarges']
    if abbrev in infoDict:
        nR['totalDischarges'] = {'value': infoDict[abbrev]['N Obs'],
                                 'prov': totDischProvStr}
        newRecs.append(nR)
        totTransfers = sum([d['count']['value'] for d in rec['totalTransfersOut']])
        print '%s : %s were transferred' % (abbrev,
                                            float(totTransfers)/float(infoDict[abbrev]['N Obs']))

print '%d of %d records modified' % (len(newRecs), len(recs))
yaml_tools.save_all(os.path.join(modelDir, 'facilityfactsUpdated'), newRecs)

    
    

