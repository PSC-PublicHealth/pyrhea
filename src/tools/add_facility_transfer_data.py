#! /usr/bin/env python

"""
This recreates a set of facilityfacts YAML files with transfer count data added.  The
source of the data is a csv file in a particular format.
"""

import os.path
import yaml
import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools

modelDir = '/home/welling/git/pyrhea/models/ChicagoLand'
transferCSVFileName = 'Matrices_LOS_08092016_Transfer3day.csv'
losCSVFileName = 'HistogramTable_ActualLOS_PROTECT_082516.csv'

allKeySet, recs = yaml_tools.parse_all(os.path.join(modelDir, 'facilityfacts'))
facDict = {r['abbrev']:r for r in recs}

with open(os.path.join(modelDir, transferCSVFileName)) as f:
    totTransKeys, totTransRecs = csv_tools.parseCSV(f)

with open(os.path.join(modelDir, losCSVFileName)) as f:
    losHistoKeys, losHistoRecs = csv_tools.parseCSV(f)

cumTotDict = {}
for rec in losHistoRecs:
    abbrev = rec['UNIQUE_ID']
    if abbrev in cumTotDict:
        cumTotDict[abbrev] += int(rec['cumulative_total'])
    else:
        cumTotDict[abbrev] = int(rec['cumulative_total'])

transferOutDictDict = {}
for rec in totTransRecs:
    fmAbbrev = rec['UNIQUE_ID']
    assert fmAbbrev not in transferOutDictDict, 'two entries for %s' % fmAbbrev
    innerD = {}
    totTransfersOut = 0
    for key in rec:
        if key.startswith('To_'):
            toAbbrev = key[3:]
            innerD[toAbbrev] = int(rec[key])
            totTransfersOut += int(rec[key])
    transferOutDictDict[fmAbbrev] = innerD

transferInDict = {}
for fmAbbrev, toDict in transferOutDictDict.items():
    for toAbbrev, toCt in toDict.items():
        if toAbbrev in transferInDict:
            transferInDict[toAbbrev] += toCt
        else:
            transferInDict[toAbbrev] = toCt
            
for rec in totTransRecs:
    abbrev = rec['UNIQUE_ID']
    totDischarges = (cumTotDict[abbrev] if abbrev in cumTotDict else None)
    totTransfersIn = (transferInDict[abbrev] if abbrev in transferInDict else None)
    totTransfersOut = sum([v for v in transferOutDictDict[abbrev].values()])
    print ('%s: %d total transfers out, %s total transfers in, %s LOS entries' %
           (abbrev, totTransfersOut, totTransfersIn, totDischarges))

dischProvStr = 'HistogramTable_ActualLOS_PROTECT_082516.xlsx ra8b2c0a7b count by UNIQUE_ID'
transProvStr = 'Matrices_LOS_08092016.xlsx:Transfer3day r0134f4795 grouped by category'
newRecs = []
for abbrev, rec in facDict.items():
    nR = rec.copy()
    if abbrev in cumTotDict:
        nR['totalDischarges'] = {'value': cumTotDict[abbrev],
                                 'prov': dischProvStr}
    if abbrev in transferInDict:
        nR['totalTransfersIn'] = {'value': transferInDict[abbrev],
                                  'prov': transProvStr}
    if abbrev in transferOutDictDict:
        byCatDict = {}
        for toAbbrev, ct in transferOutDictDict[abbrev].items():
            if ct != 0:
                if toAbbrev not in facDict:
                    print ('Transfer of %d patients from %s to excluded location %s' %
                           (ct, abbrev, toAbbrev))
                else:
                    category = facDict[toAbbrev]['category']
                    if category in byCatDict:
                        byCatDict[category] += ct
                    else:
                        byCatDict[category] = ct
        nR['totalTransfersOut'] = []
        for category, ct in byCatDict.items():
            nR['totalTransfersOut'].append({'category': category,
                                            'count': {'value': ct,
                                                      'prov': transProvStr}})
    newRecs.append(nR)

print '%d of %d records modified' % (len(newRecs), len(recs))
yaml_tools.save_all(os.path.join(modelDir, 'facilityfactsUpdated'), newRecs)

print 'writing the full transfer matrix in yaml format:'
cleanDictDict = {}
for fmAbbrev, toDict in transferOutDictDict.items():
    cleanToDict = {toAbbrev: ct for toAbbrev, ct in toDict.items() if ct != 0}
    cleanDictDict[fmAbbrev] = cleanToDict
with open('transfer_counts.yaml', 'w') as f:
    yaml.safe_dump(cleanDictDict, f,
                   default_flow_style=False, indent=4, encoding='utf-8',
                   width=130, explicit_start=True)

    
    

