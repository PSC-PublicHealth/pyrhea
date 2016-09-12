#! /usr/bin/env python

"""
This assumes that address data has been manually added to the yaml files, and uses Google's geocoding
API to add latitude and longitude.
"""

import json
import yaml_tools
import csv_tools
import sys

allKeySet, recs = yaml_tools.parse_all('facilityfacts9')
facDict = {r['abbrev']:r for r in recs}

with open('/home/welling/workspace/pyRHEA/models/OrangeCounty/gatherdata5.csv') as f:
    totDischKeys, totDischRecs = csv_tools.parseCSV(f)

with open('/home/welling/workspace/pyRHEA/models/OrangeCounty/gatherdata6.csv') as f:
    totTransKeys, totTransRecs = csv_tools.parseCSV(f)

nameDict = {r['name'].upper(): r['abbrev'] for r in recs}

transferInDict = {}
for r in totTransRecs:
    for k,v in r.items():
        if k.startswith('To_'):
            assert k[3:] in facDict
            abbrev = k[3:]
            if abbrev not in transferInDict:
                transferInDict[abbrev] = 0
            transferInDict[abbrev] += v
            
transferOutDict = {}
for r in totTransRecs:
    fmAbbrev = r['']
    if fmAbbrev not in transferOutDict:
        transferOutDict[fmAbbrev] = {}
    for k,v in r.items():
        if k.startswith('To_'):
            toAbbrev = k[3:]
            toCat = facDict[toAbbrev]['category']
            if toCat not in transferOutDict[fmAbbrev]:
                transferOutDict[fmAbbrev][toCat] = 0
            transferOutDict[fmAbbrev][toCat] += v
        

dischargeDict = {}
for r in totDischRecs:
    if r['name'] in nameDict:
        dischargeDict[nameDict[r['name']]] = r['count']
    else:
        print '%s: no corresponding facility' % r['name']

dischProvStr = '2007_Annual_Admissions__Bed_Sizes_TA DKsmb.xls:AnnualAdmissions $B'
transProvStr = '2007_OC_Nursing_Homes_and_Hospitals_Patient_Flow_UPS_RHEA.xlsx:Full_Grid_direct_SelfR:$B2:$CX102'
newRecs = []
for rec in recs:
    nR = rec.copy()
    abbrev = rec['abbrev']
    nR['totalDischarges'] = {'value': dischargeDict[abbrev],
                             'prov': dischProvStr}
    if abbrev in transferInDict:
        nR['totalTransfersIn'] = {'value': transferInDict[abbrev],
                                  'prov': transProvStr}
    if abbrev in transferOutDict:
        nR['totalTransfersOut'] = []
        toTot = 0
        for k, v in transferOutDict[abbrev].items():
            nR['totalTransfersOut'].append({'category': k,
                                            'count': {'value': v,
                                                      'prov': transProvStr}})
            toTot += v
        print '%s: discharges %s, transfers in: %s out: %s' % (abbrev, dischargeDict[abbrev],
                                                               transferInDict[abbrev],
                                                               toTot)
    newRecs.append(nR)

print '%d records modified' % len(newRecs)
yaml_tools.save_all('facilityfacts10', newRecs)
