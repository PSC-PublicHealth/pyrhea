#! /usr/bin/env python

"""
This recreates a set of facilityfacts YAML files with bed count info added.  The
sources of the data are csv files.
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

def hospFracICUOp(rec):
    return float(rec['Avg_Daily_Census_ICU'])/float(rec['Avg_Daily_Census'])

modelDir = '/home/welling/git/pyrhea/models/ChicagoLand'

allKeySet, recs = yaml_tools.parse_all(os.path.join(modelDir, 'facilityfacts'))
facDict = {r['abbrev']:r for r in recs}

outRecs = []
for abbrev, rec in facDict.items():
    oR = {'abbrev': abbrev}
    if rec['category'] == 'HOSPITAL':
        allBeds = rec['nBeds']['value']
        icuBeds = rec['nBedsICU']['value']
        nonICU = allBeds - icuBeds
        oR['HOSP'] = nonICU
        oR['ICU'] = icuBeds
    elif rec['category'] == 'LTACH':
        oR['LTACH'] = rec['nBeds']['value']
    elif rec['category'] == 'SNF':
        oR['NURSING'] = rec['nBeds']['value']
    elif rec['category'] == 'VSNF':
        nBeds = int(rec['nBeds']['value'])
        nBedsVent = int(nBeds * rec['fracVentBeds']['value'])
        nBedsSkil = int(nBeds * rec['fracSkilledBeds']['value'])
        nBedsOther = nBeds - (nBedsVent + nBedsSkil)
        oR['VENT'] = nBedsVent
        oR['SKILLED'] = nBedsSkil
        oR['NURSING'] = nBedsOther
    elif rec['category'] == 'COMMUNITY':
        pass
    else:
        raise RuntimeError('unknown category %s' % rec['category'])
    outRecs.append(oR)

with open('bed_tbl.csv', 'w') as f:
    csv_tools.writeCSV(f, ['abbrev', 'NURSING', 'SKILLED', 'VENT', 'LTACH', 'HOSP', 'ICU'],
                       outRecs)