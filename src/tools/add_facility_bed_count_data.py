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

modelDir = '/home/welling/git/pyrhea/models/ChicagoLand'

allKeySet, recs = yaml_tools.parse_all(os.path.join(modelDir, 'facilityfacts'))
facDict = {r['abbrev']:r for r in recs}

illinoisFName = 'Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_v082116_BedLOS.csv'
ilInfo = loadCSVByAbbrev(modelDir, illinoisFName, key='UNIQUE_ID')
ilNBedsProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_v082116.xlsx:BedLOS:$S (includes ICU)"
ilNBedsICUProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_v082116.xlsx:BedLOS:$Q"
ilMeanPopProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_v082116.xlsx:BedLOS:$T (includes ICU)"
ilMeanPopICUProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_v082116.xlsx:BedLOS:$R"

remoteFName = 'PROTECT_WI_IN_Facility_BedSizes_091316_Facility.csv'
remInfo = loadCSVByAbbrev(modelDir, remoteFName, key='UNIQUE_ID')
remNBedsProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility:$O (includes ICU)'
remNBedsICUProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility:$N'
remMeanPopProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility Average_Daily_Census'
remMeanPopICUProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility Average_Daily_Census_ICU'

newRecs = []
for abbrev, rec in facDict.items():
    nR = rec.copy()
    changed = False
    
    for key in ['nBeds', 'meanPop', 'meanLOSICU']:
        if key in nR and not typeCheck(nR[key]['value']):
            del nR[key]
            changed = True
    
    if rec['category'] in ['SNF', 'VSNF', 'NURSINGHOME', 'LTACH']:
        if 'nBeds' in nR and typeCheck(nR['nBeds']['value']):
            oldVal = nR['nBeds']['value']
            if abbrev in ilInfo and typeCheck(ilInfo[abbrev]['Peak_Beds']):
                assert ilInfo[abbrev]['Peak_Beds'] == oldVal, ("%s has value mismatch with %s" %
                                                               (abbrev, illinoisFName))
            elif abbrev in remInfo:
                assert remInfo[abbrev]['Peak_Beds'] == oldVal, ('%s has a value mismatch with %s' %
                                                                (abbrev, remoteFName))
            else:
                pass
        elif abbrev in ilInfo and typeCheck(ilInfo[abbrev]['Peak_Beds']):
            nR['nBeds'] = {'value': ilInfo[abbrev]['Peak_Beds'],
                           'prov': ilNBedsProvStr}
            changed = True
        elif abbrev in remInfo and typeCheck(remInfo[abbrev]['Peak_Beds']):
            nR['nBeds'] = {'value': remInfo[abbrev]['Peak_Beds'],
                           'prov': remNBedsProvStr}
            changed = True
                
        if 'meanPop' in nR and typeCheck(nR['meanPop']['value']):
            oldVal = nR['meanPop']['value']
            if abbrev in ilInfo:
                assert ilInfo[abbrev]['Avg_Daily_Census'] == oldVal, ("%s has value mismatch with %s" %
                                                                      (abbrev, illinoisFName))
            elif abbrev in remInfo:
                assert remInfo[abbrev]['Avg_Daily_Census'] == oldVal, ('%s has a value mismatch with %s' %
                                                                       (abbrev, remoteFName))
            else:
                pass
        elif abbrev in ilInfo and typeCheck(ilInfo[abbrev]['Avg_Daily_Census']):
            nR['meanPop'] = {'value': ilInfo[abbrev]['Avg_Daily_Census'],
                             'prov': ilMeanPopProvStr}
            changed = True
        elif (abbrev in remInfo and 'Avg_Daily_Census' in remInfo[abbrev] and
              typeCheck(remInfo[abbrev]['Avg_Daily_Census'])):
            nR['meanPop'] = {'value': remInfo[abbrev]['Avg_Daily_Census'],
                             'prov': remMeanPopProvStr}
            changed = True
                
    elif rec['category'] in ['HOSPITAL']:
        if 'nBeds' in nR and typeCheck(nR['nBeds']['value']):
            oldVal = nR['nBeds']['value']
            if abbrev in ilInfo and typeCheck(ilInfo[abbrev]['Peak_Beds']):
                assert ilInfo[abbrev]['Peak_Beds'] == oldVal, ("%s has value mismatch with %s" %
                                                               (abbrev, illinoisFName))
            elif abbrev in remInfo:
                assert remInfo[abbrev]['Peak_Beds'] == oldVal, ('%s has a value mismatch with %s' %
                                                                (abbrev, remoteFName))
            else:
                pass
        elif abbrev in ilInfo and typeCheck(ilInfo[abbrev]['Peak_Beds']):
            nR['nBeds'] = {'value': ilInfo[abbrev]['Peak_Beds'],
                           'prov': ilNBedsProvStr}
            changed = True
        elif abbrev in remInfo and typeCheck(remInfo[abbrev]['Peak_Beds']):
            nR['nBeds'] = {'value': remInfo[abbrev]['Peak_Beds'],
                           'prov': remNBedsProvStr}
            changed = True

        if 'nBedsICU' in nR and typeCheck(nR['nBedsICU']['value']):
            oldVal = nR['nBedsICU']['value']
            if abbrev in ilInfo and typeCheck(ilInfo[abbrev]['Peak_Beds_ICU']):
                assert ilInfo[abbrev]['Peak_Beds_ICU'] == oldVal, ("%s has value mismatch with %s" %
                                                                   (abbrev, illinoisFName))
            elif abbrev in remInfo:
                assert remInfo[abbrev]['Peak_Beds_ICU'] == oldVal, ('%s has a value mismatch with %s' %
                                                                    (abbrev, remoteFName))
            else:
                pass
        elif abbrev in ilInfo and typeCheck(ilInfo[abbrev]['Peak_Beds_ICU']):
            nR['nBedsICU'] = {'value': ilInfo[abbrev]['Peak_Beds_ICU'],
                              'prov': ilNBedsICUProvStr}
            changed = True
        elif abbrev in remInfo and typeCheck(remInfo[abbrev]['Peak_Beds_ICU']):
            nR['nBedsICU'] = {'value': remInfo[abbrev]['Peak_Beds_ICU'],
                              'prov': remNBedsICUProvStr}
            changed = True
                
        if 'meanPop' in nR and typeCheck(nR['meanPop']['value']):
            oldVal = nR['meanPop']['value']
            if abbrev in ilInfo:
                assert ilInfo[abbrev]['Avg_Daily_Census'] == oldVal, ("%s has value mismatch with %s" %
                                                                      (abbrev, illinoisFName))
            elif abbrev in remInfo:
                assert remInfo[abbrev]['Avg_Daily_Census'] == oldVal, ('%s has a value mismatch with %s' %
                                                                       (abbrev, remoteFName))
            else:
                pass
        elif abbrev in ilInfo and typeCheck(ilInfo[abbrev]['Avg_Daily_Census']):
            nR['meanPop'] = {'value': ilInfo[abbrev]['Avg_Daily_Census'],
                             'prov': ilMeanPopProvStr}
            changed = True
        elif (abbrev in remInfo and 'Avg_Daily_Census' in remInfo[abbrev] and
              typeCheck(remInfo[abbrev]['Avg_Daily_Census'])):
            nR['meanPop'] = {'value': remInfo[abbrev]['Avg_Daily_Census'],
                             'prov': remMeanPopProvStr}
            changed = True

        if 'meanPopICU' in nR and typeCheck(nR['meanPopICU']['value']):
            oldVal = nR['meanPopICU']['value']
            if abbrev in ilInfo:
                assert ilInfo[abbrev]['Avg_Daily_Census_ICU'] == oldVal, ("%s has value mismatch with %s" %
                                                                          (abbrev, illinoisFName))
            elif abbrev in remInfo and 'meanPopICU' in remInfo[abbrev]:
                assert remInfo[abbrev]['Avg_Daily_Census_ICU'] == oldVal, ('%s has a value mismatch with %s' %
                                                                           (abbrev, remoteFName))
            else:
                pass
        elif abbrev in ilInfo and typeCheck(ilInfo[abbrev]['Avg_Daily_Census_ICU']):
            nR['meanPopICU'] = {'value': ilInfo[abbrev]['Avg_Daily_Census_ICU'],
                                'prov': ilMeanPopICUProvStr}
            changed = True
        elif (abbrev in remInfo and 'Avg_Daily_Census_ICU' in remInfo[abbrev] and
              typeCheck(remInfo[abbrev]['Avg_Daily_Census_ICU'])):
            nR['meanPopICU'] = {'value': remInfo[abbrev]['Avg_Daily_Census_ICU'],
                                'prov': remMeanPopICUProvStr}
            changed = True
    
    else:
        sys.exit('%s has unknown category %s' % (abbrev, rec['category']))
        
    if changed:
        newRecs.append(nR)
        
                
print '%d of %d records modified' % (len(newRecs), len(recs))
yaml_tools.save_all(os.path.join(modelDir, 'facilityfactsUpdated'), newRecs)

    
    

