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

illinoisFName = 'Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016_BedLOS.csv'
ilInfo = loadCSVByAbbrev(modelDir, illinoisFName, key='UNIQUE_ID')
ilNBedsProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS:$S (includes ICU) 14ba183e"
ilNBedsICUProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS:$Q 14ba183e"
ilMeanPopProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS:$T (includes ICU) 14ba183e"
ilMeanPopICUProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS:$R 14ba183e"
ilMeanLOSICUProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS:$W 14ba183e"
ilHospFracICUProvStr = '%s 14ba183e Avg_Daily_Census_ICU/Avg_Daily_Census' % illinoisFName
ilMeanLOSProvStr = '%s 14ba183e LOS' % illinoisFName



remoteFName = 'PROTECT_WI_IN_Facility_BedSizes_091316_Facility.csv'
remInfo = loadCSVByAbbrev(modelDir, remoteFName, key='UNIQUE_ID')
remNBedsProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility:$O (includes ICU) 14ba183e'
remNBedsICUProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility:$N 14ba183e'
remMeanPopProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility Average_Daily_Census 14ba183e'
remMeanPopICUProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility Average_Daily_Census_ICU 14ba183e'
remMeanLOSICUProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility LOS_ICU 14ba183e'

matLosFName = 'Matrices_LOS_09292016_Facilities_LOS.csv'
matLosInfo = loadCSVByAbbrev(modelDir, matLosFName, key='UNIQUE_ID')
mlMeanLOSProvStr = '%s rev 14ba183e Mean' % matLosFName


hospSet = set(['HOSPITAL'])
nonHospSet = set(['SNF', 'VSNF', 'NURSINGHOME', 'LTACH'])
anySet = hospSet.union(nonHospSet)

"""
Tuple format is (restrictionSet, key, srcKey, srcProv, srcDict, srcFName)
"""
tupleList = [(anySet, 'nBeds', [('Peak_Beds', ilNBedsProvStr, ilInfo, illinoisFName),
                                ('Peak_Beds', remNBedsProvStr, remInfo, remoteFName)]),
             (anySet, 'meanPop', [('Avg_Daily_Census', ilMeanPopProvStr, ilInfo, illinoisFName),            
                                  ('Avg_Daily_Census', remMeanPopProvStr, remInfo, remoteFName)]),
             (hospSet, 'nBedsICU', [('Peak_Beds_ICU', ilNBedsICUProvStr, ilInfo, illinoisFName),          
                                    ('Peak_Beds_ICU', remNBedsICUProvStr, remInfo, remoteFName)]),          
             (hospSet, 'meanPopICU', [('Avg_Daily_Census_ICU', ilMeanPopICUProvStr, ilInfo, illinoisFName),          
                                      ('Avg_Daily_Census_ICU', remMeanPopICUProvStr, remInfo, remoteFName)]),          
             (hospSet, 'meanLOSICU', [('LOS_ICU', ilMeanLOSICUProvStr, ilInfo, illinoisFName),          
                                      ('LOS_ICU', remMeanLOSICUProvStr, remInfo, remoteFName)]),
             (hospSet, 'fracAdultPatientDaysICU', [(hospFracICUOp, ilHospFracICUProvStr, 
                                                    ilInfo, illinoisFName)]),
             (anySet, 'meanLOS', [('LOS', ilMeanLOSProvStr, ilInfo, illinoisFName),
                                  ('Mean', mlMeanLOSProvStr, matLosInfo, matLosFName)])
             ]

eraseInvalidsList = ['nBeds', 'meanPop', 'meanLOSICU', 'meanLOS', 'meanPopICU', 'nBedsICU',
                     'fracAdultPatientDaysICU']

eraseUnconditionallyList = []

#updateUnconditionallyList = ['RESU_500_S', 'THC_2544_L', 'RESU_2380_S', 'BETH_5025_H', 'MERC_901_H']
updateUnconditionallyList = ['ELMH_200_H', 'MERC_901_H', 'NORT_660_H', 'RESU_500_S', 'EDWA_801_H',
                             'THC_2544_L', 'BETH_5025_H', 'RESU_2380_S']

#doNotUpdateList = ['MIDW_8540_S', 'MCAL_18300_S', '']
doNotUpdateList = []

newRecs = []
for abbrev, rec in facDict.items():
    nR = rec.copy()
    changed = False

    for key in eraseInvalidsList:
#         if key in nR and not typeCheck(nR[key]['value']):
        if key in nR and ((not typeCheck(nR[key]['value']))
                          or 'v082116' in nR[key]['prov']
                          or 'v100716' in nR[key]['prov']
                          or 'v102016' in nR[key]['prov']):
            print 'ping %s %s' % (abbrev, key)
            del nR[key]
            changed = True

    for key in eraseUnconditionallyList:
        if key in nR:
            del nR[key]
            changed = True
    
    for restrictionSet, key, srcList in tupleList:
        if key == 'meanLOS' and key in nR and 'Avg_Daily_Census' in nR[key]['prov']:
            print 'Nonsense %s' % abbrev
            del nR[key]
            changed = True

        if rec['category'] not in restrictionSet:
            continue
        oldVal = None
        oldProv = None
        newVal = None
        newProv = None
        srcFound = False
        for srcKey, srcProv, info, srcFName in srcList:
            if key in nR and typeCheck(nR[key]['value']):
                oldVal = nR[key]['value']
                oldProv = nR[key]['prov']
            if abbrev in info:
                if isinstance(srcKey, (types.FunctionType, types.LambdaType)):
                    try:
                        newVal = srcKey(info[abbrev])
                        srcFound = True
                    except ValueError, e:
                        print '%s %s %s: %s' % (abbrev, key, srcFName, e)
                elif (srcKey in info[abbrev] and typeCheck(info[abbrev][srcKey]) and
                      info[abbrev][srcKey] != 0):
                    newVal = info[abbrev][srcKey]
                    srcFound = True
                if srcFound:
                    newProv = srcProv
                    break
        if abbrev == 'MCAL_18300_S': print '%s %s %s %s' % (abbrev, key, oldVal, newVal)
        if oldVal:
            if newVal:
                assert newProv, ('Update %s %s -> %s for %s has no prov string?' %
                                 (key, oldVal, newVal, abbrev))
                if newVal == oldVal:
                    if newProv != oldProv:
                        nR[key]['prov'] = newProv
                        changed = True
                else:
                    if abbrev in doNotUpdateList:
                        print ('skipping update of %s %s -> %s for %s' %
                               (key, oldVal, newVal, abbrev))
                    elif abbrev in updateUnconditionallyList:
                        nR[key] = {'value': newVal, 'prov': newProv}
                        changed = True
                    else:
                        errStr = ('%s has value mismatch for %s with %s entry %s: %s vs %s' %
                                  (abbrev, key, srcFName, srcKey, oldVal, newVal))
                        raise RuntimeError(errStr)
            if not srcFound:
                for testStr in ['median', 'mean', 'average', 'regression', 'interpolat', 'line ']:
                    if testStr in nR[key]['prov'].lower():
                        srcFound = True
                        break
            if not srcFound:
                print 'No valid src for %s %s' % (abbrev, key)
        elif newVal:
            assert newProv, ('Update %s %s -> %s for %s has no prov string?' %
                             (key, oldVal, newVal, abbrev))
            nR[key] = {'value': newVal, 'prov': newProv}
            changed = True

    if changed:
        newRecs.append(nR)

print '%d of %d records modified' % (len(newRecs), len(recs))
yaml_tools.save_all(os.path.join(modelDir, 'facilityfactsUpdated'), newRecs)
