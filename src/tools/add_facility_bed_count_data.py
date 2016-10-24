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

illinoisFName = 'Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016_BedLOS.csv'
ilInfo = loadCSVByAbbrev(modelDir, illinoisFName, key='UNIQUE_ID')
ilNBedsProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS:$S (includes ICU)"
ilNBedsICUProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS:$Q"
ilMeanPopProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS:$T (includes ICU)"
ilMeanPopICUProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS:$R"
ilMeanLOSICUProvStr = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS:$W"

remoteFName = 'PROTECT_WI_IN_Facility_BedSizes_091316_Facility.csv'
remInfo = loadCSVByAbbrev(modelDir, remoteFName, key='UNIQUE_ID')
remNBedsProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility:$O (includes ICU)'
remNBedsICUProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility:$N'
remMeanPopProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility Average_Daily_Census'
remMeanPopICUProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility Average_Daily_Census_ICU'
remMeanLOSICUProvStr = 'PROTECT_WI_IN_Facility_BedSizes_091316.xlsx:Facility LOS_ICU'

srcList = [(ilInfo, illinoisFName), (remInfo, remoteFName)]

nonHospTuples = [('nBeds', [('Peak_Beds', ilNBedsProvStr), ('Peak_Beds', remNBedsProvStr)]),
                 ('meanPop', [('Avg_Daily_Census', ilMeanPopProvStr),
                              ('Avg_Daily_Census', remMeanPopProvStr)])]

hospTuples = ([tpl for tpl in nonHospTuples][:] +
              [('nBedsICU', [('Peak_Beds_ICU', ilNBedsICUProvStr),
                             ('Peak_Beds_ICU', remNBedsICUProvStr)]),
               ('meanPopICU', [('Avg_Daily_Census_ICU', ilMeanPopICUProvStr),
                               ('Avg_Daily_Census_ICU', remMeanPopICUProvStr)]),
               ('meanLOSICU', [('LOS_ICU', ilMeanLOSICUProvStr),
                               ('LOS_ICU', remMeanLOSICUProvStr)])
              ])

nonHospEraseInvalidsList = ['nBeds', 'meanPop', 'meanLOSICU', 'meanLOS']
hospEraseInvalidsList = ['nBeds', 'meanPop', 'meanLOSICU', 'meanLOS']

nonHospEraseUnconditionallyList = []
hospEraseUnconditionallyList = []

nonHospUpdateUnconditionallyList = ['RESU_500_S', 'THC_2544_L', 'RESU_2380_S']
hospUpdateUnconditionallyList = ['ELMH_200_H', 'MERC_901_H','NORT_660_H', 'EDWA_801_H',
                                 'BETH_5025_H']

hospDoNotUpdateList = []
nonHospDoNotUpdateList = ['MIDW_8540_S', 'MCAL_18300_S']


newRecs = []
for abbrev, rec in facDict.items():
    nR = rec.copy()
    changed = False

    if rec['category'] in ['SNF', 'VSNF', 'NURSINGHOME', 'LTACH']:
        (eraseInvalidsList, eraseUnconditionallyList, updateUnconditionallyList, 
         doNotUpdateList, tupleList) = \
            (nonHospEraseInvalidsList, nonHospEraseUnconditionallyList,
             nonHospUpdateUnconditionallyList, nonHospDoNotUpdateList, nonHospTuples)
    elif rec['category'] in ['HOSPITAL']:
        (eraseInvalidsList, eraseUnconditionallyList, updateUnconditionallyList, 
         doNotUpdateList, tupleList) = \
            (hospEraseInvalidsList, hospEraseUnconditionallyList,
             hospUpdateUnconditionallyList, hospDoNotUpdateList, hospTuples)
    else:
        sys.exit('%s has unknown category %s' % (abbrev, rec['category']))

    for key in eraseInvalidsList:
        if key in nR and not typeCheck(nR[key]['value']):
            del nR[key]
            changed = True

    for key in eraseUnconditionallyList:
        if key in nR:
            del nR[key]
            changed = True

    for key, pairList in tupleList:

        if key in nR and typeCheck(nR[key]['value']):
            oldVal = nR[key]['value']
            srcFound = False
            for (info, srcFName), (srcKey, provStr) in zip(srcList, pairList):
                if abbrev in info and srcKey in info[abbrev] and typeCheck(info[abbrev][srcKey]):
                    srcFound = True
#                     if abbrev in ['ELMH_200_H', 'MERC_901_H','NORT_660_H', 'RESU_500_S', 'EDWA_801_H',
#                                   'BETH_5025_H', 'MIDW_8540_S', 'MCAL_18300_S']:
# #                     if False:
#                         nR[key] = {'value': info[abbrev][srcKey], 'prov': provStr}
#                         changed = True
#                     else:
                    if info[abbrev][srcKey] == oldVal:
                        nR[key]['prov'] = provStr
                        changed = True
                    else:
                        if abbrev in doNotUpdateList:
                            print ('skipping update of %s %s -> %s for %s' %
                                   (key, oldVal, info[abbrev][srcKey], abbrev))
                        elif abbrev in updateUnconditionallyList:
                            nR[key] = {'value': info[abbrev][srcKey], 'prov': provStr}
                            changed = True
                        else:
                            errStr = ('%s has value mismatch for %s with %s entry %s: %s vs %s' %
                                      (abbrev, key, srcFName, srcKey, oldVal, info[abbrev][srcKey]))
                            raise RuntimeError(errStr)
            if not srcFound:
                for testStr in ['median', 'mean', 'average', 'regression', 'interpolat']:
                    if testStr in nR[key]['prov'].lower():
                        srcFound = True
                        break
            if not srcFound:
                print 'No valid src for %s %s' % (abbrev, key)
        else:
            for (info, srcFName), (srcKey, provStr) in zip(srcList, pairList):
                if abbrev in info and srcKey in info[abbrev] and typeCheck(info[abbrev][srcKey]):
                    nR[key] = {'value': info[abbrev][srcKey], 'prov': provStr}
                    changed = True

    if changed:
        newRecs.append(nR)

print '%d of %d records modified' % (len(newRecs), len(recs))
yaml_tools.save_all(os.path.join(modelDir, 'facilityfactsUpdated'), newRecs)
