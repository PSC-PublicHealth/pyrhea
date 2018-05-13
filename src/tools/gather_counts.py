#! /usr/bin/env python

###################################################################################
# Copyright   2017, Pittsburgh Supercomputing Center (PSC).  All Rights Reserved. #
# =============================================================================== #
#                                                                                 #
# Permission to use, copy, and modify this software and its documentation without #
# fee for personal use within your organization is hereby granted, provided that  #
# the above copyright notice is preserved in all copies and that the copyright    #
# and this permission notice appear in supporting documentation.  All other       #
# restrictions and obligations are defined in the GNU Affero General Public       #
# License v3 (AGPL-3.0) located at http://www.gnu.org/licenses/agpl-3.0.html  A   #
# copy of the license is also provided in the top level of the source directory,  #
# in the file LICENSE.txt.                                                        #
#                                                                                 #
###################################################################################

import csv
import sys
import os
import glob
from multiprocessing import Process,Manager,cpu_count,Pool
from optparse import OptionParser
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy.stats as st
import yaml

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))
import pyrheautils
import schemautils
import pathogenbase as pth
from facilitybase import CareTier, PthStatus
from notes_plotter import readFacFiles, scanAllFacilities, checkInputFileSchema
from notes_plotter import importNotes
from notes_plotter import SCHEMA_DIR, INPUT_SCHEMA, CARE_TIERS, FAC_TYPE_TO_CATEGORY_MAP
from time_series_plotter import mergeNotesFiles, getTimeSeriesList
from pyrhea import getLoggerConfig
import time
#import affinity

#print "Cpu_Count = {0}".format(cpu_count())
#sys.stdout.flush()
#if cpu_count() < 60:
 #   affinity.set_process_affinity_mask(0,2**cpu_count()-1)

#os.system("taskset -p 0xff %d"%os.getpid())

DEFAULT_OUT_FILE = 'counts_output.yaml'


def extractCountsFromNotes(note, abbrevList, facDict, translationDict, burninDays):
    print "note = {0}".format(note)
    returnDict = {}
    ### totalStats is assumed to be a Manager
    try:
        specialDict = mergeNotesFiles([note],False)
    except Exception as e:
        print "for file {0} there is an exception {1}".format(note,e)
        return {}
    print "here {0}".format(note)
    for abbrev in abbrevList:
        returnDict[abbrev] = {}
        tplList = getTimeSeriesList(abbrev,specialDict,'localtierpathogen')
        for dayVec, curves in tplList:
            #print dayVec
            dayIndex = np.where(dayVec==burninDays)[0][0]
            #print dayIndex
            tmpVecD = defaultdict(list)
            totVecD = {}
            for tpl, lVec in curves.items():
                tier, pthStatus = tpl
                tmpVecD[tier].append(lVec)
            for tier, lVecList in tmpVecD.items():
                totVecD[tier] = sum(lVecList)
            found = False
            for tpl, lVec in curves.items():
                tier, pthStatus = tpl
                if pthStatus == PthStatus.COLONIZED:
                    #print len(lVec[dayIndex:])
                    returnDict[abbrev][CareTier.names[tier]] = {'colonizedDays': np.sum(lVec[dayIndex:]),
                                                                'bedDays':np.sum(totVecD[tier][dayIndex:]),
                                                                'colonizedDaysTS':lVec[dayIndex:],
                                                                'bedDaysTS':totVecD[tier][dayIndex:]}

                    for key in translationDict.keys():
                        returnDict[abbrev][CareTier.names[tier]][key] = 0.0

    for key,noteKey in translationDict.items():
        print "key = {0}".format(key)
        for abbrev in abbrevList:
            tplList = getTimeSeriesList(abbrev, specialDict, noteKey)
            for dayVec, curves in tplList:
                dayIndex = np.where(dayVec==(burninDays+1))[0][0]
                for tpl, curve in curves.items():
                    returnDict[abbrev][CareTier.names[tpl]][key] = np.sum(curve[dayIndex:])
                    returnDict[abbrev][CareTier.names[tpl]]["{0}TS".format(key)] = curve[dayIndex:]

    return returnDict


def combineTimeSeries(tsArray,runDays, tsIncrement=1):
    def ci(x):
        ci_tuple = st.t.interval(0.95,len(x)-1, loc = np.mean(x), scale=st.sem(x))
        return pd.Series(ci_tuple)
#
    #DataFrame({'5per':ci_tuple[0],'95per':ci_tuple[1]})

    returnDict = {'mean':np.zeros(runDays),'median':np.zeros(runDays),
                  'stdv':np.zeros(runDays),'5%CI':np.zeros(runDays),'95%CI':np.zeros(runDays)}


    df = pd.DataFrame(tsArray).transpose()

    #df2 = df.mean(axis=1)

    #print tsArray
    returnDict['mean'] = df.mean(axis=1)
    returnDict['median'] = df.median(axis=1)
    returnDict['stdv'] = df.std(axis=1)

    n = len(df.columns)
    rf = pd.DataFrame({'sem':(df.std(axis=1)),'mean':df.mean(axis=1)})
    i95 = st.t._ppf(1.95/2.0, n-1)
    rf['c95_lower'] = rf['mean'] - (rf['sem'] * i95)
    rf['c95_upper'] = rf['mean'] + (rf['sem'] * i95)

    returnDict['5%CI'] = rf['c95_lower']
    returnDict['95%CI'] = rf['c95_upper']
    #ciTmp = df.apply(lambda x: pd.Series(sms.DescrStatsW(x).tconfint_mean(), index=['ci5','ci95']), axis=1)
    #returnDict['5%CI'] = np.zeciTmp.ci5
    #returnDict['95%CI'] = ciTmp.ci95
    
#     semArray = np.apply_along_axis(st.sem, 0, tsArray)
#     #[0.0 for i in range(0,runDays)]
#     #returnDict['95%CI'] = [0.0 for i in range(0,runDays)]
#     lenArray= np.array(tsArray.shape[1])
#     
#     print "len = {0}".format(lenArray)
#     for i in range(0,runDays):
#         #print "array = {0}".format([x[i] for x in tsArray])
#         #print "mean = {0}".format(np.mean([x[i] for x in tsArray]))
# #        returnDict['mean'].append(np.mean([x[i] for x in tsArray]))
# #        returnDict['median'].append(np.median([x[i] for x in tsArray]))
# #        returnDict['stdv'].append(np.std([x[i] for x in tsArray]))
#         ci = st.t.interval(0.95, lenArray - 1, loc=returnDict['mean'][i], scale=semArray[i])
#         returnDict['5%CI'].append(ci[0])
#         returnDict['95%CI'].append(ci[1])
#         
#     #print "returnDict mean = {0}".format(returnDict['mean'])

    return returnDict


def pool_helper(args):
    return extractCountsFromNotes(*args)

def main():
    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl [--out outname.yaml] run_descr.yaml
    """)
    
    parser.add_option('-n', '--notes', action='append', type='string',
                     help="Notes filename - may be repeated")
    parser.add_option('-o', '--out', action='store', type='string',
                     help="Output file (defaults to %s)" % DEFAULT_OUT_FILE)
    parser.add_option('-g','--glob', action='store_true',
                      help=("Apply filename globbing for notes files."
                            "  (Remember to protect the filename string from the shell!)"))
    parser.add_option('-x','--xdroyamlfile',type='string', default=None,
                      help="specify the xdro scenario yaml file if there is one, if not, XDRO won't be costed")
    parser.add_option('-m','--nprocs',type='int',default=1,
                      help='number of cpus to run the costing model over')
    parser.add_option('-t','--producetimeseries',action='store_true',default=False)


    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A run file is required.  Try "-h"')

    inputDict = checkInputFileSchema(args[0], os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    if opts.out:
        outFileName = opts.out
    else:
        outFileName = DEFAULT_OUT_FILE

    if 'trackedFacilities' not in inputDict or not len(inputDict['trackedFacilities']):
        raise RuntimeError('Run description file does not list any tracked facilities')

    pyrheautils.prepPathTranslations(inputDict)
    if 0:
        modelDir = inputDict['modelDir']
        pyrheautils.PATH_STRING_MAP['MODELDIR'] = os.path.abspath(modelDir)
        implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])
        pyrheautils.PATH_STRING_MAP['IMPLDIR'] = implDir
        if 'pathTranslations' in inputDict:
            for elt in inputDict['pathTranslations']:
                pyrheautils.PATH_STRING_MAP[elt['key']] = elt['value']

    facDirList = [pyrheautils.pathTranslate(fPth) for fPth in inputDict['facilityDirs']]
    facDict = readFacFiles(facDirList)

    ### get the abbreviations within a 13 mile radius
    facIn13Path = pyrheautils.pathTranslate("$(CONSTANTS)/facilities_in_13miles.yaml")
    with open(facIn13Path, "rb") as f:
        facilitiesWithin13Miles = yaml.load(f)['facilitiesWithin13Miles']['locAbbrevList']

    facInCookPath = pyrheautils.pathTranslate("$(CONSTANTS)/facilities_in_CookCounty.yaml")
    with open(facInCookPath, "rb") as f:
        facilitiesWithinCookCounty = yaml.load(f)['facilitiesWithinCookCounty']['locAbbrevList']

    burninDays = int(inputDict['burnInDays'])
    print "burninDays = {0}".format(burninDays)
    runDays = int(inputDict['runDurationDays'])
    print "runDays = {0}".format(runDays)

    ### Translation Dict
    valuesToGatherList = ['newColonized', 'creArrivals', 'arrivals', 'contactPrecautionDays',
                          'creBundlesHandedOut','creSwabsUsed','newPatientsOnCP','passiveCPDays',
                          'swabCPDays','xdroCPDays','otherCPDays']

    valuesToGather = {'newColonized':'localtiernewcolonized',
                      'creArrivals':'localtiercrearrivals',
                      'arrivals':'localtierarrivals',
                      'contactPrecautionDays': 'localtierCP',
                      'creBundlesHandedOut': 'localtierCREBundle',
                      'creSwabsUsed':'localtierCRESwabs',
                      'newPatientsOnCP': 'localtierpatientsOnCP',
                      'passiveCPDays':'localtierpassiveCP',
                      'swabCPDays':'localtierswabCP',
                      'xdroCPDays':'localtierxdroCP',
                      'otherCPDays':'localtierotherCP',
                      }

    tableHeadings = {'newColonized':'Newly Colonized',
                      'creArrivals':'CRE Colonized Patients Admissions',
                      'arrivals':'Patient Admissions',
                      'contactPrecautionDays': 'Contact Precaution Days',
                      'creBundlesHandedOut': 'CRE Baths Given Out',
                      'creSwabsUsed':'CRE Swabs Used',
                      'newPatientsOnCP':'Number of Patients Put on CP',
                      'passiveCPDays':'CRE CP Days due to passive surveillance',
                      'swabCPDays':'CRE CP Days due to acitve surveillance',
                      'xdroCPDays':'CRE CP Days due to xdro registry',
                      'otherCPDays':'CP Days for other reasons'
                      }

    notes = []
    if opts.glob:
        notes = glob.glob('{0}'.format(opts.notes[0]))
    else:
        notes = opts.notes

    abbrevList = inputDict['trackedFacilities']
    nprocs = opts.nprocs
    print "nprocs = {0}".format(nprocs)
    print "notes = {0}".format(notes)

    argsList = [(notes[i], abbrevList, facDict, valuesToGather, burninDays)
                for i in range(0, len(notes))]

    xdroAbbrevs = []
    if opts.xdroyamlfile:
        with open("{0}".format(opts.xdroyamlfile),"rb") as f:
            xdroparams = yaml.load(f)
            xdroAbbrevs = xdroparams['locationsImplementingScenario']['locAbbrevList']

    print "XDRO Facilities = {0}".format(xdroAbbrevs)

    if nprocs > 1:
        p = Pool(nprocs)
        totalStats = p.map(pool_helper, argsList)
        p.close()
    else:
        print 'Warning: using serial processing because nprocs == 1'
        totalStats = []
        for args in argsList:
            totalStats.append(pool_helper(args))
    print 'Finished scanning notes'

    numNotesFiles = len(totalStats)
    tsDFL = []
    scalarDL = []
    for idx, dct in enumerate(totalStats):
        for abbrev, subD in dct.items():
            for tier, subSubD in subD.items():
                
                tsD = {}
                scalarD = {}
                for key, val in subSubD.items():
                    if key.endswith('TS'):
                        tsD[key[:-2]] = val
                    else:
                        scalarD[key] = val
                maxLen = max([len(val) for val in tsD.values()])
                newTSD = {}
                for key, vec in tsD.items():
                    if len(vec) < maxLen:
                        newTSD[key] = np.pad(vec, (maxLen - len(vec), 0),
                                             'constant', constant_values=(0, 0))
                    else:
                        newTSD[key] = vec
                df = pd.DataFrame.from_dict(newTSD)
                df['abbrev'] = abbrev
                df['tier'] = tier
                df['run'] = idx
                df.index.name = 'day'
                tsDFL.append(df)
                scalarD['abbrev'] = abbrev
                scalarD['tier'] = tier
                scalarD['run'] = idx
                scalarDL.append(scalarD)
    tsDF = pd.concat(tsDFL).reset_index()
    scalarDF = pd.DataFrame(scalarDL)
    #print tsDF
    #print scalarDF


    print "totalStats = {0}".format(totalStats)
    tsDF.to_msgpack('/tmp/total_counts.mpz')

    ### By Tier
    ### Each of these should be the same in terms of the abbrevs and tiers, so we can use the first to index the rest

    print "Processing Outputs"

    sys.stdout.flush()
    statsByTier = {}
    statsByAbbrev = {}
    time1Counter = 0.0
    time2Counter = 0.0
    for abbrev,tD in totalCounts[0].items():
        if abbrev not in statsByTier.keys():
            statsByTier[abbrev] = {}
            statsByAbbrev[abbrev] = {}
            statDict = {k:{'value': np.array([0.0 for x in range(0,numNotesFiles)]),
                           'ts': [[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]}
                           for k in valuesToGather.keys()}
            #statDictTS = {k:{'value': np.array([0.0 for x in range(0,numNotesFiles)]),
            #               'ts': [[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]} for k in valuesToGather.keys()}
            #statDictTS = {k:[[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]}
            cDAs = np.array([0.0 for x in range(0,numNotesFiles)])
            bDAs = np.array([0.0 for x in range(0,numNotesFiles)])
            prevAs = np.array([0.0 for x in range(0,numNotesFiles)])
            
            
            if opts.producetimeseries:
                cTAs = [[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]
                bTAs = [[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]
          
        for tier,d in tD.items():
            time1 = time.time()
            if tier not in statsByTier[abbrev].keys():
                statsByTier[abbrev][tier] = {}
                #statsByTierTS[abbrev][tier] = {}
                statTierDict = {k:{'value': np.array([0.0 for x in range(0,numNotesFiles)]),
                                   'ts': [[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]} for k in valuesToGather.keys()}
                #statTierDict = {k:np.array([0.0 for x in range(0,numNotesFiles)]) for k in valuesToGather.keys()}
                #statTierDictTS = {k:[[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]} 
            cDList = []
            bDList = []
            if opts.producetimeseries:
                cTs = []
                bTs = []
            for k in valuesToGather.keys():
                statTierDict[k]['value'] = []
                if opts.producetimeseries:
                    statTierDict[k]['ts'] = []

            for x in totalCounts:
                try:
                    cDList.append(float(x[abbrev][tier]['colonizedDays']))
                    bDList.append(float(x[abbrev][tier]['bedDays']))
                    if opts.producetimeseries:
                        cTs.append(x[abbrev][tier]['colonizedDaysTS'])
                        bTs.append(x[abbrev][tier]['bedDaysTS'])

                    for k in valuesToGather.keys():
                        statTierDict[k]['value'].append(float(x[abbrev][tier][k]))
                        if opts.producetimeseries:
                            statTierDict[k]['ts'].append(x[abbrev][tier]["{0}TS".format(k)])
                        
                except:
                    cDList.append(0.0)
                    bDList.append(0.0)
                    if opts.producetimeseries:
                        cTList.append(0)
                        bTList.append(0)
                    for k in valuesToGather.keys():
                        statTierDict[k]['value'].append(0.0)
                        if opts.producetimeseries:
                            statTierDict[k]['ts'].append(0)
 
            cDs = np.array(cDList)
            bDs = np.array(bDList)
            for k in valuesToGather.keys():
                statTierDict[k]['value'] = np.array(statTierDict[k]['value'])


#            cDs = np.array([float(x[abbrev][tier]['colonizedDays']) for x in totalCounts])
#            bDs = np.array([float(x[abbrev][tier]['bedDays']) for x in totalCounts])
#            if opts.producetimeseries:
#                cTs = [ x[abbrev][tier]['colonizedDaysTS'] for x in totalCounts ]
#                bTs = [ x[abbrev][tier]['bedDaysTS'] for x in totalCounts ]
#                #print "bTs = {0}".format(bTs)
#            
#            for k in valuesToGather.keys():
#                #if k == "newColonized":
#                #    print np.array([float(x[abbrev][tier][k]) for x in totalCounts])
#                statTierDict[k]['value'] = np.array([float(x[abbrev][tier][k]) for x in totalCounts])
#                if opts.producetimeseries:
#                    statTierDict[k]['ts'] = [ x[abbrev][tier]["{0}TS".format(k)] for x in totalCounts ]

            #ncDs = np.array([float(x[abbrev][tier]['newColonized']) for x in totalCounts])
            #caDs = np.array([float(x[abbrev][tier]['creArrivals']) for x in totalCounts])
            #aDs = np.array([float(x[abbrev][tier]['arrivals']) for x in totalCounts])
            prevs = []
            for i in range(0,len(cDs)):
                prevs.append((cDs[i]/bDs[i]))
                #print cTs
                if opts.producetimeseries:
                    for j in range(0,runDays):
    #                    print j
    #                    print cTs[i][j]
                        cTAs[i][j] += cTs[i][j]
                        bTAs[i][j] += bTs[i][j]
                cDAs[i] += cDs[i]
                bDAs[i] += bDs[i]
                for k in valuesToGather.keys():
                    statDict[k]['value'][i] += statTierDict[k]['value'][i]
                    if opts.producetimeseries:
                        for j in range(0,runDays):
                            statDict[k]['ts'][i][j] += statTierDict[k]['ts'][i][j]
#                 ncAs[i] += ncDs[i]
#                 caDAs[i] += caDs[i]
#                 aDAs[i] += aDs[i]
#                 
            #print "{0} {1} {2}".format(abbrev,tier,cDs)
            #print statsByTier  
            
            time2 = time.time()
            time1Counter += time2-time1
            
            time1 = time.time()
            statsByTier[abbrev][tier]['colonizedDays'] = {'mean':np.mean(cDs),
                                                          'median':np.median(cDs),
                                                          'stdv':np.std(cDs),
                                                          '5%CI':st.t.interval(0.95,len(cDs)-1, loc=np.mean(cDs),scale=st.sem(cDs))[0],
                                                          '95%CI':st.t.interval(0.95,len(cDs)-1, loc=np.mean(cDs),scale=st.sem(cDs))[1]}
            
            statsByTier[abbrev][tier]['bedDays'] = {'mean':np.mean(bDs),
                                                    'median':np.median(bDs),
                                                    'stdv':np.std(bDs),
                                                    '5%CI':st.t.interval(0.95,len(bDs)-1, loc=np.mean(bDs),scale=st.sem(bDs))[0],
                                                    '95%CI':st.t.interval(0.95,len(bDs)-1, loc=np.mean(bDs),scale=st.sem(bDs))[1]}
            
            statsByTier[abbrev][tier]['prevalence'] = {'mean':np.mean(prevs),
                                                       'median':np.median(prevs),
                                                       'stdv':np.std(prevs),
                                                       '5%CI':st.t.interval(0.95,len(prevs)-1, loc=np.mean(prevs),scale=st.sem(prevs))[0],
                                                       '95%CI':st.t.interval(0.95,len(prevs)-1, loc=np.mean(prevs),scale=st.sem(prevs))[1]}
            if opts.producetimeseries:
                statsByTier[abbrev][tier]['colonizedDaysTS'] = combineTimeSeries(cTs,runDays,1)
                statsByTier[abbrev][tier]['bedDaysTS'] = combineTimeSeries(bTs,runDays,1)
            for k in valuesToGather.keys():
                #print "Key 2 = {0}".format(k)
                statsByTier[abbrev][tier][k] = {'mean':np.mean(statTierDict[k]['value']),
                                                       'median':np.median(statTierDict[k]['value']),
                                                       'stdv':np.std(statTierDict[k]['value']),
                                                       '5%CI':st.t.interval(0.95, len(statTierDict[k]['value']) - 1, loc=np.mean(statTierDict[k]['value']), scale=st.sem(statTierDict[k]['value']))[0],
                                                       '95%CI':st.t.interval(0.95, len(statTierDict[k]['value']) - 1, loc=np.mean(statTierDict[k]['value']), scale=st.sem(statTierDict[k]['value']))[1]}
                if opts.producetimeseries:
                    statsByTier[abbrev][tier]["{0}TS".format(k)] = combineTimeSeries(statTierDict[k]['ts'],runDays,1)
            
            if abbrev in xdroAbbrevs:
                statsByTier[abbrev][tier]["xdroAdmissions"] = {'mean':np.mean(statTierDict['arrivals']['value']),
                                                               'median':np.median(statTierDict['arrivals']['value']),
                                                               'stdv':np.std(statTierDict['arrivals']['value']),
                                                               '5%CI':st.t.interval(0.95, len(statTierDict['arrivals']['value']) - 1, loc=np.mean(statTierDict['arrivals']['value']), scale=st.sem(statTierDict['arrivals']['value']))[0],
                                                               '95%CI':st.t.interval(0.95, len(statTierDict['arrivals']['value']) - 1, loc=np.mean(statTierDict['arrivals']['value']), scale=st.sem(statTierDict['arrivals']['value']))[1]}
                if opts.producetimeseries:
                    statsByTier[abbrev][tier]['xdroAdmissionsTS'] = combineTimeSeries(statTierDict['arrivals']['ts'],runDays,1)
            
            else:
                statsByTier[abbrev][tier]["xdroAdmissions"] = {'mean':0.0,
                                                               'median':0.0,
                                                               'stdv':0.0,
                                                               '5%CI':0.0,
                                                               '95%CI':0.0}
                if opts.producetimeseries:
                    statsByTier[abbrev][tier]['xdroAdmissionsTS'] = {'mean':[0.0 for x in range(0,runDays)],
                                                                     'median':[0.0 for x in range(0,runDays)],
                                                                     'stdv':[0.0 for x in range(0,runDays)],
                                                                     '5%CI':[0.0 for x in range(0,runDays)],
                                                                     '95%CI':[0.0 for x in range(0,runDays)]}
            time2 = time.time()
            time2Counter += time2-time1
            
        for i in range(0,len(cDAs)):
            prevAs[i] = (cDAs[i]/bDAs[i])
        
        statsByAbbrev[abbrev]['colonizedDays'] = {'mean':np.mean(cDAs),
                                                  'median':np.median(cDAs),
                                                  'stdv':np.std(cDAs),
                                                  '5%CI':st.t.interval(0.95,len(cDAs)-1, loc=np.mean(cDAs),scale=st.sem(cDAs))[0],
                                                  '95%CI':st.t.interval(0.95,len(cDAs)-1, loc=np.mean(cDAs),scale=st.sem(cDAs))[1]}
        
        statsByAbbrev[abbrev]['bedDays'] = {'mean':np.mean(bDAs),
                                          'median':np.median(bDAs),
                                          'stdv':np.std(bDAs),
                                          '5%CI':st.t.interval(0.95,len(bDAs)-1, loc=np.mean(bDAs),scale=st.sem(bDAs))[0],
                                          '95%CI':st.t.interval(0.95,len(bDAs)-1, loc=np.mean(bDAs),scale=st.sem(bDAs))[1]}
        statsByAbbrev[abbrev]['prevalence'] = {'mean':np.mean(prevAs),
                                             'median':np.median(prevAs),
                                             'stdv':np.std(prevAs),
                                             '5%CI':st.t.interval(0.95,len(prevAs)-1, loc=np.mean(prevAs),scale=st.sem(prevAs))[0],
                                             '95%CI':st.t.interval(0.95,len(prevAs)-1, loc=np.mean(prevAs),scale=st.sem(prevAs))[1]}
        
        if opts.producetimeseries:
            statsByAbbrev[abbrev]['colonizedDaysTS'] = combineTimeSeries(cTAs, runDays)
            statsByAbbrev[abbrev]['bedDaysTS'] = combineTimeSeries(bTAs,runDays)
        
        for k in valuesToGather.keys():
            statsByAbbrev[abbrev][k] = {'mean':np.mean(statDict[k]['value']),
                                        'median':np.median(statDict[k]['value']),
                                        'stdv':np.std(statDict[k]['value']),
                                        '5%CI':st.t.interval(0.95,len(statDict[k]['value'])-1, loc=np.mean(statDict[k]['value']),scale=st.sem(statDict[k]['value']))[0],
                                        '95%CI':st.t.interval(0.95,len(statDict[k]['value'])-1, loc=np.mean(statDict[k]['value']),scale=st.sem(statDict[k]['value']))[1]}
            if opts.producetimeseries:
                statsByAbbrev[abbrev]["{0}TS".format(k)] = combineTimeSeries(statDict[k]['ts'], runDays, 1)
            
        if abbrev in xdroAbbrevs:
            statsByAbbrev[abbrev]["xdroAdmissions"] = {'mean':np.mean(statDict['arrivals']['value']),
                                                           'median':np.median(statDict['arrivals']['value']),
                                                           'stdv':np.std(statDict['arrivals']['value']),
                                                           '5%CI':st.t.interval(0.95, len(statDict['arrivals']['value']) - 1, loc=np.mean(statDict['arrivals']['value']), scale=st.sem(statDict['arrivals']['value']))[0],
                                                           '95%CI':st.t.interval(0.95, len(statDict['arrivals']['value']) - 1, loc=np.mean(statDict['arrivals']['value']), scale=st.sem(statDict['arrivals']['value']))[1]}
            if opts.producetimeseries:
                statsByAbbrev[abbrev]['xdroAdmissionsTS'] = combineTimeSeries(statDict['arrivals']['ts'],runDays,1)
        
        else:
            statsByAbbrev[abbrev]["xdroAdmissions"] = {'mean':0.0,
                                                       'median':0.0,
                                                       'stdv':0.0,
                                                       '5%CI':0.0,
                                                       '95%CI':0.0}
            if opts.producetimeseries:
                statsByAbbrev[abbrev]['xdroAdmissionsTS'] = {'mean':[0.0 for x in range(0,runDays)],
                                                            'median':[0.0 for x in range(0,runDays)],
                                                            'stdv':[0.0 for x in range(0,runDays)],
                                                            '5%CI':[0.0 for x in range(0,runDays)],
                                                            '95%CI':[0.0 for x in range(0,runDays)]}
    
    
        
    print "time 1 = {0}".format(time1Counter)
    print "time 2 = {0}".format(time2Counter)
    print "Writing Files"
    sys.stdout.flush()        
    with open("{0}_stats_by_tier.csv".format(outFileName),"wb") as f:
        csvWriter = csv.writer(f)
        headingRow = ['Facility Abbrev','Tier Of Care','Colonized Patient Days','Paitent Bed Days','Prevalence']
        for k in valuesToGatherList:
            headingRow.append(tableHeadings[k])
        if len(xdroAbbrevs) > 0:
            headingRow.append("XDRO Admissions")
         
        csvWriter.writerow(headingRow)
        #csvWriter.writerow(['Facility Abbrev','Tier Of Care','Colonized Patient Days',','New Colonizations','CRE Colonized Admissions','Admissions'])
        for abbrev,tD in statsByTier.items():
            for tier,d in tD.items():
                entryRow = [abbrev,tier,
                            d['colonizedDays']['mean'],
                            d['bedDays']['mean'],
                            d['prevalence']['mean']]
                for k in valuesToGatherList:
                    entryRow.append(d[k]['mean'])
                
                if len(xdroAbbrevs) > 0:
                    entryRow.append(d['xdroAdmissions']['mean'])
                    
                csvWriter.writerow(entryRow)
                
    with open("{0}_stats_by_abbrev.csv".format(outFileName),"wb") as f:
        csvWriter = csv.writer(f)
        headingRow = ['Facility Abbrev','Colonized Patient Days','Paitent Bed Days','Prevalence']
        for k in valuesToGatherList:
            headingRow.append(tableHeadings[k])
        if len(xdroAbbrevs) > 0:
            headingRow.append("XDRO Admissions")
            
        csvWriter.writerow(headingRow)
        for abbrev,d in statsByAbbrev.items():
            entryRow = [abbrev,
                        d['colonizedDays']['mean'],
                        d['bedDays']['mean'],
                        d['prevalence']['mean']]
            for k in valuesToGatherList:
                entryRow.append(d[k]['mean'])
            if len(xdroAbbrevs) > 0:
                entryRow.append(d['xdroAdmissions']['mean'])
            csvWriter.writerow(entryRow)

    with open("{0}_stats_intervals_by_abbrev.csv".format(outFileName),"wb") as f:
        csvWriter = csv.writer(f)
        csvWriter.writerow(['Facility Abbrev','Prevalence Mean','Prevalence Median',
                            'Prevalence St. Dev.','Prevalence 5% CI','Prevalence 95% CI'])
        for abbrev,d in statsByAbbrev.items():
            csvWriter.writerow([abbrev,
                                d['prevalence']['mean'],
                                d['prevalence']['median'],
                                d['prevalence']['stdv'],
                                d['prevalence']['5%CI'],
                                d['prevalence']['95%CI']])

    

    with open("{0}_stats_intervals_by_abbrev.csv".format(outFileName),"wb") as f:
        csvWriter = csv.writer(f)
        csvWriter.writerow(['Facility Abbrev','Prevalence Mean','Prevalence Median',
                            'Prevalence St. Dev.','Prevalence 5% CI','Prevalence 95% CI'])
        for abbrev,d in statsByAbbrev.items():
            csvWriter.writerow([abbrev,
                                d['prevalence']['mean'],
                                d['prevalence']['median'],
                                d['prevalence']['stdv'],
                                d['prevalence']['5%CI'],
                                d['prevalence']['95%CI']])
            
    if opts.producetimeseries:
        with open("{0}_prevalence_per_day_by_abbrev.csv".format(outFileName),"wb") as f:
            csvWriter = csv.writer(f)
            headRow = ['Day']
            abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
            for abbrev in abbrevsSorted:
                headRow.append("{0}".format(abbrev))
                
            csvWriter.writerow(headRow)
            
            for i in range(0,runDays):
                entryRow = ['{0}'.format(i)]
                for abbrev in abbrevsSorted:
                    dayPrev = statsByAbbrev[abbrev]['colonizedDaysTS']['mean'][i]/statsByAbbrev[abbrev]['bedDaysTS']['mean'][i]
                    entryRow.append('{0}'.format(dayPrev))
                
                csvWriter.writerow(entryRow)
        
        with open("{0}_prevalence_and_incidence_per_day_13mile.csv".format(outFileName),"wb") as f:
            csvWriter = csv.writer(f)
            headRow = ['Day','Prev within 13','Prev outside 13','Prev within Cook','Prev outside Cook',
                       'Prev target','Prev nonTarget','Prev regionWide',
                       'Inc within 13','Inc outside 13','Inc within Cook','Inc outside Cook','Inc target','Inc nonTarget','Inc regionWide']
            csvWriter.writerow(headRow)
            
            for i in range(0,runDays):
                colWithin = sum([statsByAbbrev[x]['colonizedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if x in facilitiesWithin13Miles])
                bedWithin = sum([statsByAbbrev[x]['bedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if x in facilitiesWithin13Miles])
                ncolsWithin = sum([statsByAbbrev[x]['newColonizedTS']['mean'][i] for x in statsByAbbrev.keys() if x in facilitiesWithin13Miles])
                
                colWithout = sum([statsByAbbrev[x]['colonizedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if x not in facilitiesWithin13Miles])
                bedWithout = sum([statsByAbbrev[x]['bedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if x not in facilitiesWithin13Miles])
                ncolsWithout = sum([statsByAbbrev[x]['newColonizedTS']['mean'][i] for x in statsByAbbrev.keys() if x not in facilitiesWithin13Miles])
                
                colWithinC = sum([statsByAbbrev[x]['colonizedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if x in facilitiesWithinCookCounty])
                bedWithinC = sum([statsByAbbrev[x]['bedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if x in facilitiesWithinCookCounty])
                ncolsWithinC = sum([statsByAbbrev[x]['newColonizedTS']['mean'][i] for x in statsByAbbrev.keys() if x in facilitiesWithinCookCounty])

                colWithoutC = sum([statsByAbbrev[x]['colonizedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if x not in facilitiesWithinCookCounty])
                bedWithoutC = sum([statsByAbbrev[x]['bedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if x not in facilitiesWithinCookCounty])
                ncolsWithoutC = sum([statsByAbbrev[x]['newColonizedTS']['mean'][i] for x in statsByAbbrev.keys() if x not in facilitiesWithinCookCounty])

                colTarget = sum([statsByAbbrev[x]['colonizedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if statsByAbbrev[x]['creBundlesHandedOut']['mean'] != 0 or statsByAbbrev[x]['xdroAdmissions']['mean'] != 0])
                bedTarget = sum([statsByAbbrev[x]['bedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if statsByAbbrev[x]['creBundlesHandedOut']['mean'] != 0 or statsByAbbrev[x]['xdroAdmissions']['mean'] != 0])
                ncolsTarget = sum([statsByAbbrev[x]['newColonizedTS']['mean'][i] for x in statsByAbbrev.keys() if statsByAbbrev[x]['creBundlesHandedOut']['mean'] != 0 or statsByAbbrev[x]['xdroAdmissions']['mean'] != 0])
                
                colNonTarget = sum([statsByAbbrev[x]['colonizedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if statsByAbbrev[x]['creBundlesHandedOut']['mean'] == 0 and statsByAbbrev[x]['xdroAdmissions']['mean'] == 0])
                bedNonTarget = sum([statsByAbbrev[x]['bedDaysTS']['mean'][i] for x in statsByAbbrev.keys() if statsByAbbrev[x]['creBundlesHandedOut']['mean'] == 0 and statsByAbbrev[x]['xdroAdmissions']['mean'] == 0])
                ncolsNonTarget = sum([statsByAbbrev[x]['newColonizedTS']['mean'][i] for x in statsByAbbrev.keys() if statsByAbbrev[x]['creBundlesHandedOut']['mean'] == 0 and statsByAbbrev[x]['xdroAdmissions']['mean'] == 0])
                
                colTotal = sum([statsByAbbrev[x]['colonizedDaysTS']['mean'][i] for x in statsByAbbrev.keys()])
                bedTotal = sum([statsByAbbrev[x]['bedDaysTS']['mean'][i] for x in statsByAbbrev.keys()])
                ncolsTotal = sum([statsByAbbrev[x]['newColonizedTS']['mean'][i] for x in statsByAbbrev.keys()])
                
                prevWithin = 0.0
                if bedWithin > 0.0:
                    prevWithin = colWithin/bedWithin
                
                prevWithout = 0.0
                if bedWithin > 0.0:
                    prevWithout = colWithout/bedWithout
                
                prevWithinC = 0.0
                if bedWithinC > 0.0:
                    prevWithinC = colWithinC/bedWithinC

                prevWithoutC = 0.0
                if bedWithoutC > 0.0:
                    prevWithoutC = colWithoutC/bedWithoutC

                prevTarget = 0.0
                if bedTarget > 0.0:
                    prevTarget = colTarget/bedTarget
                    
                prevNonTarget = 0.0
                if bedNonTarget > 0.0:
                    prevNonTarget = colNonTarget/bedNonTarget
                    
                prevTotal = 0.0
                if bedTotal > 0.0:
                    prevTotal = colTotal/bedTotal
                    
                entryRow = ['{0}'.format(i),
                            prevWithin,
                            prevWithout,
                            prevWithinC,
                            prevWithoutC,
                            prevTarget,
                            prevNonTarget,
                            prevTotal,
                            ncolsWithin,
                            ncolsWithout,
                            ncolsWithinC,
                            ncolsWithoutC,
                            ncolsTarget,
                            ncolsNonTarget,
                            ncolsTotal
                            ]
                csvWriter.writerow(entryRow)
                
#         with open("{0}_prevalence_and_incidence_per_day_target.csv".format(outFileName),"wb") as f:
#             csvWriter = csv.writer(f)
#             headRow = ['Day',,'Inc target','Inc nonTarget','Prev region']
#             csvWriter.writerow(headRow)  
            
                   
        with open("{0}_colonized_patients_per_day_by_abbrev.csv".format(outFileName),"wb") as f:
            csvWriter = csv.writer(f)
            headRow = ['Day']
            abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
            for abbrev in abbrevsSorted:
                headRow.append("{0}".format(abbrev))
                #print "{0}: {1}".format(abbrev,statsByAbbrev[abbrev]['colonizedDaysTS']['mean'])
            
            csvWriter.writerow(headRow)
            
            for i in range(0,runDays):
                entryRow = ["{0}".format(i)]
                for abbrev in abbrevsSorted:
                    entryRow.append("{0}".format(statsByAbbrev[abbrev]['colonizedDaysTS']['mean'][i]))
                csvWriter.writerow(entryRow)
                
        with open("{0}_bedsfilled_per_day_by_abbrev.csv".format(outFileName),"wb") as f:
            csvWriter = csv.writer(f)
            headRow = ['Day']
            abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
            for abbrev in abbrevsSorted:
                headRow.append("{0}".format(abbrev))
                #print "{0}: {1}".format(abbrev,statsByAbbrev[abbrev]['bedDaysTS']['mean'])
            
            csvWriter.writerow(headRow)
            
            for i in range(0,runDays):
                entryRow = ["{0}".format(i)]
                for abbrev in abbrevsSorted:
                    entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['mean'][i]))
                csvWriter.writerow(entryRow)            
        
        with open("{0}_bedsfilled_per_day_by_abbrev_full_stats.csv".format(outFileName),"wb") as f:
            csvWriter = csv.writer(f)
            headRow = ['Day']
            abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
            for abbrev in abbrevsSorted:
                headRow.append("{0}_mean".format(abbrev))
                headRow.append("{0}_median".format(abbrev))
                headRow.append("{0}_stdev".format(abbrev))
                headRow.append("{0}_5%CI".format(abbrev))
                headRow.append("{0}_95%CI".format(abbrev))
                
                #print "{0}: {1}".format(abbrev,statsByAbbrev[abbrev]['colonizedDaysTS']['mean'])
            
            csvWriter.writerow(headRow)
            
            for i in range(0,runDays):
                entryRow = ["{0}".format(i)]
                for abbrev in abbrevsSorted:
                    entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['mean'][i]))
                    entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['median'][i]))
                    entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['stdv'][i]))
                    entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['5%CI'][i]))
                    entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['95%CI'][i]))
                
                csvWriter.writerow(entryRow)
        
        
        
        valuesToWrite =[x for x in valuesToGatherList]
        if len(xdroAbbrevs) > 0:
            valuesToWrite.append("xdroAdmissions")
        for key in valuesToWrite:
                with open("{0}_{1}_per_day_by_abbrev.csv".format(outFileName,key),"wb") as f:
                    csvWriter = csv.writer(f)
                    headRow = ['Day']
                    abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
                    for abbrev in abbrevsSorted:
                        headRow.append("{0}".format(abbrev))
                        #print "{0}: {1}".format(abbrev,statsByAbbrev[abbrev]['bedDaysTS']['mean'])
                    
                    csvWriter.writerow(headRow)
                    
                    for i in range(0,runDays):
                        entryRow = ["{0}".format(i)]
                        for abbrev in abbrevsSorted:
                            entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['mean'][i]))
                        csvWriter.writerow(entryRow)            
                
                with open("{0}_{1}_per_day_by_abbrev_full_stats.csv".format(outFileName,key),"wb") as f:
                    csvWriter = csv.writer(f)
                    headRow = ['Day']
                    abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
                    for abbrev in abbrevsSorted:
                        headRow.append("{0}_mean".format(abbrev))
                        headRow.append("{0}_median".format(abbrev))
                        headRow.append("{0}_stdev".format(abbrev))
                        headRow.append("{0}_5%CI".format(abbrev))
                        headRow.append("{0}_95%CI".format(abbrev))
                        
                        #print "{0}: {1}".format(abbrev,statsByAbbrev[abbrev]['colonizedDaysTS']['mean'])
                    
                    csvWriter.writerow(headRow)
                    
                    for i in range(0,runDays):
                        entryRow = ["{0}".format(i)]
                        for abbrev in abbrevsSorted:
                            entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['mean'][i]))
                            entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['median'][i]))
                            entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['stdv'][i]))
                            entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['5%CI'][i]))
                            entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['95%CI'][i]))
                        
                        csvWriter.writerow(entryRow)
    
    
                
                      
    with open("{0}_prev_by_cat.csv".format(outFileName),"wb") as f:
        csvWriter = csv.writer(f)
        headingRow = ['Facility Type','Prevalence Mean']
        for k in valuesToGatherList:
            headingRow.append(tableHeadings[k])
        csvWriter.writerow(headingRow)
        typeDict = {}
        for abbrev,tD in statsByTier.items():
            for tier,d in tD.items():
                if tier not in typeDict.keys():
                    typeDict[tier] = {'prev':[]}
                    for k in valuesToGatherList:
                        typeDict[tier][k] = []
                        
                if not np.isnan(d['prevalence']['mean']):
                    #print d['prevalence']['mean']
                    typeDict[tier]['prev'].append(d['prevalence']['mean'])
                    for k in valuesToGatherList:
                        typeDict[tier][k].append(d[k]['mean'])
            
        for tier,p in typeDict.items():
            entryRow = [tier,np.mean(typeDict[tier]['prev'])]
            for k in valuesToGatherList:
                entryRow.append(np.mean(typeDict[tier][k]))
            
            csvWriter.writerow(entryRow)

if __name__ == "__main__":
    main()


