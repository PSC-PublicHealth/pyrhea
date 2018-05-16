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
from collections import defaultdict, OrderedDict
import types

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

DEFAULT_OUT_FILE = 'counts_output'


def extractCountsFromNotes(note, abbrevList, translationDict, burninDays):
    """Convert the time series contents of a notes file to a Pandas DataFrame"""
    print "note = {0}".format(note)
    ### totalStats is assumed to be a Manager
    bigDF = pd.DataFrame(columns=['day', 'abbrev', 'tier'])
    try:
        specialDict = mergeNotesFiles([note],False)
        print "here {0}".format(note)
        for key, noteKey in translationDict.items():
            print "key = {0}".format(key)
            entryDF = None
            for abbrev in abbrevList:
                facDF = None
                tplList = getTimeSeriesList(abbrev, specialDict, noteKey)
                if not tplList:
                    continue
                for dayVec, curves in tplList:
                    dayIndex = np.where(dayVec==burninDays)[0][0]
                    for innerKey, lVec in curves.items():
                        thisD = {'day': dayVec[dayIndex:]}
                        if isinstance(innerKey, (int, long)):
                            tier = CareTier.names[innerKey]
                            thisD['tier'] = tier
                            thisD[key] = lVec[dayIndex:]
                        elif isinstance(innerKey, types.TupleType) and len(innerKey) == 2:
                            tier, pthStatus = innerKey
                            tier = CareTier.names[tier]
                            pthStatus = PthStatus.names[pthStatus]
                            thisD['tier'] = tier
                            thisD[pthStatus] = lVec[dayIndex:]
                        else:
                            raise RuntimeError('Unknown time series key format {0} for {1} {2}'
                                               .format(innerKey, abbrev, noteKey))
                        thisDF = pd.DataFrame(thisD)
                        if facDF is None:
                            facDF = thisDF
                        else:
                            facDF = pd.merge(facDF, thisDF, how='outer')
                facDF['abbrev'] = abbrev
                if entryDF is None:
                    entryDF = facDF
                else:
                    entryDF = pd.concat((entryDF, facDF))
            entryDF.reset_index(drop=True)
            bigDF = pd.merge(bigDF, entryDF, how='outer', suffixes=['', '_' + key])
    except Exception as e:
        print "for file {0} there is an exception {1}".format(note,e)

    return bigDF



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

    return returnDict


def pool_helper(args):
    return extractCountsFromNotes(*args)


def computeConfidenceIntervals(df, fieldsOfInterest):
    """This is used via pd.Groupby.apply() to provide CIs for data columns"""
    rsltL = []
    rsltIdx = []
    for fld in fieldsOfInterest:
        valV = df[fld]
        lo, hi = st.t.interval(0.95,len(valV)-1, loc=np.mean(valV),scale=st.sem(valV))
        rsltL.extend([lo, hi])
        rsltIdx.extend([fld+'_5%CI', fld+'_95%CI'])
    return pd.Series(rsltL, index=rsltIdx)


def addStatColumns(baseGp, fieldsOfInterest):
    """
    This creates a new DataFrame containing _median, _mean, _stdv, _5%CT, and _95%CT for
    each of the columns in fieldsOfInterest
    """
    statDF = baseGp.mean()
    dropL = [key + '_mean' for key in statDF.columns if key not in fieldsOfInterest]
    statDF = statDF.add_suffix('_mean').reset_index().drop(dropL, axis=1)
    suffix = '_median'
    df = baseGp.median().add_suffix(suffix).reset_index()
    kwargs = {st : df[st] for st in ['%s%s'% (fld, suffix) for fld in fieldsOfInterest]}
    statDF = statDF.assign(**kwargs)
    suffix = '_stdv'
    df = baseGp.std().add_suffix(suffix).reset_index()
    kwargs = {st : df[st] for st in ['%s%s'% (fld, suffix) for fld in fieldsOfInterest]}
    statDF = statDF.assign(**kwargs)
    df = baseGp.apply(computeConfidenceIntervals, fieldsOfInterest).reset_index()
    for suffix in ['_5%CI', '_95%CI']:
        kwargs = {st : df[st] for st in ['%s%s'% (fld, suffix) for fld in fieldsOfInterest]}
        statDF = statDF.assign(**kwargs)
    return statDF


def generateXDROAdmissions(df, xdroAbbrevs):
    admV = df['arrivals']
    flagV = [(abbrev in xdroAbbrevs) for abbrev in df['abbrev']]
    rsltV = np.where(flagV, admV, 0.0)
    return rsltV


def main():
    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl [--out outname.yaml] run_descr.yaml
    """)

    parser.add_option('-n', '--notes', action='append', type='string',
                     help="Notes filename - may be repeated")
    parser.add_option('-o', '--out', action='store', type='string',
                     help="Output file prefix (defaults to %s)" % DEFAULT_OUT_FILE)
    parser.add_option('-g','--glob', action='store_true',
                      help=("Apply filename globbing for notes files."
                            "  (Remember to protect the filename string from the shell!)"))
    parser.add_option('-x','--xdroyamlfile',type='string', default=None,
                      help=("specify the xdro scenario yaml file if there is one,"
                            " if not, XDRO won't be costed"))
    parser.add_option('-m','--nprocs',type='int',default=1,
                      help='number of cpus to run the costing model over')
    parser.add_option('-t','--producetimeseries',action='store_true',default=False)
    parser.add_option('--producetable', action='store_true', default=False,
                      help='write the entire dataset as a Pandas DataFrame in .mpz format')


    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A run file is required.  Try "-h"')

    inputDict = checkInputFileSchema(args[0], os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    if opts.out:
        outFileName = opts.out
    else:
        outFileName = DEFAULT_OUT_FILE
    dfOutFileName = os.path.splitext(outFileName)[0] + '_dataframe.mpz'

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
    scenarioWaitDays = int(inputDict['scenarioWaitDays'])

    ### Translation Dict
    # Format for the tuple value is (notesKey, offsetDays, tableHeading)
    valuesToGather = OrderedDict([('newColonized', ('localtiernewcolonized', 1,
                                                    'Newly Colonized')),
                                  ('creArrivals', ('localtiercrearrivals', 1,
                                                   'CRE Colonized Patients Admissions')),
                                  ('arrivals', ('localtierarrivals', 1,
                                                'Patient Admissions')),
                                  ('contactPrecautionDays', ('localtierCP', 1,
                                                             'Contact Precaution Days')),
                                  ('creBundlesHandedOut', ('localtierCREBundle', 1,
                                                           'CRE Baths Given Out')),
                                  ('creSwabsUsed', ('localtierCRESwabs', 1,
                                                    'CRE Swabs Used')),
                                  ('newPatientsOnCP', ('localtierpatientsOnCP', 1,
                                                       'Number of Patients Put on CP')),
                                  ('passiveCPDays', ('localtierpassiveCP', 1,
                                                     'CRE CP Days due to passive surveillance')),
                                  ('swabCPDays', ('localtierswabCP', 1,
                                                  'CRE CP Days due to acitve surveillance')),
                                  ('xdroCPDays', ('localtierxdroCP', 1,
                                                  'CRE CP Days due to xdro registry')),
                                  ('otherCPDays', ('localtierotherCP', 1,
                                                   'CP Days for other reasons')),
                                  ('pthStatus', ('localtierpathogen', 0, None))
                                  ])

    valuesToGatherList = valuesToGather.keys()

    tableHeadings = {key: tpl[2] for key, tpl in valuesToGather.items()}

    notes = []
    if opts.glob:
        notes = glob.glob('{0}'.format(opts.notes[0]))
    else:
        notes = opts.notes

    abbrevList = inputDict['trackedFacilities']
    nprocs = opts.nprocs
    print "nprocs = {0}".format(nprocs)
    print "notes = {0}".format(notes)

    argsList = [(notes[i], abbrevList,
                 {key: tpl[0] for key, tpl in valuesToGather.items()}, burninDays)
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
    for idx, df in enumerate(totalStats):
        df['run'] = idx
        tsDFL.append(df)
    tsDF = pd.concat(tsDFL).reset_index()
    #print tsDF

    #print "totalStats = {0}".format(totalStats)
    if opts.producetable:
        tsDF.to_msgpack(dfOutFileName)

    print "Processing Outputs"
    sys.stdout.flush()

    dayOffset = valuesToGather['pthStatus'][1]
    tsSumByTier = (tsDF[tsDF['day'] >= burninDays + scenarioWaitDays + dayOffset]
                   .groupby(['abbrev', 'tier', 'run'])
                   .sum().add_suffix('_sum').reset_index())
    cols = [nm + '_sum' for nm in PthStatus.names.values()]
    tierDF = tsSumByTier[['abbrev', 'tier', 'run']].copy()
    tierDF = tierDF.assign(colonizedDays=tsSumByTier['COLONIZED_sum'],
                           bedDays=tsSumByTier[cols].sum(axis=1))
    tierDF = tierDF.assign(prevalence=tierDF['colonizedDays'].divide(tierDF['bedDays']))

    fieldsOfInterest = ['colonizedDays', 'bedDays', 'prevalence']
    for key in valuesToGather.keys():
        if key == 'pthStatus':
            continue  # Did this one as a special case
        if dayOffset != valuesToGather[key][1]:
            # We need to regenerate the intermediate DataFrame
            dayOffset = valuesToGather[key][1]
            tsSumByTier = (tsDF[tsDF['day'] >= burninDays + scenarioWaitDays + dayOffset]
                          .groupby(['abbrev', 'tier', 'run'])
                          .sum().add_suffix('_sum').reset_index())
        tierDF = tierDF.assign(**{key: tsSumByTier[key + '_sum']})
        fieldsOfInterest.append(key)
    if xdroAbbrevs:
        fieldsOfInterest.append('xdroAdmissions')
        tierDF = tierDF.assign(xdroAdmissions=generateXDROAdmissions(tierDF, xdroAbbrevs))

    tierAbbrevStatDF = addStatColumns(tierDF.groupby(['abbrev', 'tier']), fieldsOfInterest)
    abbrevStatDF = addStatColumns(tierDF.groupby(['abbrev']), fieldsOfInterest)
    tierStatDF = addStatColumns(tierDF.groupby(['tier']), fieldsOfInterest)

    print tierAbbrevStatDF.columns
    #print abbrevStatDF.columns
    #print tierStatDF.columns

    print "Writing Files"
    sys.stdout.flush()

    # Means by tier
    headingRow = ['Facility Abbrev', 'Tier Of Care', 'Colonized Patient Days', 'Patient Bed Days',
                  'Prevalence']
    entries = ['abbrev', 'tier', 'colonizedDays_mean', 'bedDays_mean', 'prevalence_mean']
    for key, tpl in valuesToGather.items():
        if key != 'pthStatus':
            headingRow.append(tpl[2])
            entries.append(key + '_mean')
    if xdroAbbrevs:
        headingRow.append("XDRO Admissions")
        entries.append('xdroAdmissions_mean')
    tierAbbrevStatDF.to_csv("{0}_stats_by_tier.csv".format(outFileName), index=False,
                            columns=entries, header=headingRow)

    # Means by abbrev
    headingRow = ['Facility Abbrev','Colonized Patient Days','Patient Bed Days','Prevalence']
    entries = ['abbrev', 'colonizedDays_mean', 'bedDays_mean', 'prevalence_mean']
    for key, tpl in valuesToGather.items():
        if key != 'pthStatus':
            headingRow.append(tpl[2])
            entries.append(key + '_mean')
    if xdroAbbrevs > 0:
        headingRow.append("XDRO Admissions")
        entries.append('xdroAdmissions_mean')
    abbrevStatDF.to_csv("{0}_stats_by_abbrev.csv".format(outFileName), index=False,
                        columns=entries, header=headingRow)

    # Prevalence intervals by abbrev
    headingRow = ['Facility Abbrev','Prevalence Mean','Prevalence Median',
                  'Prevalence St. Dev.','Prevalence 5% CI','Prevalence 95% CI']
    entries = ['abbrev', 'prevalence_mean', 'prevalence_median', 'prevalence_stdv',
               'prevalence_5%CI', 'prevalence_95%CI']
    abbrevStatDF.to_csv("{0}_stats_intervals_by_abbrev.csv".format(outFileName), index=False,
                        columns=entries, header=headingRow)

    # Prevalence by tier - naming suggests it was originally by facility category?
    headingRow = ['Facility Type','Prevalence Mean']
    entries = ['tier', 'prevalence_mean']
    for key, tpl in valuesToGather.items():
        if key != 'pthStatus':
            headingRow.append(tpl[2])
            entries.append(key + '_mean')
    tierStatDF.to_csv("{0}_prev_by_cat.csv".format(outFileName), index=False,
                      columns=entries, header=headingRow)

    sys.exit('done')

    ### By Tier
    ### Each of these should be the same in terms of the abbrevs and tiers, so we can use the first to index the rest

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
#     with open("{0}_stats_by_tier.csv".format(outFileName),"wb") as f:
#         csvWriter = csv.writer(f)
#         headingRow = ['Facility Abbrev','Tier Of Care','Colonized Patient Days','Paitent Bed Days','Prevalence']
#         for k in valuesToGatherList:
#             headingRow.append(tableHeadings[k])
#         if len(xdroAbbrevs) > 0:
#             headingRow.append("XDRO Admissions")
#          
#         csvWriter.writerow(headingRow)
#         #csvWriter.writerow(['Facility Abbrev','Tier Of Care','Colonized Patient Days',','New Colonizations','CRE Colonized Admissions','Admissions'])
#         for abbrev,tD in statsByTier.items():
#             for tier,d in tD.items():
#                 entryRow = [abbrev,tier,
#                             d['colonizedDays']['mean'],
#                             d['bedDays']['mean'],
#                             d['prevalence']['mean']]
#                 for k in valuesToGatherList:
#                     entryRow.append(d[k]['mean'])
#                 
#                 if len(xdroAbbrevs) > 0:
#                     entryRow.append(d['xdroAdmissions']['mean'])
#                     
#                 csvWriter.writerow(entryRow)
                
#     with open("{0}_stats_by_abbrev.csv".format(outFileName),"wb") as f:
#         csvWriter = csv.writer(f)
#         headingRow = ['Facility Abbrev','Colonized Patient Days','Paitent Bed Days','Prevalence']
#         for k in valuesToGatherList:
#             headingRow.append(tableHeadings[k])
#         if len(xdroAbbrevs) > 0:
#             headingRow.append("XDRO Admissions")
#             
#         csvWriter.writerow(headingRow)
#         for abbrev,d in statsByAbbrev.items():
#             entryRow = [abbrev,
#                         d['colonizedDays']['mean'],
#                         d['bedDays']['mean'],
#                         d['prevalence']['mean']]
#             for k in valuesToGatherList:
#                 entryRow.append(d[k]['mean'])
#             if len(xdroAbbrevs) > 0:
#                 entryRow.append(d['xdroAdmissions']['mean'])
#             csvWriter.writerow(entryRow)

#     with open("{0}_stats_intervals_by_abbrev.csv".format(outFileName),"wb") as f:
#         csvWriter = csv.writer(f)
#         csvWriter.writerow(['Facility Abbrev','Prevalence Mean','Prevalence Median',
#                             'Prevalence St. Dev.','Prevalence 5% CI','Prevalence 95% CI'])
#         for abbrev,d in statsByAbbrev.items():
#             csvWriter.writerow([abbrev,
#                                 d['prevalence']['mean'],
#                                 d['prevalence']['median'],
#                                 d['prevalence']['stdv'],
#                                 d['prevalence']['5%CI'],
#                                 d['prevalence']['95%CI']])
# 
#     
# 
#     with open("{0}_stats_intervals_by_abbrev.csv".format(outFileName),"wb") as f:
#         csvWriter = csv.writer(f)
#         csvWriter.writerow(['Facility Abbrev','Prevalence Mean','Prevalence Median',
#                             'Prevalence St. Dev.','Prevalence 5% CI','Prevalence 95% CI'])
#         for abbrev,d in statsByAbbrev.items():
#             csvWriter.writerow([abbrev,
#                                 d['prevalence']['mean'],
#                                 d['prevalence']['median'],
#                                 d['prevalence']['stdv'],
#                                 d['prevalence']['5%CI'],
#                                 d['prevalence']['95%CI']])
            
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
    
    
                
                      
#     with open("{0}_prev_by_cat.csv".format(outFileName),"wb") as f:
#         csvWriter = csv.writer(f)
#         headingRow = ['Facility Type','Prevalence Mean']
#         for k in valuesToGatherList:
#             headingRow.append(tableHeadings[k])
#         csvWriter.writerow(headingRow)
#         typeDict = {}
#         for abbrev,tD in statsByTier.items():
#             for tier,d in tD.items():
#                 if tier not in typeDict.keys():
#                     typeDict[tier] = {'prev':[]}
#                     for k in valuesToGatherList:
#                         typeDict[tier][k] = []
#                         
#                 if not np.isnan(d['prevalence']['mean']):
#                     #print d['prevalence']['mean']
#                     typeDict[tier]['prev'].append(d['prevalence']['mean'])
#                     for k in valuesToGatherList:
#                         typeDict[tier][k].append(d[k]['mean'])
#             
#         for tier,p in typeDict.items():
#             entryRow = [tier,np.mean(typeDict[tier]['prev'])]
#             for k in valuesToGatherList:
#                 entryRow.append(np.mean(typeDict[tier][k]))
#             
#             csvWriter.writerow(entryRow)

if __name__ == "__main__":
    main()


