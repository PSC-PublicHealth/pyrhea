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

import sys
import os
import glob
from multiprocessing import Pool
from optparse import OptionParser
from collections import OrderedDict
import types
import time

import numpy as np
import pandas as pd
import scipy.stats as st
import yaml

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))
import pyrheautils
from facilitybase import CareTier, PthStatus
from time_series_plotter import mergeNotesFiles, getTimeSeriesList
from pyrhea import getLoggerConfig
import tools_util as tu
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
        print "finished opening {0}".format(note)
        for key, noteKey in translationDict.items():
            print "  key = {0}".format(key)
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
                    entryDF = pd.concat((entryDF, facDF), sort=True)
            entryDF = entryDF.reset_index(drop=True)
            bigDF = pd.merge(bigDF, entryDF, how='outer', suffixes=['', '_' + key])
    except Exception as e:
        print "for file {0} there is an exception {1}".format(note,e)

    return bigDF


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


def convertToStats(baseGp, fieldsOfInterest):
    """
    This creates a new DataFrame containing _median, _mean, _stdv, _Q1, _Q3, _5%CT, and _95%CT for
    each of the columns in fieldsOfInterest
    """
    dropL = [key for key in ['day', 'run'] if key in baseGp.columns]
    #dropL = []
    statSr = pd.concat([baseGp.mean().drop(dropL).add_suffix('_mean'),
                        baseGp.median().drop(dropL).add_suffix('_median'),
                        baseGp.std().drop(dropL).add_suffix('_stdv'),
                        baseGp.quantile(0.25).drop(dropL).add_suffix('_Q1'),
                        baseGp.quantile(0.75).drop(dropL).add_suffix('_Q3'),
                        computeConfidenceIntervals(baseGp, fieldsOfInterest)
                        ], axis=0)
    return statSr


def generateXDROAdmissions(df, xdroAbbrevs):
    admV = df['arrivals']
    flagV = [(abbrev in xdroAbbrevs) for abbrev in df['abbrev']]
    rsltV = np.where(flagV, admV, 0.0)
    return rsltV


def addBedStats(df):
    cols = PthStatus.names.values()
    return df.assign(colonizedDays=df['COLONIZED'], bedDays=df[cols].sum(axis=1))


def sumDaysWithOffsets(fullDF, valuesToGather, baseDays, fieldsOfInterest):
    """
    Sum relevant values over the duration of the intervention.
    Different attributes have different day offsets, so we need to do this manually;
    we are preserving the functionality of the original gatherCounts.py without really
    understanding why it was implemented in this way.
    baseDays would normally be burninDays + scenarioWaitDays
    """
    rslt = pd.Series()

    dayOffset = valuesToGather['pthStatus'][1]
    cols = [nm for nm in PthStatus.names.values()]
    sr = fullDF[fullDF['day'] >= baseDays + dayOffset][cols].sum()
    rslt['colonizedDays'] = sr['COLONIZED']
    rslt['bedDays'] = sr[cols].sum()
    for key in fieldsOfInterest:
        if key in valuesToGather:
            dayOffset = valuesToGather[key][1]
            rslt[key] = fullDF[fullDF['day'] >= baseDays + dayOffset][key].sum()

    if 'xdroAdmissions' in fieldsOfInterest:
        assert 'xdroAdmissions' in fullDF.columns, 'xdroAdmissions was not added to the df?'
        dayOffset = valuesToGather['arrivals'][1]
        sr = fullDF[fullDF['day'] >= baseDays + dayOffset]['xdroAdmissions'].sum()
        rslt['xdroAdmissions'] = sr
    return rslt


def sumSelectedGroups(inDF, facilitiesWithin13Miles, facilitiesWithinCookCounty, targetFacilitySet):
    fullSr = inDF.sum(numeric_only=True).drop(['index', 'day', 'run'])
    fullSr['prevalence'] = fullSr['colonizedDays'] / fullSr['bedDays']
    fullSr = fullSr.add_suffix('_regn')
    in13mi = inDF['abbrev'].isin(facilitiesWithin13Miles)
    out13mi = np.logical_not(in13mi)
    inCook = inDF['abbrev'].isin(facilitiesWithinCookCounty)
    outCook = np.logical_not(inCook)
    inTarget = inDF['abbrev'].isin(targetFacilitySet)
    outTarget = np.logical_not(inTarget)
    for sel, sfx in [(in13mi, '_in13mi'), (out13mi, '_out13mi'), (inCook, '_inCook'),
                     (outCook, '_outCook'), (inTarget, '_inTarget'), (outTarget, '_outTarget')]:
        thisSr = inDF[sel].sum(numeric_only=True).drop(['index', 'day', 'run'])
        thisSr['prevalence'] = thisSr['colonizedDays'] / thisSr['bedDays']
        fullSr = pd.concat([fullSr, thisSr.add_suffix(sfx)], axis=0)
    return fullSr


def main():
    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl [--out outname.yaml] run_descr.yaml
    """)

    parser.add_option('-n', '--notes', action='append', type='string',
                     help="Notes filename - may be repeated; may be either .mpz or .pkl")
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
    parser.add_option('--producetable', action='store_true', default=False,
                      help='write the entire dataset as a Pandas DataFrame in .mpz format')
    parser.add_option('--nocsv', action='store_true', default=False,
                      help='do not write the output csv files.  Implies --producetable.')
    parser.add_option('--nostats', action='store_true', default=False,
                      help=('suppress calculation of stats (mean, median, etc.) leaving only'
                            'raw csv output'))

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A run file is required.  Try "-h"')
    if opts.nocsv:
        opts.producetable = True

    inputDict = tu.readModelInputs(args[0])
    if opts.out:
        outFileName = opts.out
    else:
        outFileName = DEFAULT_OUT_FILE
    dfOutFileName = os.path.splitext(outFileName)[0] + '_dataframe.mpz'

    if 'trackedFacilities' not in inputDict or not inputDict['trackedFacilities']:
        raise RuntimeError('Run description file does not list any tracked facilities')

    pyrheautils.prepPathTranslations(inputDict)
    facDict = tu.getFacDict(inputDict)

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
                                                  'CRE CP Days due to active surveillance')),
                                  ('xdroCPDays', ('localtierxdroCP', 1,
                                                  'CRE CP Days due to xdro registry')),
                                  ('otherCPDays', ('localtierotherCP', 1,
                                                   'CP Days for other reasons')),
                                  ('pthStatus', ('localtierpathogen', 0, None))
                                  ])

    notes = []
    if opts.glob:
        notes = glob.glob('{0}'.format(opts.notes[0]))
    else:
        notes = opts.notes

    abbrevList = inputDict['trackedFacilities']
    nprocs = opts.nprocs
    print "nprocs = {0}".format(nprocs)
    print "notes = {0}".format(notes)

    xdroAbbrevs = []
    if opts.xdroyamlfile:
        with open("{0}".format(opts.xdroyamlfile),"rb") as f:
            xdroparams = yaml.load(f)
            xdroAbbrevs = xdroparams['locationsImplementingScenario']['locAbbrevList']

    print "XDRO Facilities = {0}".format(xdroAbbrevs)

    pklNotes = [nf for nf in notes if nf.endswith('.pkl')]
    mpzNotes = [nf for nf in notes if (nf.endswith('.mpk') or nf.endswith('.mpz'))]
    assert len(pklNotes) + len(mpzNotes) == len(notes), 'some input notes are in unknown format'

    # pkl files get handled in parallel
    argsList = [(pklNotes[i], abbrevList,
                 {key: tpl[0] for key, tpl in valuesToGather.items()}, burninDays)
                for i in range(0, len(pklNotes))]
    if nprocs > 1:
        pool = Pool(nprocs)
        totalStats = pool.map(pool_helper, argsList)
        pool.close()
    else:
        print 'Warning: using serial processing because nprocs == 1'
        totalStats = []
        for args in argsList:
            totalStats.append(pool_helper(args))
    print 'Finished scanning pickled notes'

    # read any mpz files in line
    if mpzNotes:
        print 'reading pre-parsed notes'
        for mFname in mpzNotes:
            totalStats.append(pd.read_msgpack(mFname))
    print 'Finished scanning all notes'

    tsDFL = []
    for idx, df in enumerate(totalStats):
        df['run'] = idx
        tsDFL.append(df)
    tsDF = pd.concat(tsDFL).reset_index()
    tsDFL = []
    #print tsDF

    if xdroAbbrevs:
        tsDF = tsDF.assign(xdroAdmissions=generateXDROAdmissions(tsDF, xdroAbbrevs))

    #print "totalStats = {0}".format(totalStats)
    if opts.producetable:
        print 'Writing overall sample table'
        tsDF.to_msgpack(dfOutFileName)

    print "Processing Outputs"

    fieldsOfInterest = ['colonizedDays', 'bedDays', 'prevalence']
    for key in valuesToGather.keys():
        if key != 'pthStatus':
            fieldsOfInterest.append(key)
    if xdroAbbrevs:
        fieldsOfInterest.append('xdroAdmissions')
    if opts.nocsv:
        print 'CSV files not written'
    else:
        print "Writing Files"

        # Means by tier.  No 'raw' output for this one, so skip it if opts.nostats
        if not opts.nostats:
            print 'Summarizing by abbrev and tier'
            sumDF = tsDF.groupby(['abbrev', 'tier', 'run']).apply(sumDaysWithOffsets, valuesToGather,
                                                                  burninDays + scenarioWaitDays,
                                                                  fieldsOfInterest)
            sumDF = sumDF.assign(prevalence=sumDF['colonizedDays'].divide(sumDF['bedDays']))
            headingRow = ['Facility Abbrev', 'Tier Of Care', 'Colonized Patient Days',
                          'Patient Bed Days', 'Prevalence']
            entries = ['abbrev', 'tier', 'colonizedDays_mean', 'bedDays_mean', 'prevalence_mean']
            for key, tpl in valuesToGather.items():
                if key != 'pthStatus':
                    headingRow.append(tpl[2])
                    entries.append(key + '_mean')
            if xdroAbbrevs:
                headingRow.append("XDRO Admissions")
                entries.append('xdroAdmissions_mean')
            statDF = sumDF.groupby(['abbrev', 'tier']).apply(convertToStats,
                                                             fieldsOfInterest).reset_index()
            statDF.to_csv("{0}_stats_by_tier.csv".format(outFileName), index=False,
                          columns=entries, header=headingRow)

        # Means by abbrev
        print 'Summarizing by abbrev'
        sumDF = tsDF.groupby(['abbrev', 'run']).apply(sumDaysWithOffsets, valuesToGather,
                                                      burninDays + scenarioWaitDays,
                                                      fieldsOfInterest)
        sumDF = sumDF.assign(prevalence=sumDF['colonizedDays'].divide(sumDF['bedDays']))
        if not opts.nostats:
            headingRow = ['Facility Abbrev','Colonized Patient Days','Patient Bed Days',
                          'Prevalence']
            entries = ['abbrev', 'colonizedDays_mean', 'bedDays_mean', 'prevalence_mean']
            for key, tpl in valuesToGather.items():
                if key != 'pthStatus':
                    headingRow.append(tpl[2])
                    entries.append(key + '_mean')
            if xdroAbbrevs:
                headingRow.append("XDRO Admissions")
                entries.append('xdroAdmissions_mean')
            statDF = sumDF.groupby(['abbrev']).apply(convertToStats, fieldsOfInterest).reset_index()
            statDF.to_csv("{0}_stats_by_abbrev.csv".format(outFileName), index=False,
                          columns=entries, header=headingRow)

        # raw records by abbrev
        headingRow = ['Facility Abbrev','run', 'Colonized Patient Days','Patient Bed Days']
        entries = ['abbrev', 'run', 'colonizedDays', 'bedDays']
        for key, tpl in valuesToGather.items():
            if key != 'pthStatus':
                headingRow.append(tpl[2])
                entries.append(key)
        if xdroAbbrevs:
            headingRow.append('XDRO Admissions')
            entries.append('xdroAdmissions')
        sumDF = sumDF.reset_index()
        sumDF.to_csv("{0}_raw_by_abbrev.csv".format(outFileName), index=False,
                     columns=entries, header=headingRow)

        # Prevalence intervals by abbrev
        if not opts.nostats:
            headingRow = ['Facility Abbrev','Prevalence Mean','Prevalence Median',
                          'Prevalence St. Dev.','Prevalence 5% CI','Prevalence 95% CI']
            entries = ['abbrev', 'prevalence_mean', 'prevalence_median', 'prevalence_stdv',
                       'prevalence_5%CI', 'prevalence_95%CI']
            # the statDF needed here is conveniently the same as above
            statDF.to_csv("{0}_stats_intervals_by_abbrev.csv".format(outFileName), index=False,
                          columns=entries, header=headingRow)

        # Prevalence by tier - naming suggests it was originally by facility category?
        print 'Summarizing by tier'
        sumDF = tsDF.groupby(['tier', 'run']).apply(sumDaysWithOffsets, valuesToGather,
                                                    burninDays + scenarioWaitDays,
                                                    fieldsOfInterest)
        sumDF = sumDF.assign(prevalence=sumDF['colonizedDays'].divide(sumDF['bedDays']))
        if not opts.nostats:
            headingRow = ['Facility Type','Prevalence Mean']
            entries = ['tier', 'prevalence_mean']
            for key, tpl in valuesToGather.items():
                if key != 'pthStatus':
                    headingRow.append(tpl[2])
                    entries.append(key + '_mean')
            statDF = sumDF.groupby(['tier']).apply(convertToStats, fieldsOfInterest).reset_index()
            statDF.to_csv("{0}_prev_by_cat.csv".format(outFileName), index=False,
                          columns=entries, header=headingRow)

        # raw records by tier
        headingRow = ['Tier of Care','run', 'Colonized Patient Days','Patient Bed Days']
        entries = ['tier', 'run', 'colonizedDays', 'bedDays']
        if xdroAbbrevs:
            headingRow.append('XDRO Admissions')
            entries.append('xdroAdmissions')
        sumDF = sumDF.reset_index()
        sumDF.to_csv("{0}_raw_by_tier.csv".format(outFileName), index=False,
                     columns=entries, header=headingRow)

        # Prevalence and incidence per day 13 miles, etc.
        print 'Summarizing by group over days'
        targetAbbrevS = set(xdroAbbrevs)
        df = tsDF[['abbrev', 'creBundlesHandedOut']].groupby('abbrev').sum().reset_index()
        targetAbbrevS = targetAbbrevS.union(df[df['creBundlesHandedOut'] != 0]['abbrev'])
        fullDF = addBedStats(tsDF).groupby(['day', 'run']).apply(sumSelectedGroups,
                                                                 facilitiesWithin13Miles,
                                                                 facilitiesWithinCookCounty,
                                                                 targetAbbrevS)
        fullDF = fullDF.reset_index()
        #fullDF.to_msgpack('test_fulldf.mpz')

        allFlds = []
        for fld in fieldsOfInterest:
            allFlds.extend([fld + sfx for sfx in ['_regn', '_in13mi', '_out13mi', '_inCook',
                                                  '_outCook', '_inTarget', '_outTarget']])

        statHeadingRow = ['Day']
        statEntries = ['day']
        rawHeadingRow = ['Day', 'run']
        rawEntries = ['day', 'run']
        for hd, ent in zip(['Prev within 13','Prev outside 13','Prev within Cook',
                            'Prev outside Cook', 'Prev target','Prev nonTarget',
                            'Prev regionWide',
                            'Inc within 13','Inc outside 13','Inc within Cook',
                            'Inc outside Cook', 'Inc target', 'Inc nonTarget','Inc regionWide'],
                            ['prevalence_in13mi', 'prevalence_out13mi', 'prevalence_inCook',
                             'prevalence_outCook', 'prevalence_inTarget', 'prevalence_outTarget',
                             'prevalence_regn',
                             'newColonized_in13mi', 'newColonized_out13mi', 'newColonized_inCook',
                             'newColonized_outCook', 'newColonized_inTarget',
                             'newColonized_outTarget', 'newColonized_regn']):
            rawHeadingRow.append(hd)
            rawEntries.append(ent)
            for hSfx, eSfx in zip([' mean', ' 5% CI', ' 95% CI'], ['_mean', '_5%CI', '_95%CI']):
                statHeadingRow.append(hd + hSfx)
                statEntries.append(ent + eSfx)

        if not opts.nostats:
            fullStatDF = fullDF.groupby(['day']).apply(convertToStats, allFlds).reset_index()
            #print fullStatDF
            #fullStatDF.to_msgpack('test_fullstatsdf.mpz')
            fullStatDF.to_csv("{0}_prevalence_and_incidence_per_day_13mile.csv".format(outFileName),
                              index=False, columns=statEntries, header=statHeadingRow)
        fullDF.to_csv("{0}_raw_prevalence_and_incidence_pre_day_13mile.csv".format(outFileName),
                      index=False, columns=rawEntries, header=rawHeadingRow)

#     numNotesFiles = len(totalStats)
#     ### By Tier
#     ### Each of these should be the same in terms of the abbrevs and tiers, so we can use the first to index the rest
# 
#     statsByTier = {}
#     statsByAbbrev = {}
#     time1Counter = 0.0
#     time2Counter = 0.0
#     for abbrev,tD in totalCounts[0].items():
#         if abbrev not in statsByTier.keys():
#             statsByTier[abbrev] = {}
#             statsByAbbrev[abbrev] = {}
#             statDict = {k:{'value': np.array([0.0 for x in range(0,numNotesFiles)]),
#                            'ts': [[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]}
#                            for k in valuesToGather.keys()}
#             #statDictTS = {k:{'value': np.array([0.0 for x in range(0,numNotesFiles)]),
#             #               'ts': [[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]} for k in valuesToGather.keys()}
#             #statDictTS = {k:[[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]}
#             cDAs = np.array([0.0 for x in range(0,numNotesFiles)])
#             bDAs = np.array([0.0 for x in range(0,numNotesFiles)])
#             prevAs = np.array([0.0 for x in range(0,numNotesFiles)])
#             
#             
#             if opts.producetimeseries:
#                 cTAs = [[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]
#                 bTAs = [[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]
#           
#         for tier,d in tD.items():
#             time1 = time.time()
#             if tier not in statsByTier[abbrev].keys():
#                 statsByTier[abbrev][tier] = {}
#                 #statsByTierTS[abbrev][tier] = {}
#                 statTierDict = {k:{'value': np.array([0.0 for x in range(0,numNotesFiles)]),
#                                    'ts': [[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]} for k in valuesToGather.keys()}
#                 #statTierDict = {k:np.array([0.0 for x in range(0,numNotesFiles)]) for k in valuesToGather.keys()}
#                 #statTierDictTS = {k:[[0.0 for x in range(0,runDays)] for y in range(0,numNotesFiles)]} 
#             cDList = []
#             bDList = []
#             if opts.producetimeseries:
#                 cTs = []
#                 bTs = []
#             for k in valuesToGather.keys():
#                 statTierDict[k]['value'] = []
#                 if opts.producetimeseries:
#                     statTierDict[k]['ts'] = []
# 
#             for x in totalCounts:
#                 try:
#                     cDList.append(float(x[abbrev][tier]['colonizedDays']))
#                     bDList.append(float(x[abbrev][tier]['bedDays']))
#                     if opts.producetimeseries:
#                         cTs.append(x[abbrev][tier]['colonizedDaysTS'])
#                         bTs.append(x[abbrev][tier]['bedDaysTS'])
# 
#                     for k in valuesToGather.keys():
#                         statTierDict[k]['value'].append(float(x[abbrev][tier][k]))
#                         if opts.producetimeseries:
#                             statTierDict[k]['ts'].append(x[abbrev][tier]["{0}TS".format(k)])
#                         
#                 except:
#                     cDList.append(0.0)
#                     bDList.append(0.0)
#                     if opts.producetimeseries:
#                         cTList.append(0)
#                         bTList.append(0)
#                     for k in valuesToGather.keys():
#                         statTierDict[k]['value'].append(0.0)
#                         if opts.producetimeseries:
#                             statTierDict[k]['ts'].append(0)
#  
#             cDs = np.array(cDList)
#             bDs = np.array(bDList)
#             for k in valuesToGather.keys():
#                 statTierDict[k]['value'] = np.array(statTierDict[k]['value'])
# 
# 
# #            cDs = np.array([float(x[abbrev][tier]['colonizedDays']) for x in totalCounts])
# #            bDs = np.array([float(x[abbrev][tier]['bedDays']) for x in totalCounts])
# #            if opts.producetimeseries:
# #                cTs = [ x[abbrev][tier]['colonizedDaysTS'] for x in totalCounts ]
# #                bTs = [ x[abbrev][tier]['bedDaysTS'] for x in totalCounts ]
# #                #print "bTs = {0}".format(bTs)
# #            
# #            for k in valuesToGather.keys():
# #                #if k == "newColonized":
# #                #    print np.array([float(x[abbrev][tier][k]) for x in totalCounts])
# #                statTierDict[k]['value'] = np.array([float(x[abbrev][tier][k]) for x in totalCounts])
# #                if opts.producetimeseries:
# #                    statTierDict[k]['ts'] = [ x[abbrev][tier]["{0}TS".format(k)] for x in totalCounts ]
# 
#             #ncDs = np.array([float(x[abbrev][tier]['newColonized']) for x in totalCounts])
#             #caDs = np.array([float(x[abbrev][tier]['creArrivals']) for x in totalCounts])
#             #aDs = np.array([float(x[abbrev][tier]['arrivals']) for x in totalCounts])
#             prevs = []
#             for i in range(0,len(cDs)):
#                 prevs.append((cDs[i]/bDs[i]))
#                 #print cTs
#                 if opts.producetimeseries:
#                     for j in range(0,runDays):
#     #                    print j
#     #                    print cTs[i][j]
#                         cTAs[i][j] += cTs[i][j]
#                         bTAs[i][j] += bTs[i][j]
#                 cDAs[i] += cDs[i]
#                 bDAs[i] += bDs[i]
#                 for k in valuesToGather.keys():
#                     statDict[k]['value'][i] += statTierDict[k]['value'][i]
#                     if opts.producetimeseries:
#                         for j in range(0,runDays):
#                             statDict[k]['ts'][i][j] += statTierDict[k]['ts'][i][j]
# #                 ncAs[i] += ncDs[i]
# #                 caDAs[i] += caDs[i]
# #                 aDAs[i] += aDs[i]
# #                 
#             #print "{0} {1} {2}".format(abbrev,tier,cDs)
#             #print statsByTier
# 
#             time2 = time.time()
#             time1Counter += time2-time1
# 
#             time1 = time.time()
# 
# 
#             statsByTier[abbrev][tier]['colonizedDays'] = {'mean':np.mean(cDs),
#                                                           'median':np.median(cDs),
#                                                           'stdv':np.std(cDs),
#                                                           '5%CI':st.t.interval(0.95,len(cDs)-1, loc=np.mean(cDs),scale=st.sem(cDs))[0],
#                                                           '95%CI':st.t.interval(0.95,len(cDs)-1, loc=np.mean(cDs),scale=st.sem(cDs))[1]}
#             
#             statsByTier[abbrev][tier]['bedDays'] = {'mean':np.mean(bDs),
#                                                     'median':np.median(bDs),
#                                                     'stdv':np.std(bDs),
#                                                     '5%CI':st.t.interval(0.95,len(bDs)-1, loc=np.mean(bDs),scale=st.sem(bDs))[0],
#                                                     '95%CI':st.t.interval(0.95,len(bDs)-1, loc=np.mean(bDs),scale=st.sem(bDs))[1]}
#             
#             statsByTier[abbrev][tier]['prevalence'] = {'mean':np.mean(prevs),
#                                                        'median':np.median(prevs),
#                                                        'stdv':np.std(prevs),
#                                                        '5%CI':st.t.interval(0.95,len(prevs)-1, loc=np.mean(prevs),scale=st.sem(prevs))[0],
#                                                        '95%CI':st.t.interval(0.95,len(prevs)-1, loc=np.mean(prevs),scale=st.sem(prevs))[1]}
#             if opts.producetimeseries:
#                 statsByTier[abbrev][tier]['colonizedDaysTS'] = combineTimeSeries(cTs,runDays,1)
#                 statsByTier[abbrev][tier]['bedDaysTS'] = combineTimeSeries(bTs,runDays,1)
#             for k in valuesToGather.keys():
#                 #print "Key 2 = {0}".format(k)
#                 statsByTier[abbrev][tier][k] = {'mean':np.mean(statTierDict[k]['value']),
#                                                        'median':np.median(statTierDict[k]['value']),
#                                                        'stdv':np.std(statTierDict[k]['value']),
#                                                        '5%CI':st.t.interval(0.95, len(statTierDict[k]['value']) - 1, loc=np.mean(statTierDict[k]['value']), scale=st.sem(statTierDict[k]['value']))[0],
#                                                        '95%CI':st.t.interval(0.95, len(statTierDict[k]['value']) - 1, loc=np.mean(statTierDict[k]['value']), scale=st.sem(statTierDict[k]['value']))[1]}
#                 if opts.producetimeseries:
#                     statsByTier[abbrev][tier]["{0}TS".format(k)] = combineTimeSeries(statTierDict[k]['ts'],runDays,1)
#             
#             if abbrev in xdroAbbrevs:
#                 statsByTier[abbrev][tier]["xdroAdmissions"] = {'mean':np.mean(statTierDict['arrivals']['value']),
#                                                                'median':np.median(statTierDict['arrivals']['value']),
#                                                                'stdv':np.std(statTierDict['arrivals']['value']),
#                                                                '5%CI':st.t.interval(0.95, len(statTierDict['arrivals']['value']) - 1, loc=np.mean(statTierDict['arrivals']['value']), scale=st.sem(statTierDict['arrivals']['value']))[0],
#                                                                '95%CI':st.t.interval(0.95, len(statTierDict['arrivals']['value']) - 1, loc=np.mean(statTierDict['arrivals']['value']), scale=st.sem(statTierDict['arrivals']['value']))[1]}
#                 if opts.producetimeseries:
#                     statsByTier[abbrev][tier]['xdroAdmissionsTS'] = combineTimeSeries(statTierDict['arrivals']['ts'],runDays,1)
#             
#             else:
#                 statsByTier[abbrev][tier]["xdroAdmissions"] = {'mean':0.0,
#                                                                'median':0.0,
#                                                                'stdv':0.0,
#                                                                '5%CI':0.0,
#                                                                '95%CI':0.0}
#                 if opts.producetimeseries:
#                     statsByTier[abbrev][tier]['xdroAdmissionsTS'] = {'mean':[0.0 for x in range(0,runDays)],
#                                                                      'median':[0.0 for x in range(0,runDays)],
#                                                                      'stdv':[0.0 for x in range(0,runDays)],
#                                                                      '5%CI':[0.0 for x in range(0,runDays)],
#                                                                      '95%CI':[0.0 for x in range(0,runDays)]}
#             time2 = time.time()
#             time2Counter += time2-time1
#             
#         for i in range(0,len(cDAs)):
#             prevAs[i] = (cDAs[i]/bDAs[i])
#         
#         statsByAbbrev[abbrev]['colonizedDays'] = {'mean':np.mean(cDAs),
#                                                   'median':np.median(cDAs),
#                                                   'stdv':np.std(cDAs),
#                                                   '5%CI':st.t.interval(0.95,len(cDAs)-1, loc=np.mean(cDAs),scale=st.sem(cDAs))[0],
#                                                   '95%CI':st.t.interval(0.95,len(cDAs)-1, loc=np.mean(cDAs),scale=st.sem(cDAs))[1]}
#         
#         statsByAbbrev[abbrev]['bedDays'] = {'mean':np.mean(bDAs),
#                                           'median':np.median(bDAs),
#                                           'stdv':np.std(bDAs),
#                                           '5%CI':st.t.interval(0.95,len(bDAs)-1, loc=np.mean(bDAs),scale=st.sem(bDAs))[0],
#                                           '95%CI':st.t.interval(0.95,len(bDAs)-1, loc=np.mean(bDAs),scale=st.sem(bDAs))[1]}
#         statsByAbbrev[abbrev]['prevalence'] = {'mean':np.mean(prevAs),
#                                              'median':np.median(prevAs),
#                                              'stdv':np.std(prevAs),
#                                              '5%CI':st.t.interval(0.95,len(prevAs)-1, loc=np.mean(prevAs),scale=st.sem(prevAs))[0],
#                                              '95%CI':st.t.interval(0.95,len(prevAs)-1, loc=np.mean(prevAs),scale=st.sem(prevAs))[1]}
#         
#         if opts.producetimeseries:
#             statsByAbbrev[abbrev]['colonizedDaysTS'] = combineTimeSeries(cTAs, runDays)
#             statsByAbbrev[abbrev]['bedDaysTS'] = combineTimeSeries(bTAs,runDays)
#         
#         for k in valuesToGather.keys():
#             statsByAbbrev[abbrev][k] = {'mean':np.mean(statDict[k]['value']),
#                                         'median':np.median(statDict[k]['value']),
#                                         'stdv':np.std(statDict[k]['value']),
#                                         '5%CI':st.t.interval(0.95,len(statDict[k]['value'])-1, loc=np.mean(statDict[k]['value']),scale=st.sem(statDict[k]['value']))[0],
#                                         '95%CI':st.t.interval(0.95,len(statDict[k]['value'])-1, loc=np.mean(statDict[k]['value']),scale=st.sem(statDict[k]['value']))[1]}
#             if opts.producetimeseries:
#                 statsByAbbrev[abbrev]["{0}TS".format(k)] = combineTimeSeries(statDict[k]['ts'], runDays, 1)
#             
#         if abbrev in xdroAbbrevs:
#             statsByAbbrev[abbrev]["xdroAdmissions"] = {'mean':np.mean(statDict['arrivals']['value']),
#                                                            'median':np.median(statDict['arrivals']['value']),
#                                                            'stdv':np.std(statDict['arrivals']['value']),
#                                                            '5%CI':st.t.interval(0.95, len(statDict['arrivals']['value']) - 1, loc=np.mean(statDict['arrivals']['value']), scale=st.sem(statDict['arrivals']['value']))[0],
#                                                            '95%CI':st.t.interval(0.95, len(statDict['arrivals']['value']) - 1, loc=np.mean(statDict['arrivals']['value']), scale=st.sem(statDict['arrivals']['value']))[1]}
#             if opts.producetimeseries:
#                 statsByAbbrev[abbrev]['xdroAdmissionsTS'] = combineTimeSeries(statDict['arrivals']['ts'],runDays,1)
#         
#         else:
#             statsByAbbrev[abbrev]["xdroAdmissions"] = {'mean':0.0,
#                                                        'median':0.0,
#                                                        'stdv':0.0,
#                                                        '5%CI':0.0,
#                                                        '95%CI':0.0}
#             if opts.producetimeseries:
#                 statsByAbbrev[abbrev]['xdroAdmissionsTS'] = {'mean':[0.0 for x in range(0,runDays)],
#                                                             'median':[0.0 for x in range(0,runDays)],
#                                                             'stdv':[0.0 for x in range(0,runDays)],
#                                                             '5%CI':[0.0 for x in range(0,runDays)],
#                                                             '95%CI':[0.0 for x in range(0,runDays)]}
#     
#     
#         
#     print "time 1 = {0}".format(time1Counter)
#     print "time 2 = {0}".format(time2Counter)
#             
#     if opts.producetimeseries:
#         with open("{0}_prevalence_per_day_by_abbrev.csv".format(outFileName),"wb") as f:
#             csvWriter = csv.writer(f)
#             headRow = ['Day']
#             abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
#             for abbrev in abbrevsSorted:
#                 headRow.append("{0}".format(abbrev))
#                 
#             csvWriter.writerow(headRow)
#             
#             for i in range(0,runDays):
#                 entryRow = ['{0}'.format(i)]
#                 for abbrev in abbrevsSorted:
#                     dayPrev = statsByAbbrev[abbrev]['colonizedDaysTS']['mean'][i]/statsByAbbrev[abbrev]['bedDaysTS']['mean'][i]
#                     entryRow.append('{0}'.format(dayPrev))
#                 
#                 csvWriter.writerow(entryRow)
#         
#                 
#                    
#         with open("{0}_colonized_patients_per_day_by_abbrev.csv".format(outFileName),"wb") as f:
#             csvWriter = csv.writer(f)
#             headRow = ['Day']
#             abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
#             for abbrev in abbrevsSorted:
#                 headRow.append("{0}".format(abbrev))
#                 #print "{0}: {1}".format(abbrev,statsByAbbrev[abbrev]['colonizedDaysTS']['mean'])
#             
#             csvWriter.writerow(headRow)
#             
#             for i in range(0,runDays):
#                 entryRow = ["{0}".format(i)]
#                 for abbrev in abbrevsSorted:
#                     entryRow.append("{0}".format(statsByAbbrev[abbrev]['colonizedDaysTS']['mean'][i]))
#                 csvWriter.writerow(entryRow)
#                 
#         with open("{0}_bedsfilled_per_day_by_abbrev.csv".format(outFileName),"wb") as f:
#             csvWriter = csv.writer(f)
#             headRow = ['Day']
#             abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
#             for abbrev in abbrevsSorted:
#                 headRow.append("{0}".format(abbrev))
#                 #print "{0}: {1}".format(abbrev,statsByAbbrev[abbrev]['bedDaysTS']['mean'])
#             
#             csvWriter.writerow(headRow)
#             
#             for i in range(0,runDays):
#                 entryRow = ["{0}".format(i)]
#                 for abbrev in abbrevsSorted:
#                     entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['mean'][i]))
#                 csvWriter.writerow(entryRow)            
#         
#         with open("{0}_bedsfilled_per_day_by_abbrev_full_stats.csv".format(outFileName),"wb") as f:
#             csvWriter = csv.writer(f)
#             headRow = ['Day']
#             abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
#             for abbrev in abbrevsSorted:
#                 headRow.append("{0}_mean".format(abbrev))
#                 headRow.append("{0}_median".format(abbrev))
#                 headRow.append("{0}_stdev".format(abbrev))
#                 headRow.append("{0}_5%CI".format(abbrev))
#                 headRow.append("{0}_95%CI".format(abbrev))
#                 
#                 #print "{0}: {1}".format(abbrev,statsByAbbrev[abbrev]['colonizedDaysTS']['mean'])
#             
#             csvWriter.writerow(headRow)
#             
#             for i in range(0,runDays):
#                 entryRow = ["{0}".format(i)]
#                 for abbrev in abbrevsSorted:
#                     entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['mean'][i]))
#                     entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['median'][i]))
#                     entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['stdv'][i]))
#                     entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['5%CI'][i]))
#                     entryRow.append("{0}".format(statsByAbbrev[abbrev]['bedDaysTS']['95%CI'][i]))
#                 
#                 csvWriter.writerow(entryRow)
#         
#         
#         
#         valuesToWrite =[x for x in valuesToGatherList]
#         if len(xdroAbbrevs) > 0:
#             valuesToWrite.append("xdroAdmissions")
#         for key in valuesToWrite:
#                 with open("{0}_{1}_per_day_by_abbrev.csv".format(outFileName,key),"wb") as f:
#                     csvWriter = csv.writer(f)
#                     headRow = ['Day']
#                     abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
#                     for abbrev in abbrevsSorted:
#                         headRow.append("{0}".format(abbrev))
#                         #print "{0}: {1}".format(abbrev,statsByAbbrev[abbrev]['bedDaysTS']['mean'])
#                     
#                     csvWriter.writerow(headRow)
#                     
#                     for i in range(0,runDays):
#                         entryRow = ["{0}".format(i)]
#                         for abbrev in abbrevsSorted:
#                             entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['mean'][i]))
#                         csvWriter.writerow(entryRow)            
#                 
#                 with open("{0}_{1}_per_day_by_abbrev_full_stats.csv".format(outFileName,key),"wb") as f:
#                     csvWriter = csv.writer(f)
#                     headRow = ['Day']
#                     abbrevsSorted = sorted([x for x in statsByAbbrev.keys()])
#                     for abbrev in abbrevsSorted:
#                         headRow.append("{0}_mean".format(abbrev))
#                         headRow.append("{0}_median".format(abbrev))
#                         headRow.append("{0}_stdev".format(abbrev))
#                         headRow.append("{0}_5%CI".format(abbrev))
#                         headRow.append("{0}_95%CI".format(abbrev))
#                         
#                         #print "{0}: {1}".format(abbrev,statsByAbbrev[abbrev]['colonizedDaysTS']['mean'])
#                     
#                     csvWriter.writerow(headRow)
#                     
#                     for i in range(0,runDays):
#                         entryRow = ["{0}".format(i)]
#                         for abbrev in abbrevsSorted:
#                             entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['mean'][i]))
#                             entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['median'][i]))
#                             entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['stdv'][i]))
#                             entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['5%CI'][i]))
#                             entryRow.append("{0}".format(statsByAbbrev[abbrev]['{0}TS'.format(key)]['95%CI'][i]))
#                         
#                         csvWriter.writerow(entryRow)

if __name__ == "__main__":
    main()


