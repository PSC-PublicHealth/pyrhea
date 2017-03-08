import csv
import sys,os
import glob
import numpy as np
import scipy.stats as st
import yaml
from multiprocessing import Process,Manager,cpu_count,Pool
from optparse import OptionParser

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))
import pyrheautils
import schemautils
import pathogenbase as pth
from facilitybase import CareTier, PthStatus
from notes_plotter import readFacFiles, scanAllFacilities, checkInputFileSchema
from notes_plotter import importNotes
from notes_plotter import SCHEMA_DIR, INPUT_SCHEMA, CARE_TIERS, FAC_TYPE_TO_CATEGORY_MAP
from pyrhea import getLoggerConfig
import print_counts
from collections import defaultdict
#import affinity

#print "Cpu_Count = {0}".format(cpu_count())
#sys.stdout.flush()
#affinity.set_process_affinity_mask(0,2**50-1)

#os.system("taskset -p 0xff %d"%os.getpid())

DEFAULT_OUT_FILE = 'counts_output.yaml'


def extractCountsFromNotes(note, abbrevList, facDict, translationDict, burninDays):
    print "note = {0}".format(note)
    returnDict = {}
    ### totalStats is assumed to be a Manager
    specialDict = print_counts.mergeNotesFiles([note],False)
    for abbrev in abbrevList:
        returnDict[abbrev] = {}
        tplList = print_counts.getTimeSeriesList(abbrev,specialDict,'localtierpathogen')
        tierAxisD = {}
        tAOffset = 0
        for dayVec, curves in tplList:
            dayIndex = np.where(dayVec==burninDays)[0][0]
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
                    returnDict[abbrev][CareTier.names[tier]] = {'colonizedDays': np.sum(lVec[dayIndex:]),
                                                                'bedDays':np.sum(totVecD[tier][dayIndex:]),
                                                                'colonizedDaysTS':lVec[dayIndex:],
                                                                'bedDaysTS':lVec[dayIndex:]}
                    
                    for key in translationDict.keys():
                        returnDict[abbrev][CareTier.names[tier]][key] = 0.0
                        
    for key,noteKey in translationDict.items():
        for abbrev in abbrevList:
            tplList = print_counts.getTimeSeriesList(abbrev,specialDict,noteKey)
            for dayVec,curves in tplList:
                dayIndex = np.where(dayVec==burninDays)[0][0]
                for tpl, curve in curves.items():
                    returnDict[abbrev][CareTier.names[tpl]][key] = np.sum(curve[dayIndex:])
                    returnDict[abbrev][CareTier.names[tpl]]["{0}TS".format(key)] = curve[dayIndex:]

    return returnDict

def combineTimeSeries(tsArray,runDays, tsIncrement=1):
    returnDict = {'mean':[],'median':[],'stdv':[],'5%CI':[],'95%CI':[]}
    
    for i in range(0,runDays):
        #print "array = {0}".format([x[i] for x in tsArray])
        #print "mean = {0}".format(np.mean([x[i] for x in tsArray]))
        returnDict['mean'].append(np.mean([x[i] for x in tsArray]))
        returnDict['median'].append(np.median([x[i] for x in tsArray]))
        returnDict['stdv'].append(np.std([x[i] for x in tsArray]))
        returnDict['5%CI'].append(st.t.interval(0.95, len([x[i] for x in tsArray]) - 1, loc=np.mean([x[i] for x in tsArray]), scale=st.sem([x[i] for x in tsArray]))[0])
        returnDict['95%CI'].append(st.t.interval(0.95, len([x[i] for x in tsArray]) - 1, loc=np.mean([x[i] for x in tsArray]), scale=st.sem([x[i] for x in tsArray]))[1])
        
    
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
    
    inputDict = checkInputFileSchema(args[0], os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    if opts.out:
        outFileName = opts.out
    else:
        outFileName = DEFAULT_OUT_FILE
        
    if 'trackedFacilities' not in inputDict or not len(inputDict['trackedFacilities']):
        raise RuntimeError('Run description file does not list any tracked facilities')

    modelDir = inputDict['modelDir']
    pyrheautils.PATH_STRING_MAP['MODELDIR'] = os.path.abspath(modelDir)
    implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])
    pyrheautils.PATH_STRING_MAP['IMPLDIR'] = implDir

    facDirList = [pyrheautils.pathTranslate(fPth) for fPth in inputDict['facilityDirs']]
    facDict = readFacFiles(facDirList)
    
    burninDays = int(inputDict['burnInDays'])
    print "burninDays = {0}".format(burninDays)
    runDays = int(inputDict['runDurationDays'])
    print "runDays = {0}".format(runDays)
    
    ### Translation Dict
    valuesToGatherList = ['newColonized', 'creArrivals', 'arrivals', 'contactPrecautionDays', 'creBundlesHandedOut']
    valuesToGather = {'newColonized':'localtiernewcolonized',
                      'creArrivals':'localtiercrearrivals',
                      'arrivals':'localtierarrivals',
                      'contactPrecautionDays': 'localtierCP',
                      'creBundlesHandedOut': 'localtierCREBundle'
                      }
    tableHeadings = {'newColonized':'Newly Colonized',
                      'creArrivals':'CRE Colonized Patients Admissions',
                      'arrivals':'Patient Admissions',
                      'contactPrecautionDays': 'Contact Precaution Days',
                      'creBundlesHandedOut': 'CRE Bundles Given Out'
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
    totalCounts = [{} for x in range(0,len(notes))]
    
    argsList = [(notes[i],abbrevList,facDict, valuesToGather, burninDays) for i in range(0,len(notes))]

    xdroAbbrevs = []
    if opts.xdroyamlfile:
        with open("{0}".format(opts.xdroyamlfile),"rb") as f:
            xdroparams = yaml.load(f)
            xdroAbbrevs = xdroparams['locationsImplementingScenario']['locAbbrevList']
    
    print "XDRO Facilities = {0}".format(xdroAbbrevs)    
        
    p = Pool(nprocs)
    totalStats = p.map(pool_helper,argsList)
    p.close()
    
    '''
    for i in range(0,len(notes),nprocs):
        print "i = {0}".format(i)
        manager = Manager()
        totalStats = manager.dict()
        #newColsReturn = manager.dict()
        print "HERE"
        jobs = []
        end = i + nprocs
        if end > len(notes):
            end = len(notes)

        for j in range(i,end):
            print "j = {0}".format(j)
            #p = Process(target=shit,args=())
            p = Process(target=extractCountsFromNotes, args=(notes[j],abbrevList, facDict, j,totalStats))
            jobs.append(p)
            p.start()

        for proc in jobs:
            proc.join()
    '''
    #print totalStats
    
    for i in range(0,len(totalStats)):
        #print "K = {0} v={1}".format(k,v)
        totalCounts[i] = totalStats[i]

    #print "totalCounts = {0}".format(totalCounts)
    ### By Tier
    ### Each of these should be the same in terms of the abbrevs and tiers, so we can use the first to index the rest
    
    print "Processing Outputs"
    sys.stdout.flush()
    statsByTier = {}
    statsByAbbrev = {}
    for abbrev,tD in totalCounts[0].items():
        if abbrev not in statsByTier.keys():
            statsByTier[abbrev] = {}
            statsByAbbrev[abbrev] = {}
            statDict = {k:{'value': np.array([0.0 for x in range(0,len(totalCounts))]),
                           'ts': [[0.0 for x in range(0,runDays)] for y in range(0,len(totalCounts))]} for k in valuesToGather.keys()}
            #statDictTS = {k:{'value': np.array([0.0 for x in range(0,len(totalCounts))]),
            #               'ts': [[0.0 for x in range(0,runDays)] for y in range(0,len(totalCounts))]} for k in valuesToGather.keys()}
            #statDictTS = {k:[[0.0 for x in range(0,runDays)] for y in range(0,len(totalCounts))]}
            cDAs = np.array([0.0 for x in range(0,len(totalCounts))])
            bDAs = np.array([0.0 for x in range(0,len(totalCounts))])
            prevAs = np.array([0.0 for x in range(0,len(totalCounts))])
            
            
            if opts.producetimeseries:
                cTAs = [[0.0 for x in range(0,runDays)] for y in range(0,len(totalCounts))]
                bTAs = [[0.0 for x in range(0,runDays)] for y in range(0,len(totalCounts))]
            
        for tier,d in tD.items():
            if tier not in statsByTier[abbrev].keys():
                statsByTier[abbrev][tier] = {}
                #statsByTierTS[abbrev][tier] = {}
                statTierDict = {k:{'value': np.array([0.0 for x in range(0,len(totalCounts))]),
                                   'ts': [[0.0 for x in range(0,runDays)] for y in range(0,len(totalCounts))]} for k in valuesToGather.keys()}
                #statTierDict = {k:np.array([0.0 for x in range(0,len(totalCounts))]) for k in valuesToGather.keys()}
                #statTierDictTS = {k:[[0.0 for x in range(0,runDays)] for y in range(0,len(totalCounts))]} 
            cDs = np.array([float(x[abbrev][tier]['colonizedDays']) for x in totalCounts])
            bDs = np.array([float(x[abbrev][tier]['bedDays']) for x in totalCounts])
            if opts.producetimeseries:
                cTs = [ x[abbrev][tier]['colonizedDaysTS'] for x in totalCounts ]
                bTs = [ x[abbrev][tier]['bedDaysTS'] for x in totalCounts ]
            
            for k in valuesToGather.keys():
                #if k == "newColonized":
                #    print np.array([float(x[abbrev][tier][k]) for x in totalCounts])
                statTierDict[k]['value'] = np.array([float(x[abbrev][tier][k]) for x in totalCounts])
                if opts.producetimeseries:
                    statTierDict[k]['ts'] = [ x[abbrev][tier]["{0}TS".format(k)] for x in totalCounts ]
                
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
                
            '''
            statsByTier[abbrev][tier]['newColonized'] = {'mean':np.mean(ncDs),
                                                       'median':np.median(ncDs),
                                                       'stdv':np.std(ncDs),
                                                       '5%CI':st.t.interval(0.95, len(ncDs) - 1, loc=np.mean(ncDs), scale=st.sem(ncDs))[0],
                                                       '95%CI':st.t.interval(0.95, len(ncDs) - 1, loc=np.mean(ncDs), scale=st.sem(ncDs))[1]}
            statsByTier[abbrev][tier]['creArrivals'] = {'mean':np.mean(caDs),
                                                       'median':np.median(caDs),
                                                       'stdv':np.std(caDs),
                                                       '5%CI':st.t.interval(0.95, len(caDs) - 1, loc=np.mean(caDs), scale=st.sem(caDs))[0],
                                                       '95%CI':st.t.interval(0.95, len(caDs) - 1, loc=np.mean(caDs), scale=st.sem(caDs))[1]}
            statsByTier[abbrev][tier]['arrivals'] = {'mean':np.mean(aDs),
                                                     'median':np.median(aDs),
                                                     'stdv':np.std(aDs),
                                                     '5%CI':st.t.interval(0.95, len(aDs) - 1, loc=np.mean(aDs), scale=st.sem(aDs))[0],
                                                     '95%CI':st.t.interval(0.95, len(aDs) - 1, loc=np.mean(aDs), scale=st.sem(aDs))[1]}
            '''
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
            
    '''
        statsByAbbrev[abbrev]['newColonized'] = {'mean':np.mean(ncAs),
                                             'median':np.median(ncAs),
                                             'stdv':np.std(ncAs),
                                             '5%CI':st.t.interval(0.95,len(ncAs)-1, loc=np.mean(ncAs),scale=st.sem(ncAs))[0],
                                             '95%CI':st.t.interval(0.95,len(ncAs)-1, loc=np.mean(ncAs),scale=st.sem(ncAs))[1]}
        
        statsByAbbrev[abbrev]['creArrivals'] = {'mean':np.mean(caDAs),
                                             'median':np.median(caDAs),
                                             'stdv':np.std(caDAs),
                                             '5%CI':st.t.interval(0.95,len(caDAs)-1, loc=np.mean(caDAs),scale=st.sem(caDAs))[0],
                                             '95%CI':st.t.interval(0.95,len(caDAs)-1, loc=np.mean(caDAs),scale=st.sem(caDAs))[1]}
        
        statsByAbbrev[abbrev]['arrivals'] = {'mean':np.mean(aDAs),
                                             'median':np.median(aDAs),
                                             'stdv':np.std(aDAs),
                                             '5%CI':st.t.interval(0.95,len(aDAs)-1, loc=np.mean(aDAs),scale=st.sem(aDAs))[0],
                                             '95%CI':st.t.interval(0.95,len(aDAs)-1, loc=np.mean(aDAs),scale=st.sem(aDAs))[1]}
        '''  
    
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
                    typeDict[tier]['prev'].append(d['prevalence']['mean'])
                    for k in valuesToGatherList:
                        typeDict[tier][k].append(d[k]['mean'])
            
        for tier,p in typeDict.items():
            entryRow = [tier,np.sum(typeDict[tier]['prev'])]
            for k in valuesToGatherList:
                entryRow.append(np.mean(typeDict[tier][k]))
            
            csvWriter.writerow(entryRow)

if __name__ == "__main__":
    main()


