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
                    returnDict[abbrev][CareTier.names[tier]] = {'colonizedDays': np.sum(lVec[burninDays:]),
                                                                'bedDays':np.sum(totVecD[tier][burninDays:])}
                    
                    for key in translationDict.keys():
                        returnDict[abbrev][CareTier.names[tier]][key] = 0.0
                        
#                                                                 'newColonized': 0.0,
#                                                                 'creArrivals':0.0,
#                                                                 'arrivals': 0.0,
#                                                                 'contactPrecautionDays':0.0,
#                                                                 'creBundlesHandedOut':0.0 }
    
    for key,noteKey in translationDict.items():
        for abbrev in abbrevList:
            tplList = print_counts.getTimeSeriesList(abbrev,specialDict,noteKey)
            for dayVec,curves in tplList:
                for tpl, curve in curves.items():
                    returnDict[abbrev][CareTier.names[tpl]][key] = np.sum(curve[burninDays:])
    
#     for abbrev in abbrevList:
#         tplList = print_counts.getTimeSeriesList(abbrev,specialDict,'localtiercrearrivals')
#         for dayVec,curves in tplList:
#             for tpl, curve in curves.items():
#                 returnDict[abbrev][CareTier.names[tpl]]['creArrivals'] = np.sum(curve[burninDays:])
#         
#     for abbrev in abbrevList:
#         tplList = print_counts.getTimeSeriesList(abbrev,specialDict,'localtierarrivals')
#         for dayVec,curves in tplList:
#             for tpl, curve in curves.items():
#                 returnDict[abbrev][CareTier.names[tpl]]['arrivals'] = np.sum(curve[burninDays:])
#                 
#     for abbrev in abbrevList:
#         tplList = print_counts.getTimeSeriesList(abbrev,specialDict,'localtierCP')
#         for dayVec,curves in tplList:
#             for tpl, curve in curves.items():
#                 returnDict[abbrev][CareTier.names[tpl]]['contactPrecautionDays'] = np.sum(curve[burninDays:])
#                 
#     for abbrev in abbrevList:
#         tplList = print_counts.getTimeSeriesList(abbrev,specialDict,'localtierCREBundle')
#         for dayVec,curves in tplList:
#             for tpl, curve in curves.items():
#                 returnDict[abbrev][CareTier.names[tpl]]['creBundlesHandedOut'] = np.sum(curve[burninDays:])    
                
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
    
    for v in totalStats:
        #print "K = {0} v={1}".format(k,v)
        totalCounts[totalStats.index(v)] = v

    #print "totalCounts = {0}".format(totalCounts)
    ### By Tier
    ### Each of these should be the same in terms of the abbrevs and tiers, so we can use the first to index the rest
        
    statsByTier = {}
    statsByAbbrev = {}
    for abbrev,tD in totalCounts[0].items():
        if abbrev not in statsByTier.keys():
            statsByTier[abbrev] = {}
            statsByAbbrev[abbrev] = {}
            statDict = {k:np.array([0.0 for x in range(0,len(totalCounts))]) for k in valuesToGather.keys()}
            cDAs = np.array([0.0 for x in range(0,len(totalCounts))])
            bDAs = np.array([0.0 for x in range(0,len(totalCounts))])
            prevAs = np.array([0.0 for x in range(0,len(totalCounts))])
            
        for tier,d in tD.items():
            if tier not in statsByTier[abbrev].keys():
                statsByTier[abbrev][tier] = {}
                statTierDict = {k:np.array([0.0 for x in range(0,len(totalCounts))]) for k in valuesToGather.keys()}
            cDs = np.array([float(x[abbrev][tier]['colonizedDays']) for x in totalCounts])
            bDs = np.array([float(x[abbrev][tier]['bedDays']) for x in totalCounts])
            
            for k in valuesToGather.keys():
                statTierDict[k] = np.array([float(x[abbrev][tier][k]) for x in totalCounts])
            #ncDs = np.array([float(x[abbrev][tier]['newColonized']) for x in totalCounts])
            #caDs = np.array([float(x[abbrev][tier]['creArrivals']) for x in totalCounts])
            #aDs = np.array([float(x[abbrev][tier]['arrivals']) for x in totalCounts])
            prevs = []
            for i in range(0,len(cDs)):
                prevs.append((cDs[i]/bDs[i]))
                cDAs[i] += cDs[i]
                bDAs[i] += bDs[i]
                for k in valuesToGather.keys():
                    statDict[k] += statTierDict[k]
#                 ncAs[i] += ncDs[i]
#                 caDAs[i] += caDs[i]
#                 aDAs[i] += aDs[i]
#                 
                
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
            for k in valuesToGather.keys():
                statsByTier[abbrev][tier][k] = {'mean':np.mean(statTierDict[k]),
                                                       'median':np.median(statTierDict[k]),
                                                       'stdv':np.std(statTierDict[k]),
                                                       '5%CI':st.t.interval(0.95, len(statTierDict[k]) - 1, loc=np.mean(statTierDict[k]), scale=st.sem(statTierDict[k]))[0],
                                                       '95%CI':st.t.interval(0.95, len(statTierDict[k]) - 1, loc=np.mean(statTierDict[k]), scale=st.sem(statTierDict[k]))[1]}
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
        
        for k in valuesToGather.keys():
            statsByAbbrev[abbrev][k] = {'mean':np.mean(statDict[k]),
                                        'median':np.median(statDict[k]),
                                        'stdv':np.std(statDict[k]),
                                        '5%CI':st.t.interval(0.95,len(statDict[k])-1, loc=np.mean(statDict[k]),scale=st.sem(statDict[k]))[0],
                                        '95%CI':st.t.interval(0.95,len(statDict[k])-1, loc=np.mean(statDict[k]),scale=st.sem(statDict[k]))[1]}
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
            
    with open("{0}_stats_by_tier.csv".format(outFileName),"wb") as f:
        csvWriter = csv.writer(f)
        headingRow = ['Facility Abbrev','Tier Of Care','Colonized Patient Days','Paitent Bed Days','Prevalence']
        for k in valuesToGatherList:
            headingRow.append(tableHeadings[k])
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
                
                csvWriter.writerow(entryRow)
                
    with open("{0}_stats_by_abbrev.csv".format(outFileName),"wb") as f:
        csvWriter = csv.writer(f)
        headingRow = ['Facility Abbrev','Colonized Patient Days','Paitent Bed Days','Prevalence']
        for k in valuesToGatherList:
            headingRow.append(tableHeadings[k])
        csvWriter.writerow(headingRow)
        for abbrev,d in statsByAbbrev.items():
            entryRow = [abbrev,
                        d['colonizedDays']['mean'],
                        d['bedDays']['mean'],
                        d['prevalence']['mean']]
            for k in valuesToGatherList:
                entryRow.append(d[k]['mean'])
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
            entryRow = [tier,np.mean(typeDict[tier]['prev'])]
            for k in valuesToGatherList:
                entryRow.append(np.mean(typeDict[tier][k]))
            
            csvWriter.writerow(entryRow)

if __name__ == "__main__":
    main()


