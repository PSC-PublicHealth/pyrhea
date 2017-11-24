import sys, os
from optparse import OptionParser
import numpy as np
import yaml
import pyrhea
import pyrheautils
import phacsl.utils.formats.yaml_tools as yaml_tools

cwd = os.getcwd()
sys.path.append(os.path.join(cwd, "../tools"))

from notes_plotter import readFacFiles, checkInputFileSchema
import collect_time_samples

#TODO: In all functions, change all static facilities to dynamically read from run description
def calcPopRatios(timeSamples, facDict):
    popRatios={}
    for ent in timeSamples:
        abbrev = ent['abbrev']
        valV = np.asfarray(ent['samples']['occupancy'])
        popRatios[abbrev] = float(np.mean(valV)) / facDict[abbrev]['meanPop']['value']
    return popRatios

def loadInitialWeightDict(modelDir):
    weightSetDict = {}
    for key, path in [('vsnf', 'com_to_vsnf_inv_sep.yaml'),
                     ('snf', 'com_to_snf_inv_sep.yaml'),
                     ('hospital', 'com_to_hosp_inv_sep.yaml'),
                     ('ltach', 'com_to_ltach_inv_sep.yaml')]:
        with open(os.path.join(modelDir, path), 'rU') as f:
            weightSetDict[key] = yaml.load(f)
    return weightSetDict

def loadWeightDict(weightDir, suffix):
    weightSetDict = {}
    for key, path in [('vsnf', 'com_to_vsnf_inv_sep_{0}.yaml'.format(suffix)),
                     ('snf', 'com_to_snf_inv_sep_{0}.yaml'.format(suffix)),
                     ('hospital', 'com_to_hosp_inv_sep_{0}.yaml'.format(suffix)),
                     ('ltach', 'com_to_ltach_inv_sep_{0}.yaml'.format(suffix))]:
        with open(os.path.join(modelDir, path), 'rU') as f:
            weightSetDict[key] = yaml.load(f)
    return weightSetDict

def rescaleWeights(weightSetDict, popRatios):
    missingSet = set()
    newWeightSetDict = {}
    for key in ['vsnf', 'snf', 'hospital', 'ltach']:
        newWeightSetDict[key] = {}
        for comAbbrev, destD in weightSetDict[key].items():
            newWeightSetDict[key][comAbbrev] = destD.copy()
            for destAbbrev, wt in destD.items():
                if destAbbrev in popRatios:
                    pR = popRatios[destAbbrev]
                    if pR > 10.0:
                        # This is clearly an anomaly, like a facility open only a
                        # small fraction of the year.  Do not correct it.
                        pR = 1.0
                    newWt = wt / pR
                    newWeightSetDict[key][comAbbrev][destAbbrev] = newWt
                elif destAbbrev in missingSet:
                    pass
                else:
                    print '######## missing %s' % destAbbrev
                    missingSet.add(destAbbrev)
        print 'finished %s' % key
    return newWeightSetDict

def writeSepFiles(newWeightSetDict, suffix, writeDir):
    if not os.path.exists(writeDir):
        os.makedirs(writeDir)
    for key, path in [('vsnf', 'com_to_vsnf_inv_sep_%s.yaml' % suffix),
                     ('snf', 'com_to_snf_inv_sep_%s.yaml' % suffix),
                     ('hospital', 'com_to_hosp_inv_sep_%s.yaml' % suffix),
                     ('ltach', 'com_to_ltach_inv_sep_%s.yaml' % suffix)]:
        fullPath = os.path.join(writeDir, path)
        with open(fullPath, 'w') as f:
            yaml.safe_dump(newWeightSetDict[key], f, default_flow_style=True, indent=4,
                           encoding='utf-8', width=130, explicit_start=True)
        print 'wrote %s' % fullPath

def writeWeightsToModel(newWeightSetDict, modelDir):
    for key, path in [('vsnf', 'com_to_vsnf_inv_sep.yaml'),
                     ('snf', 'com_to_snf_inv_sep.yaml'),
                     ('hospital', 'com_to_hosp_inv_sep.yaml'),
                     ('ltach', 'com_to_ltach_inv_sep.yaml')]:
        fullPath = os.path.join(modelDir, path)
        with open(fullPath, 'w') as f:
            yaml.safe_dump(newWeightSetDict[key], f, default_flow_style=True, indent=4,
                           encoding='utf-8', width=130, explicit_start=True)
        print 'wrote %s' % fullPath   

facDict = {}
def loadTimeSamples(fname):
    with open(fname, 'r') as f:
        timeSamplesJSON = yaml.load(f)
    #for ent in timeSamplesJSON:
        #assert ent['time'] == 465
        #assert ent['abbrev'] in facDict
    #print 'loaded %s' % fname
    return timeSamplesJSON

def plotPopRatios(popRatios, facDict, title): 
    plt.xlabel('Actual mean population / meanPop')
    plt.ylabel('counts')
    statD = {}
    vecD = {}
    keyCatPairs = [('hosp', 'HOSPITAL'), ('ltac', 'LTACH'), ('snf', 'SNF'), ('vsnf', 'VSNF')]
    for key, ctg in keyCatPairs:
        vecD[key] = [v for k, v in popRatios.items() if facDict[k]['category'] == ctg]
        mean1 = np.mean(vecD[key])
        stdv1 = np.std(vecD[key], ddof=1)
        outLim = mean1 + 5.0*stdv1
        indexV = (vecD[key] <= outLim)
        nOut = len(vecD[key]) - np.count_nonzero(indexV)
        mean2 = np.mean(np.asarray(vecD[key])[indexV])
        stdv2 = np.std(np.asarray(vecD[key])[indexV], ddof=1)
        statD[key] = (mean2, stdv2, nOut)
    plt.hist([_truncateVec(vecD[key]) for key in ['hosp', 'ltac', 'snf', 'vsnf']],
            50, range=(0.0, 2.0), stacked=True,
            label=['HOSPITAL', 'LTACH', 'SNF', 'VSNF'])
    plt.legend()
    plt.title(title)
    plt.show()
    print 'Samples above max:'
    for k, v in popRatios.items():
        if v > 2.0 or v < 0.0:
            print '  %s: %s' % (k, v)
    print 'Distributions:'
    for key, ctg in keyCatPairs:
        print '  %s: %s +- %s excluding %s samples' % (ctg, statD[key][0], statD[key][1], statD[key][2])

def run_pyrhea(run_desc, outfile_name):
    # prepare sys.argv for pyrhea.main(); needs run_desc and notes output name
    stored_args = sys.argv
    sys.argv = []
    sys.argv.append("pyrhea.py")
    sys.argv.append(run_desc)
    sys.argv.append("--out")
    sys.argv.append(outfile_name)

    print(repr(sys.argv))
    pyrhea.main()

    # restore previous arg list
    sys.argv = stored_args

def time_sample_collection(run_desc, proto_file, notes_file, outfile_name):
    # set run description and --notes argv options
    stored_args = sys.argv
    sys.argv = []
    sys.argv.append("collect_time_samples.py")
    sys.argv.append(run_desc)
    sys.argv.append("--notes")
    sys.argv.append(notes_file)
    sys.argv.append("--out")
    sys.argv.append(outfile_name)
    sys.argv.append("--proto")
    sys.argv.append(proto_file)

    print(repr(sys.argv))
    collect_time_samples.main()
    
    # restore previous arg list
    sys.argv = stored_args

def main():
    parser = OptionParser(usage="""
    %prog run_description.yaml
    """)
    opts, args = parser.parse_args()
    if len(args) != 2:
        parser.error('A YAML run description & time_samples prototype is required')
    parser.destroy()

    # Init (Build FacDict) TODO: Change from static path names to dynamically read from run description
    runDesc = checkInputFileSchema(args[0], os.path.join(pyrhea.SCHEMA_DIR, pyrhea.INPUT_SCHEMA))
    modelDir = runDesc['modelDir']
    allKeySet, recs = yaml_tools.parse_all(os.path.join(modelDir, 'facilityfacts')) #facilityfactsCurrent2013
    commKeySet, commRecs = yaml_tools.parse_all(os.path.join(modelDir, 'synthCommunities'))
    facDict = {r['abbrev']:r for r in recs}
    for r in commRecs:
        facDict[r['abbrev']] = r

    prevWeightDict = loadInitialWeightDict(modelDir)
    writeSepFiles(prevWeightDict, 'init', './weights') # store initial weights in 'weights' dir
    writeWeightsToModel(prevWeightDict, modelDir) # write to modelDir

    print("Begin.")
    for i in range(5): #TODO: Make loop count cli arg; handle signals SIGINT SIGTSTP & cleanly exit
        print("> run {0}".format(i))

	## set new weights

        ## call pyrhea; provide run description and output file
        run_pyrhea(args[0], "notes_{0}.pkl".format(i))
        
        ## call collect_time_samples; provide run description and input notes file
        time_sample_collection(args[0], args[1], "notes_{0}.pkl".format(i), "time_samples_pass_{0}.yaml".format(i))
        
        ## call calcPopRatios
	timeSamples = loadTimeSamples("time_samples_pass_{0}.yaml".format(i))
        popRatios = calcPopRatios(timeSamples, facDict)
        print("Pop Ratios:")
        print(repr(popRatios))

        ## call plotPopRatios TODO: store output to PNG
        
        ## call rescaleWeights
        newWeightDict = rescaleWeights(prevWeightDict, popRatios) 
        
        ## call writeSepFiles
        writeSepFiles(newWeightDict, i, './weights')
        
        ## ref new inv_sep files in xfer_with_draw_constants yaml
        writeWeightsToModel(newWeightDict, modelDir) # write to modelDir
	
	print("< exiting run {0} ...".format(i))

if __name__ == "__main__":
	main()

