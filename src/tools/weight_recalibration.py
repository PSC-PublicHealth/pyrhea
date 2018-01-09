############################################################################################
# file:
#	weight_recalibration.py
# Purpose:
#	This program helps with evaluating and calibrating the transfer weights in pyrhea.
#	The program runs a batch of pyrhea simulations. The number of runs in a batch
# 	is set with the option --batch. After a batch is completed, time samples are collect
# 	on the full batch. Those batch time samples are used to calculate the population
# 	ratios per facility.
# Future Work:
#	This program is intended to be extended to run multiple 'trials', where a trial
#	consists of a batch run, analysis, then followed by a heuristic step that determines
#	how to update the transfer weights. This heuristic step has yet to be implemented.
# Notes:
#	Matplotlib will likely complain about the pyplot backend being set multiple times.
#	The backend needs set in this file so it can save graphs to the filesystem.
############################################################################################

import matplotlib
matplotlib.use('Agg') # Hack for no display w/ pyplot. Change pyplot backend *before* loading it.
import matplotlib.pyplot as plt # use plt.savefig(<path>) rather than plt.show()

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

# def loadWeightDict(weightDir, suffix):
#     weightSetDict = {}
#     for key, path in [('vsnf', 'com_to_vsnf_inv_sep_{0}.yaml'.format(suffix)),
#                      ('snf', 'com_to_snf_inv_sep_{0}.yaml'.format(suffix)),
#                      ('hospital', 'com_to_hosp_inv_sep_{0}.yaml'.format(suffix)),
#                      ('ltach', 'com_to_ltach_inv_sep_{0}.yaml'.format(suffix))]:
#         with open(os.path.join(modelDir, path), 'rU') as f:
#             weightSetDict[key] = yaml.load(f)
#     return weightSetDict

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
    #print 'loaded %s' % fname
    return timeSamplesJSON

def _truncateVec(vec):
    return np.maximum(np.minimum(vec, 2.0), 0.0)

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
    #plt.show()
    plt.savefig(title)
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
    if proto_file is not None:
        sys.argv.append("--proto")
        sys.argv.append(proto_file)

    print(repr(sys.argv))
    collect_time_samples.main()
    
    # restore previous arg list
    sys.argv = stored_args

def batch_time_sample_collection(run_desc, proto_file, outfile_name, batch_size, trial_num):
    # set run description and --notes argv options
    stored_args = sys.argv
    sys.argv = []
    sys.argv.append("collect_time_samples.py")
    sys.argv.append(run_desc)
    for i in range(batch_size):
        sys.argv.append("--notes")
        sys.argv.append("notes_t{}_r{}.pkl".format(trial_num,i)) # assumes nomenclature from main()
    sys.argv.append("--out")
    sys.argv.append(outfile_name)
    if proto_file is not None:
        sys.argv.append("--proto")
        sys.argv.append(proto_file)

    print(repr(sys.argv))
    collect_time_samples.main()
    
    # restore previous arg list
    sys.argv = stored_args


def main():
    parser = OptionParser(usage="""
    %prog run_description.yaml [--proto time_samples_prototype.yaml] [--batch batch_size]
    """)
    parser.add_option('--proto', action='store', type='string',
			help="Prototype YAML file to provide sampling times")
    parser.add_option('--batch', action='store', type='int',
			help="Number of simulations to run for one batch")
    opts, args = parser.parse_args()
    if len(args) < 1:
        parser.error('A YAML run description is required')
    parser.destroy()

    if opts.batch:
        batch_size = opts.batch
    else:
        batch_size = 4

    num_trials = 1 # Currently set to 1 because no heursitic step is implemented

    # check prototype file schema (so problems are known before simulations start)
    if opts.proto:
        checkInputFileSchema(opts.proto, os.path.join(pyrhea.SCHEMA_DIR, collect_time_samples.OUTPUT_SCHEMA))
    else:
        opts.proto = None

    # Init (Build FacDict)
    inputDict = checkInputFileSchema(args[0], os.path.join(pyrhea.SCHEMA_DIR, pyrhea.INPUT_SCHEMA))
    modelDir = inputDict['modelDir']

    pyrheautils.PATH_STRING_MAP['MODELDIR'] = modelDir
    implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])
    pyrheautils.PATH_STRING_MAP['IMPLDIR'] = implDir

    facDirList = [pyrheautils.pathTranslate(pth) for pth in inputDict['facilityDirs']]
    facDict = readFacFiles(facDirList)

    # get initial tranfer weights
    #prevWeightDict = loadInitialWeightDict(modelDir)
    #writeSepFiles(prevWeightDict, 'init', './weights') # store initial weights in 'weights' dir
    #writeWeightsToModel(prevWeightDict, modelDir) # write to modelDir

    print("Begin.")
    for t in range(num_trials):
        print("> trial {0}".format(t))
        for r in range(batch_size): #TODO: handle signals SIGINT SIGTSTP & cleanly exit
            print("> batch run {0}".format(r))
    
            ## set new weights
    
            ## call pyrhea; provide run description and output file
            run_pyrhea(args[0], "notes_t{}_r{}.pkl".format(t,r))
            
            print("< exiting batch run {0} ...".format(r))
            # End batch loop

        ## call collect_time_samples; provide run description and input notes file
        batch_time_sample_collection(args[0], opts.proto, "time_samples_trial_{0}.yaml".format(t), batch_size, t)
        #time_sample_collection(args[0], None, "notes_{0}.pkl".format(i), "time_samples_pass_{0}.yaml".format(i))
        
        ## call calcPopRatios
        timeSamples = loadTimeSamples("time_samples_trial_{0}.yaml".format(t))
        #timeSamples = loadTimeSamples("mytimesamples.yaml")
        popRatios = calcPopRatios(timeSamples, facDict)
        #print("Pop Ratios:")
        #print(repr(popRatios))

        ## call plotPopRatios TODO: store output to PNG
        plotPopRatios(popRatios, facDict, 'meanPop_{}'.format(t))

        ## call rescaleWeights
        #newWeightDict = rescaleWeights(prevWeightDict, popRatios) 
        
        ## call writeSepFiles
        #writeSepFiles(newWeightDict, i, './weights')
        
        ## ref new inv_sep files in xfer_with_draw_constants yaml
        #writeWeightsToModel(newWeightDict, modelDir) # write to modelDir

        print("< exiting trial {0} ...".format(t))
        # End trial loop	

if __name__ == "__main__":
	main()

