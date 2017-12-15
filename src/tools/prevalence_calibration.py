#! /usr/bin/env python

###################################################################################
# Copyright   2015, Pittsburgh Supercomputing Center (PSC).  All Rights Reserved. #
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

# file:
#	prevalence_calibration.py
# Purpose:
#	This program helps with calibrating prevalence. The program takes a list
#	of notes files as input and produces a matrix of per-facility prevalence
#	values. It also reads in target tau values from cre_constants (hardcoded)
#	to build a vector of target prevalence values. The log liklihood of that
#	target vector, given the prevalence matrix, is calulated and printed out.
# Future Work:
#	This program is intended to be extended to run multiple 'trials', where a
#	trial consists of a bach run, calculating the log liklihood, then followed
#	by a heuristic step that determines how to update the tau values to
#	maximize the log liklihood.


import os
import os.path
import sys
import logging
import logging.config
import yaml
from optparse import OptionParser
import types
import glob
from collections import defaultdict
if 'line_profiler' not in dir():
    def profile(func):
        """
        Dummy substitute for __builtin__.profile line profiler
        """
        return func

import numpy as np
from scipy.stats import norm

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))

import pyrheautils
import schemautils
from facilitybase import CareTier
from notes_plotter import readFacFiles, checkInputFileSchema
from notes_plotter import importNotes
from notes_plotter import SCHEMA_DIR, INPUT_SCHEMA
from pyrhea import getLoggerConfig

logger = None

MIN_STDV = 1.0
SHRINK = 0.99999

def buildEnumNameStr(name, ind):
    return '%s_%d' % (name, ind)

def splitEnumNameStr(enumName):
    offset = enumName.rfind('_')
    patchName = enumName[:offset]
    ind = int(enumName[offset+1:])
    return patchName, ind

def buildEnumEnumNameStr(name, ind1, ind2):
    return '%s_%d_%d' % (name, ind1, ind2)

def splitEnumEnumNameStr(enumName):
    frontStr, ind2 = splitEnumNameStr(enumName)
    patchName, ind1 = splitEnumNameStr(frontStr)
    return patchName, ind1, ind2

@profile
def mergeNotesFiles(notesPathList, globFlag=False):
    """
    Given a list of path strings pointing to notes files, produce a dict containing all the
    time series info and all 'unusual' info in those notes.
    """
    specialDict = {}
    if globFlag:
        newNotesPathList = []
        for notesFName in notesPathList:
            print '%s yields %s' % (notesFName, glob.glob(notesFName))
            newNotesPathList += glob.glob(notesFName)
        newNotesPathList.sort()
        notesPathList = newNotesPathList
    for ind, notesFName in enumerate(notesPathList):
        notesDict = importNotes(notesFName)
        for nm, dct in notesDict.items():
            if '_' not in nm or nm.startswith('Patch'):
                specialDict[buildEnumNameStr(nm, ind)] = dct
    return specialDict


@profile
def getTimeSeriesList(locKey, specialDict, specialDictKey):
    """
    key specifies the target entity, for example a loc abbrev or a facility category. These are
        typically in uppercase, but that is not enforced.
    specialDict is of the form output by mergeNotesFiles().
    specialDictKey is one of the second-level keys of specialDict, for example 'localpathogen'.
    Returns a list of tuples, of three possible forms.
    
    The first form appears if the time series contains
    entries for multiple different status levels for multiple tiers, for example
    populations sorted by CareTier at different PthLevel values.
    That form is:
    
          [(dayArray, valArrayD), (dayArray, valArrayD), ...]
          
    where dayArray is a numpy array of dates and valArrayD has the form:
    
          {(intLvl, intLvl): fracArray, (intLvl, intLvl): fracArray, ...}
          
    and the ints in the tuple represent the two indices, for example (tier, pthStatus).
          
    The second form appears if the time series contains
    entries for multiple different status levels, for example populations at different PthLevel values.
    That form is:
    
          [(dayArray, valArrayD), (dayArray, valArrayD), ...]
          
    where dayArray is a numpy array of dates and valArrayD has the form:
    
          {intLvl: fracArray, intLvl: fracArray, ...}
          
    with pathLvl being an integer status index (like a PthLevel) and fracArray being a numpy array counts at
    that level.
    
    The second form appears if the time series data is not classified by level, for example a simple population
    count.  That form is:

          [(dayArray, valArray), (dayArray, valArray), ...]
    
    """
    patchList = []
    indList = []
    for enumPatchName in specialDict:
        patchName, ind = splitEnumNameStr(enumPatchName)
        if patchName not in patchList:
            patchList.append(patchName)
        if ind not in indList:
            indList.append(ind)
    patchList.sort()
    indList.sort()
    rsltL = []
    for ind in indList:
        for patchName in patchList:
            enumPatchName = buildEnumNameStr(patchName, ind)
            if enumPatchName not in specialDict or specialDictKey not in specialDict[enumPatchName]:
                continue  # No guarantee all patches have all inds, or the field of interest
            dataList = specialDict[enumPatchName][specialDictKey]
            assert isinstance(dataList, types.ListType), \
                'Special data %s is not a list' % enumPatchName
            fields = defaultdict(list)
            for d in dataList:
                for k, v in d.items():
                    fields[k].append(v)
            assert 'day' in fields, 'Date field is missing for special data %s' % patchName
            dayV = np.array(fields['day'])
            del fields['day']
            
            subKeyL = []
            for keyStr in fields.keys():
                try:
                    base, ind1, ind2 = splitEnumEnumNameStr(keyStr)
                    subKeyL.append((ind1, ind2))
                except ValueError:
                    try:
                        base, ind = splitEnumNameStr(keyStr)
                        subKeyL.append(ind)
                    except ValueError:
                        pass  # This field apparently lacks integer keys
  
            if subKeyL:
                curves = {}
                for subKey in subKeyL:
                    if isinstance(subKey, types.TupleType):
                        ind1, ind2 = subKey
                        key = buildEnumEnumNameStr(locKey, ind1, ind2)
                    else:
                        key = buildEnumNameStr(locKey, subKey)
                    if key in fields:
                        l = fields[key]
                        assert len(l) == len(dayV), (('field %s is the wrong length in special'
                                                         ' data %s (%d vs. %d)')
                                                        % (key, patchName, len(l), len(dayV)))
                        curves[subKey] = np.array(l)
                rsltL.append((dayV, curves))
            else:
                rsltL.append((dayV, fields[locKey]))

    return rsltL

def loadFacilityConstants(constDir, fname):
    """Reads a constants file and returns a list of faility tau values"""
    sys.path.append(constDir)
    fullPath = pyrheautils.pathTranslate(os.path.join(constDir, fname))
    if os.path.exists(fullPath):
        with open(fullPath, 'r') as f:
            #return yaml_tools._simplify(yaml.safe_load(f))
            return yaml.safe_load(f)
    raise RuntimeError('Constants file %s not found' % fullPath)

def getFacilityTargetTau(abbrevList, facTausDict):
    tauDict = {}
    for abbrev in abbrevList:
        for category in facTausDict['facilityTau']:
            for abrv in category['facilities']:
                if abrv['abbrev'] == abbrev:
                    tauDict[abbrev] = {}
                    for t in abrv['tiers']:
                        tauDict[abbrev][t['tier']] = t['frac']['value']

    return tauDict

def shrink(vec):
    """Scale the vector of P values away from the scary singularities at 0 and 1"""
    return ((vec - 0.5) * SHRINK) + 0.5

def logitLink(vec):
    """Given a numpy vector of floats, convert to logit"""
    #oldSettings = np.seterr(divide='ignore')
    rslt = np.log(vec/(1.0-vec))
    #np.seterr(**oldSettings)
    return rslt

def invLogitLink(vec):
    """Given a numpy vector of floats, invert the logit transformation"""
    expTerm = np.exp(vec)
    return expTerm/(expTerm + 1.0)

def calcLogLik(etaMat, etaTVec):
    #oldSettings = np.seterr(invalid='ignore')
    etaMeans = np.nanmean(etaMat, 1)
    # printBadStuff(etaMeans, 'etaMeans')
    etaStdv = np.nanstd(etaMat, 1, ddof=1) + MIN_STDV
    # printBadStuff(etaMeans, 'etaStdvs')
    logLikV = norm.logpdf(etaTVec, loc=etaMeans, scale=etaStdv)
    # printBadStuff(logLikV, 'logLikV')
    # The following replaces any left-over nans with 0.0
    logLikV[logLikV == np.inf] = 0.0
    logLikV = np.nan_to_num(logLikV)
    logLik = np.sum(logLikV)
    #np.seterr(**oldSettings)
    return logLik

def getRatios(abbrevList, specialDict, facDict, waitTime, locAbbrevList):
    ratioList = dict()
    if locAbbrevList:
        locAbbrevSet = set(locAbbrevList)
    for abbrev in abbrevList:
        if locAbbrevList and abbrev not in locAbbrevSet:
            continue
        tplList = getTimeSeriesList(abbrev, specialDict, 'localtierarrivals')
        ratioList[abbrev] = dict()
        for dayVec, curves in tplList:
            tmpVecD = defaultdict(list)
            totVecD = {}
            for tpl, lVec in curves.items():
                tier = tpl
                tmpVecD[tier].append(lVec)
            for tier, lVecList in tmpVecD.items():
                totVecD[tier] = sum(lVecList)
            for tpl, lVec in curves.items():
                tier = tpl
                if(np.sum(totVecD[tier]) > 0.0):
                    #print '%s, %s, %s, %s' % (abbrev, CareTier.names[tier], 
                    #                          np.sum(lVec[dayVec >= waitTime]), np.sum(totVecD[tier]))
                    ratioList[abbrev][CareTier.names[tier]] = np.sum(lVec[dayVec >= waitTime]) / np.sum(totVecD[tier])
                #else:
                    #print('no prevalnce in {}'.format(abbrev))
    return ratioList

def main():
    """
    main
    """
    global logger
    logging.config.dictConfig(getLoggerConfig())
    logger = logging.getLogger(__name__)

    parser = OptionParser(usage="""
    %prog run_descr.yaml --notes notes_file.pkl
    """)

    parser.add_option('-n', '--notes', action='append', type='string',
                      help="Notes filename - may be repeated")
   
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    if not opts.notes:
        parser.error('At least one --notes option is required')

    parser.destroy()

    schemautils.setSchemaBasePath(SCHEMA_DIR)
    inputDict = checkInputFileSchema(args[0],
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    modelDir = inputDict['modelDir']
    pyrheautils.PATH_STRING_MAP['MODELDIR'] = modelDir
    implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])
    pyrheautils.PATH_STRING_MAP['IMPLDIR'] = implDir
    constDir = os.path.join(modelDir, 'constants')
    pyrheautils.PATH_STRING_MAP['CONSTDIR'] = constDir

    facDirList = [pyrheautils.pathTranslate(pth) for pth in inputDict['facilityDirs']]
    facDict = readFacFiles(facDirList)

    facConsts = loadFacilityConstants(constDir, "cre_constants.yaml")

    # Main calibration loop

    for i in range(1): #TODO: Make loop count cli arg; handle signals SIGINT SIGTSTP & cleanly exit
        #print("> run {0}".format(i))


        ### call pyrhea; provide run description and output file
        #run_pyrhea(args[0], "notes_{0}.pkl".format(i))
        

        # Read target tau values from model
        tauDict = getFacilityTargetTau(inputDict['trackedFacilities'], facConsts)

        # Build vector of target taus
        tauVec = list()
        for abbr in tauDict.keys():
           for tier in tauDict[abbr].keys():
               tauVec.append(tauDict[abbr][tier])
        targetLogitVec = logitLink( shrink(np.array(tauVec)) )


        # For all N notes files, build vector of observed taus and make it a row in matrix, m
        m = []
        for n in opts.notes:
            specialDict = mergeNotesFiles([n], False)
            locAbbrevList = None   
            assert 'trackedFacilities' in inputDict, 'No trackedFacilities?'
            waitTime = inputDict['burnInDays'] + inputDict['scenarioWaitDays']
            
            ratioDict = getRatios(inputDict['trackedFacilities'], specialDict, facDict, waitTime, locAbbrevList)

            ratioVec = list()

            # Creates vector for the common set of facilities tracked and facilities with target taus
            for abbr in ratioDict:
                if abbr in tauDict.keys():
                    for tier in ratioDict[abbr]:
                        if tier in tauDict[abbr].keys():
                            #print( ratioDict[abbr][tier] )
                            ratioVec.append(ratioDict[abbr][tier])
            m.append( logitLink( shrink(np.array(ratioVec) ) ) )
        simLogitMat = np.array(m)


        likelihood = calcLogLik(np.transpose(simLogitMat), targetLogitVec)

        print(likelihood)

        # heuristic step to maximize logliklihood
        

    # End calibration loop

if __name__ == "__main__":
    ratios = main()
