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

import os
import os.path
import sys
import logging
import logging.config
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
# import matplotlib as mpl
# mpl.use('svg')
import matplotlib.pyplot as plt

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
from time_series_plotter import getTimeSeriesList, timeSeriesListGenerator

logger = None

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
def oldGetTimeSeriesList(locKey, specialDict, specialDictKey):
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


def sortSelectedPthStatus(abbrevSet, specialDict, pthStatus, notesKey):
    dayVecD = {}
    totVecD = {}
    for tpl in timeSeriesListGenerator(specialDict, notesKey):
        notesInd, patchName, (loc, tier, pthLvl), dayV, dataV = tpl
        if pthLvl != pthStatus:
            continue
        if tier not in dayVecD:
            dayVecD[tier] = {}
        if notesInd not in dayVecD[tier]:
            dayVecD[tier][notesInd] = dayV
        assert dayVecD[tier][notesInd] is dayV, 'Inconsistent day vectors?'
        if tier not in totVecD:
            totVecD[tier] = {}
        if notesInd in totVecD[tier]:
            totVecD[tier][notesInd] += dataV
        else:
            totVecD[tier][notesInd] = np.copy(dataV)
    
    return dayVecD, totVecD


def sortTiers(abbrevSet, specialDict, notesKey):
    dayVecD = {}
    totVecD = {}
    for tpl in timeSeriesListGenerator(specialDict, notesKey):
        notesInd, patchName, (loc, tier, pthLvl), dayV, dataV = tpl
        assert pthLvl is None, "I am confused about this format"
        if tier not in dayVecD:
            dayVecD[tier] = {}
        if notesInd not in dayVecD[tier]:
            dayVecD[tier][notesInd] = dayV
        assert dayVecD[tier][notesInd] is dayV, 'Inconsistent day vectors?'
        if tier not in totVecD:
            totVecD[tier] = {}
        if notesInd in totVecD[tier]:
            totVecD[tier][notesInd] += dataV
        else:
            totVecD[tier][notesInd] = np.copy(dataV)
    
    return dayVecD, totVecD


def plotPrevalenceFig(abbrevList, specialDictA, specialDictB, facDict):
    abbrevSet = set(abbrevList)
    dayVecAD, totVecAD = sortSelectedPthStatus(abbrevSet, specialDictA, 
                                               PthStatus.COLONIZED, 'localtierpathogen')
    dayVecBD, totVecBD = sortSelectedPthStatus(abbrevSet, specialDictB,
                                               PthStatus.COLONIZED, 'localtierpathogen')

    tierSet = set(dayVecAD.keys() + dayVecBD.keys())
    tierL = list(tierSet)
    tierL.sort()

    fig, axes = plt.subplots(len(tierL),1)
    mAx = {tierL[idx]: axes[idx] for idx in xrange(len(tierL))}
    for tier in tierL:
        mAx[tier].set_xlabel('Days')
        mAx[tier].set_ylabel('Prevalence')
        mAx[tier].set_title("%s Prevalence Across All Facilities" % CareTier.names[tier])

        clr = None
        if min(len(totVecAD[tier]), len(totVecBD[tier])) <= 3:
            alpha = 1.0
        else:
            alpha = 0.2
        argD = {'alpha': alpha}
        for noteInd in totVecAD[tier]:
            pltP, = mAx[tier].plot(dayVecAD[tier][noteInd], totVecAD[tier][noteInd], **argD)
            argD['c'] = pltP.get_color()
        argD = {'alpha': alpha}
        for noteInd in totVecBD[tier]:
            pltP, = mAx[tier].plot(dayVecBD[tier][noteInd], totVecBD[tier][noteInd], **argD)
            argD['c'] = pltP.get_color()

    fig.tight_layout()
    #plt.savefig('plot1.png')


def plotCPFig(abbrevList, specialDictA, specialDictB, facDict):
    abbrevSet = set(abbrevList)
    dayVecAD, totVecAD = sortTiers(abbrevSet, specialDictA, 'localtierCP')
    dayVecBD, totVecBD = sortTiers(abbrevSet, specialDictB, 'localtierCP')

    tierSet = set(dayVecAD.keys() + dayVecBD.keys())
    tierL = list(tierSet)
    tierL.sort()

    fig, axes = plt.subplots(len(tierL),1)
    mAx = {tierL[idx]: axes[idx] for idx in xrange(len(tierL))}
    for tier in tierL:
        mAx[tier].set_xlabel('Days')
        mAx[tier].set_ylabel('Patients On CP')
        mAx[tier].set_title("%s Contact Precautions Across All Facilities" % CareTier.names[tier])

        clr = None
        if min(len(totVecAD[tier]), len(totVecBD[tier])) <= 3:
            alpha = 1.0
        else:
            alpha = 0.2
        argD = {'alpha': alpha}
        for noteInd in totVecAD[tier]:
            pltP, = mAx[tier].plot(dayVecAD[tier][noteInd], totVecAD[tier][noteInd], **argD)
            argD['c'] = pltP.get_color()
        argD = {'alpha': alpha}
        for noteInd in totVecBD[tier]:
            pltP, = mAx[tier].plot(dayVecBD[tier][noteInd], totVecBD[tier][noteInd], **argD)
            argD['c'] = pltP.get_color()

    fig.tight_layout()
    #plt.savefig('plot1.png')


def plotNewColFig(abbrevList, specialDictA, specialDictB, facDict):
    abbrevSet = set(abbrevList)
    dayVecAD, totVecAD = sortTiers(abbrevSet, specialDictA, 'localtiernewcolonized')
    dayVecBD, totVecBD = sortTiers(abbrevSet, specialDictB, 'localtiernewcolonized')

    tierSet = set(dayVecAD.keys() + dayVecBD.keys())
    tierL = list(tierSet)
    tierL.sort()

    fig, axes = plt.subplots(len(tierL),1)
    mAx = {tierL[idx]: axes[idx] for idx in xrange(len(tierL))}
    for tier in tierL:
        mAx[tier].set_xlabel('Days')
        mAx[tier].set_ylabel('New Colonizations')
        mAx[tier].set_title("%s New Colonizations Across All Facilities" % CareTier.names[tier])

        clr = None
        if min(len(totVecAD[tier]), len(totVecBD[tier])) <= 3:
            alpha = 1.0
        else:
            alpha = 0.2
        argD = {'alpha': alpha}
        for noteInd in totVecAD[tier]:
            pltP, = mAx[tier].plot(dayVecAD[tier][noteInd], totVecAD[tier][noteInd], **argD)
            argD['c'] = pltP.get_color()
        argD = {'alpha': alpha}
        for noteInd in totVecBD[tier]:
            pltP, = mAx[tier].plot(dayVecBD[tier][noteInd], totVecBD[tier][noteInd], **argD)
            argD['c'] = pltP.get_color()

    fig.tight_layout()
    #plt.savefig('plot1.png')



def plotStuff(abbrevList, specialDictA, specialDictB, facDict):
    plotPrevalenceFig(abbrevList, specialDictA, specialDictB, facDict)
    plotCPFig(abbrevList, specialDictA, specialDictB, facDict)
    plotNewColFig(abbrevList, specialDictA, specialDictB, facDict)
# 
#     cpAxis.set_xlabel('Days')
#     cpAxis.set_ylabel('Total On CP')
#     cpAxis.set_title("Patients On Contact Precautions Across All Facilities")
        
#     newColP, = newColAxis.plot(dayVec2, totNewColV)
#     prvP, = prvAxis.plot(dayVec1, totColV / totPopV)
#     cpP, = cpAxis.plot(dayVec3, totCPV)



def main():
    """
    main
    """
    global logger
    logging.config.dictConfig(getLoggerConfig())
    logger = logging.getLogger(__name__)
    
    parser = OptionParser(usage="""
    %prog --notes notes_file.pkl run_descr.yaml
    """)
    parser.add_option('-a', '--groupa', action='append', type='string',
                      help="Notes filename to add to group a- may be repeated")

    parser.add_option('-b', '--groupb', action='append', type='string',
                      help="Notes filename to add to group b- may be repeated")
    
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    if not opts.groupa:
        parser.error('At least one -a option is required')
    if not opts.groupb:
        parser.error('At least one -b option is required')

    parser.destroy()

    runDesc = args[0]
        
    schemautils.setSchemaBasePath(SCHEMA_DIR)
    inputDict = checkInputFileSchema(args[0],
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    modelDir = inputDict['modelDir']
    pyrheautils.PATH_STRING_MAP['MODELDIR'] = modelDir
    implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])
    pyrheautils.PATH_STRING_MAP['IMPLDIR'] = implDir

    facDirList = [pyrheautils.pathTranslate(pth) for pth in inputDict['facilityDirs']]
    #facDict = readFacFiles(facDirList)
    facDict = None

    specialDictA = mergeNotesFiles(opts.groupa, False)
    specialDictB = mergeNotesFiles(opts.groupb, False)

    assert 'trackedFacilities' in inputDict, 'No trackedFacilities?'
    plotStuff(inputDict['trackedFacilities'], specialDictA, specialDictB, facDict)
    plt.show()

if __name__ == "__main__":
    main()
