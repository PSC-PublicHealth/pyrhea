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
import matplotlib as mpl
mpl.use('svg')
import matplotlib.pyplot as plt

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))

import pyrheautils
import schemautils
import pathogenbase as pth
from facilitybase import CareTier
from notes_plotter import readFacFiles, scanAllFacilities, checkInputFileSchema
from notes_plotter import importNotes
from notes_plotter import SCHEMA_DIR, INPUT_SCHEMA, CARE_TIERS, FAC_TYPE_TO_CATEGORY_MAP
from pyrhea import getLoggerConfig

logger = None

_CACHED_PARSE_DICT = {}
_CACHED_PARSE_ID = None   # the id of the info from which the cached dict was generated

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


def occupancyTimeFig(specialDict, opts, meanPopByCat=None):
    dct = defaultdict(lambda: defaultdict(dict))
    seqSet = set()
    for seqID, patchName, key, dayV, dataV in timeSeriesListGenerator(specialDict, 'occupancy'):
        seqSet.add(seqID)
        subKey, ind1, ind2 = key
        assert ind1 is None and ind2 is None, "Format of occupancy data has changed?"
        if not subKey.endswith('_all'):
            dct[patchName][seqID][subKey] = (dayV, dataV)
    seqL = list(seqSet)
    seqL.sort()
        
    fig, axes = plt.subplots(nrows=1, ncols=len(dct))
    clrDict = {}
    if len(dct) == 1:
        axes = [axes]
    for patchOffset, patchName in enumerate(dct):
        for seq in seqL:
            if seq in dct[patchName]:
                for subKey, (dayV, dataV) in dct[patchName][seq].items():
                    argD = {}
                    if seq == seqL[0]:
                        argD['label'] = subKey
                    if subKey in clrDict:
                        argD['c'] = clrDict[subKey]
                    baseLine, = axes[patchOffset].plot(dayV, dataV, **argD)
                    clrDict[subKey] = baseLine.get_color()
                    if (meanPopByCat is not None and subKey in FAC_TYPE_TO_CATEGORY_MAP
                            and FAC_TYPE_TO_CATEGORY_MAP[subKey] in meanPopByCat):
                        meanPop = meanPopByCat[FAC_TYPE_TO_CATEGORY_MAP[subKey]]
                        axes[patchOffset].plot(dayV, [meanPop] * dayV.shape[0],
                                               color=baseLine.get_color(),
                                               linestyle='--')
        axes[patchOffset].set_xlabel('Days')
        axes[patchOffset].set_ylabel('Occupancy')
        axes[patchOffset].legend()
        axes[patchOffset].set_title(patchName)
    fig.tight_layout()
    fig.canvas.set_window_title("Time History of Occupancy")
    if opts.png:
        plt.savefig(os.path.join(opts.odir, 'zzzall_occupancy.png'))


def pathogenTimeFig(specialDict, opts):
    catList = []
    patchList = []
    indList = []
    for enumPatchName, data in specialDict.items():
        patchName, ind = splitEnumNameStr(enumPatchName)
        if patchName not in patchList:
            patchList.append(patchName)
        if ind not in indList:
            indList.append(ind)
        pthDataList = data['pathogen']
        assert isinstance(pthDataList, types.ListType), 'Special data %s is not a list' % patchName
        for d in pthDataList:
            for k in d.keys():
                if k != 'day':
                    cat = k.split('_')[0]
                    if cat not in catList:
                        catList.append(cat)
    patchList.sort()
    indList.sort()
    catList.sort()
    figs6, axes6 = plt.subplots(nrows=len(catList), ncols=len(patchList))
    axes6.reshape((len(catList), len(patchList)))
    if len(catList) == 1:
        axes6 = axes6[np.newaxis, :]
    if len(patchList) == 1:
        axes6 = axes6[:, np.newaxis]
    for colOff, patchName in enumerate(patchList):
        catOffDict = {cat:row for row, cat in enumerate(catList)}
        clrDict = defaultdict(lambda: None)
        for ind in indList:
            try:
                enumPatchName = buildEnumNameStr(patchName, ind)
                if enumPatchName not in specialDict:
                    continue
                pthDataList = specialDict[enumPatchName]['pathogen']
                assert isinstance(pthDataList, types.ListType), \
                    'Special data %s is not a list' % enumPatchName
                fields = {}
                for d in pthDataList:
                    for k, v in d.items():
                        if k not in fields:
                            fields[k] = []
                        fields[k].append(v)
                assert 'day' in fields, 'Date field is missing for special data %s' % patchName
                dayList = fields['day']
                dayVec = np.array(dayList)
                del fields['day']
                
                curvesThisCat = {}
                for cat in catList:
                    rowOff = catOffDict[cat]
                    for pthLvl in xrange(len(pth.PthStatus.names)):
                        key = buildEnumNameStr(cat, pthLvl)
                        if key in fields:
                            l = fields[key]
                            assert len(l) == len(dayList), (('field %s is the wrong length in special'
                                                             ' data %s (%d vs. %d)')
                                                            % (k, patchName, len(l), len(dayList)))
                            curvesThisCat[pthLvl] = np.array(l)
                    totsThisCat = sum(curvesThisCat.values())
                    for pthLvl, lVec in curvesThisCat.items():
                        argD = {}
                        if ind == indList[0]:
                            argD['label'] = k
                        if clrDict[k] is not None:
                            argD['c'] = clrDict[pthLvl]
                        if np.count_nonzero(lVec):
                            with np.errstate(divide='ignore', invalid='ignore'):
                                scaleV = np.true_divide(np.asfarray(lVec), np.asfarray(totsThisCat))
                                scaleV[scaleV == np.inf] = 0.0
                                scaleV = np.nan_to_num(scaleV)
                                thisP = axes6[rowOff, colOff].plot(dayVec, scaleV, **argD)
                            if ind == indList[0]:
                                clrDict[pthLvl] = thisP[0].get_color()
                    if ind == indList[0]:
                        axes6[rowOff, colOff].set_xlabel('Days')
                        axes6[rowOff, colOff].set_ylabel('Pathogen Prevalence')
                        #axes6[rowOff, colOff].legend()
                        axes6[rowOff, colOff].set_title(cat)
            except Exception, e:
                print 'Exception %s' % e
    #figs6.tight_layout()
    figs6.canvas.set_window_title("Time History of Infection Status")
    if opts.png:
        plt.savefig(os.path.join(opts.odir, 'zzzall_pathogen.png'))

def pathogenTierTimeFig(specialDict, opts):
    catList = []
    patchList = []
    indList = []
    for enumPatchName, data in specialDict.items():
        patchName, ind = splitEnumNameStr(enumPatchName)
        if patchName not in patchList:
            patchList.append(patchName)
        if ind not in indList:
            indList.append(ind)
        pthDataList = data['pathogen']
        assert isinstance(pthDataList, types.ListType), 'Special data %s is not a list' % patchName
        for d in pthDataList:
            for k in d.keys():
                if k != 'day':
                    cat = k.split('_')[0]
                    if cat not in catList:
                        catList.append(cat)
    patchList.sort()
    indList.sort()
    catList.sort()
    figs6, axes6 = plt.subplots(nrows=len(catList), ncols=len(patchList))
    axes6.reshape((len(catList), len(patchList)))
    if len(catList) == 1:
        axes6 = axes6[np.newaxis, :]
    if len(patchList) == 1:
        axes6 = axes6[:, np.newaxis]
    for colOff, patchName in enumerate(patchList):
        catOffDict = {cat:row for row, cat in enumerate(catList)}
        clrDict = defaultdict(lambda: None)
        for ind in indList:
            try:
                enumPatchName = buildEnumNameStr(patchName, ind)
                if enumPatchName not in specialDict:
                    continue
                pthDataList = specialDict[enumPatchName]['pathogen']
                assert isinstance(pthDataList, types.ListType), \
                    'Special data %s is not a list' % enumPatchName
                fields = {}
                for d in pthDataList:
                    for k, v in d.items():
                        if k not in fields:
                            fields[k] = []
                        fields[k].append(v)
                assert 'day' in fields, 'Date field is missing for special data %s' % patchName
                dayList = fields['day']
                dayVec = np.array(dayList)
                del fields['day']
                
                curvesThisCat = {}
                for cat in catList:
                    rowOff = catOffDict[cat]
                    for pthLvl in xrange(len(pth.PthStatus.names)):
                        key = buildEnumNameStr(cat, pthLvl)
                        if key in fields:
                            l = fields[key]
                            assert len(l) == len(dayList), (('field %s is the wrong length in special'
                                                             ' data %s (%d vs. %d)')
                                                            % (k, patchName, len(l), len(dayList)))
                            curvesThisCat[pthLvl] = np.array(l)
                    totsThisCat = sum(curvesThisCat.values())
                    for pthLvl, lVec in curvesThisCat.items():
                        argD = {}
                        if ind == indList[0]:
                            argD['label'] = k
                        if clrDict[k] is not None:
                            argD['c'] = clrDict[pthLvl]
                        if np.count_nonzero(lVec):
                            with np.errstate(divide='ignore', invalid='ignore'):
                                scaleV = np.true_divide(np.asfarray(lVec), np.asfarray(totsThisCat))
                                scaleV[scaleV == np.inf] = 0.0
                                scaleV = np.nan_to_num(scaleV)
                                thisP = axes6[rowOff, colOff].plot(dayVec, scaleV, **argD)
                            if ind == indList[0]:
                                clrDict[pthLvl] = thisP[0].get_color()
                    if ind == indList[0]:
                        axes6[rowOff, colOff].set_xlabel('Days')
                        axes6[rowOff, colOff].set_ylabel('Pathogen Prevalence')
                        #axes6[rowOff, colOff].legend()
                        axes6[rowOff, colOff].set_title(cat)
            except Exception, e:
                print 'Exception %s' % e
    #figs6.tight_layout()
    figs6.canvas.set_window_title("Time History of Infection Status")
    if opts.png:
        plt.savefig(os.path.join(opts.odir, 'zzzall_pathogen.png'))


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
    if notesPathList:
        for ind, notesFName in enumerate(notesPathList):
            notesDict = importNotes(notesFName)
            for nm, dct in notesDict.items():
                if '_' not in nm or nm.startswith('Patch'):
                    specialDict[buildEnumNameStr(nm, ind)] = dct
    return specialDict


def timeSeriesListGenerator(specialDict, specialDictKey):
    """
    This generator returns all the time series data in specialDict having the associated
    top-level key (for example 'occupancy' or 'localtierpathogen'.  For each time series,
    the data returned is:

    (seqId, (locKey, ind1, ind2), dayVec, valVec)

    where seqId is the sequence number in order of the series of notes.pkl files being
    parsed, dayVec and valVec are numpy arrays containing date and value vectors,
    and the tuple has one of the following forms:

    category, None, None         for aggregate info like 'occupancy'

    category, lvl, None          for aggregate info like 'pathogen'

    category, tier, None         for aggregate info like 'localtierarrivals'

    abbrev, lvl, None            for aggregate info like 'localpathogen'

    abbrev, tier, lvl            for aggregate info like 'localtierpathogen'

    where category is for example 'NursingHome', lvl is an integer PthLevel, and tier is
    an integer CareTier.
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
            
            for key, fld in fields.items():
                fld = np.array(fld)
                try:
                    loc, subI1, subI2 = splitEnumEnumNameStr(key)
                except ValueError:
                    try:
                        loc, subI1 = splitEnumNameStr(key)
                        subI2 = None
                    except ValueError:
                        loc = key
                        subI1 = subI2 = None
                yield (ind, patchName, (loc, subI1, subI2), dayV, fld)

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

          {(intLvl, intLvl): valArray, (intLvl, intLvl): valArray, ...}

    and the ints in the tuple represent the two indices, for example (tier, pthStatus).

    The second form appears if the time series contains
    entries for multiple different status levels, for example populations at different PthLevel values.
    That form is:

          [(dayArray, valArrayD), (dayArray, valArrayD), ...]

    where dayArray is a numpy array of dates and valArrayD has the form:

          {intLvl: valArray, intLvl: valArray, ...}

    with pathLvl being an integer status index (like a PthLevel) and fracArray being a numpy array counts at
    that level.

    The third form appears if the time series data is not classified by level, for example a simple population
    count.  That form is:

          [(dayArray, valArray), (dayArray, valArray), ...]

    """

    global _CACHED_PARSE_DICT
    global _CACHED_PARSE_ID

    if _CACHED_PARSE_ID != id(specialDict) or specialDictKey not in _CACHED_PARSE_DICT:
        _CACHED_PARSE_ID = id(specialDict)
        thisCD = _CACHED_PARSE_DICT[specialDictKey] = defaultdict(list)
        for tpl in timeSeriesListGenerator(specialDict, specialDictKey):
            notesInd, patchName, (loc, ind1, ind2), dayV, dataV = tpl
            thisCD[loc].append((notesInd, ind1, ind2, dayV, dataV))

    tplL = _CACHED_PARSE_DICT[specialDictKey][locKey]
    resultD = {}
    for tpl in tplL:
        notesInd, ind1, ind2, dayV, dataV = tpl
        if ind1 is None:
            assert notesInd not in resultD, 'Generator returned duplicate notes index?'
            resultD[notesInd] = (dayV, dataV)
        else:
            innerKey = ind1 if ind2 is None else (ind1, ind2)
            if notesInd not in resultD:
                resultD[notesInd] = (dayV, {innerKey: dataV})
            else:
                oldDayV, innerD = resultD[notesInd]
                assert dayV is oldDayV, 'Generator provided inconsistent day vectors?'
                innerD[innerKey] = dataV
    return resultD.values()

def oneFacArrivalsFig(abbrev, specialDict, opts, meanPop=None):
    print abbrev
    if abbrev == 'CHOC':
        return
    if abbrev == 'CMIS':
        return
    
    tierMode = True
    tplList = getTimeSeriesList(abbrev, specialDict, 'localtierCREBundle')
    print tplList
    if not tplList:
        tierMode = False
        tplList = getTimeSeriesList(abbrev, specialDict, 'localpathogen')
        
    if tierMode:
        tierAxisD = {}
        tAOffset = 0
        for dayVec, curves in tplList:
            print curves
            for tier in curves:
                if tier not in tierAxisD:
                    tierAxisD[tier] = tAOffset
                    tAOffset += 1
    else:
        tierAxisD = {None: 0}

def oneFacTimeFig(abbrev, specialDict, opts, meanPop=None):
    print abbrev
    if abbrev == 'CHOC':
        return
    if abbrev == 'CMIS':
        return
    
    tierMode = True
    tplList = getTimeSeriesList(abbrev, specialDict, 'localtierpathogen')
    if not tplList:
        tierMode = False
        tplList = getTimeSeriesList(abbrev, specialDict, 'localpathogen')
        
    if tierMode:
        tierAxisD = {}
        tAOffset = 0
        for dayVec, curves in tplList:
            for tier, pthStatus in curves:
                if tier not in tierAxisD:
                    tierAxisD[tier] = tAOffset
                    tAOffset += 1
    else:
        tierAxisD = {None: 0}

    fig, axes = plt.subplots(len(tierAxisD) + 1,1)
    popAxis = axes[-1]
    popAxis.set_xlabel('Days')
    popAxis.set_ylabel('Occupancy')
    popAxis.set_title("%s all tiers History of Occupancy" % abbrev)
    for tier, aOffset in tierAxisD.items():
        axes[aOffset].set_xlabel('Days')
        axes[aOffset].set_ylabel('Pathogen Prevalence')
        if tier is None:
            axes[aOffset].set_title("%s all tiers History of Infection Status" % abbrev)
        else:
            axes[aOffset].set_title("%s %s History of Infection Status" % (abbrev, CareTier.names[tier]))
        axes[aOffset].set_xlabel('Days')
        axes[aOffset].set_ylabel('Pathogen Prevalence')
        
    clrDict = defaultdict(lambda: None)
    if len(tplList) > 3:
        alpha = 0.2
    else:
        alpha = 1.0
    for idx, (dayVec, curves) in enumerate(tplList):
        tmpVecD = defaultdict(list)
        totVecD = {}
        for tpl, lVec in curves.items():
            tier, pthStatus = (tpl if tierMode else (None, tpl))
            tmpVecD[tier].append(lVec)
        for tier, lVecList in tmpVecD.items():
            totVecD[tier] = sum(lVecList)
        for tpl, lVec in curves.items():
            tier, pthStatus = (tpl if tierMode else (None, tpl))
            with np.errstate(divide='ignore', invalid='ignore'):
                scaleV = np.true_divide(np.asfarray(lVec), np.asfarray(totVecD[tier]))
                scaleV[scaleV == np.inf] = 0.0
                scaleV = np.nan_to_num(scaleV)
            if tier is None:
                lblStr = '%s' % pth.PthStatus.names[pthStatus]
            else:
                lblStr = '%s %s' % (CareTier.names[tier], pth.PthStatus.names[pthStatus])
            argD = {'alpha': alpha}
            if idx:
                argD['label'] = lblStr
            if clrDict[tpl] is not None:
                argD['c'] = clrDict[tpl]
            if np.count_nonzero(lVec):
                thisP, = axes[tierAxisD[tier]].plot(dayVec, scaleV, **argD)
                clrDict[tpl] = thisP.get_color()
    for tAOffset in tierAxisD.values():
        #axes[tAOffset].legend()
        pass

    popTplList = getTimeSeriesList(abbrev, specialDict, 'localoccupancy')
    popClr = None
    if len(popTplList) > 3:
        alpha = 0.2
    else:
        alpha = 1.0
    for dayVec, popVec in popTplList:
        argD = {'alpha': alpha}
        if popClr is not None:
            argD['c'] = popClr
        baseLine, = popAxis.plot(dayVec, popVec, **argD)
        popClr = baseLine.get_color()
    if popTplList and meanPop is not None:
        dayVec = popTplList[0][0]
        popAxis.plot(dayVec, [meanPop] * len(dayVec), c=popClr, linestyle='--')

    fig.tight_layout()
    fig.canvas.set_window_title(abbrev)
    if opts.png:
        plt.savefig(os.path.join(opts.odir, '%s.png' % abbrev))


def main():
    """
    main
    """
    global logger
    logging.config.dictConfig(getLoggerConfig())
    logger = logging.getLogger(__name__)
    
    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl] [--glob] run_descr.yaml
    """)
    parser.add_option('-n', '--notes', action='append', type='string',
                      help="Notes filename - may be repeated")
    parser.add_option('--glob', action='store_true',
                      help="Apply filename globbing to the given notes files")
    parser.add_option('--png', action='store_true',
                      help="Save plots to current directory rather than viewing interactively")
    parser.add_option('--odir', action='store', default=os.getcwd(),
                      help="Save any output to the specified directory rather than the current directory")
    
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    if not opts.notes:
        parser.error('At least one --notes option is required')

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
    facDict = readFacFiles(facDirList)

    if "ChicagoLand" in runDesc:
        allOfCategoryFacilityInfo, meanPopByCategory = scanAllFacilities(facDirList,
                                                                         facDict=facDict)  # @UnusedVariable

    specialDict = mergeNotesFiles(opts.notes, opts.glob)

    if "ChicagoLand" in runDesc:
        occupancyTimeFig(specialDict, opts, meanPopByCat=meanPopByCategory)
    else:
        occupancyTimeFig(specialDict, opts) #, meanPopByCat=meanPopByCategory)
    pathogenTimeFig(specialDict, opts)

    if 'trackedFacilities' in inputDict:
        for abbrev in inputDict['trackedFacilities']:
            if abbrev in facDict:
                if 'meanPop' in facDict[abbrev]:
                    oneFacTimeFig(abbrev, specialDict, opts, meanPop=facDict[abbrev]['meanPop']['value'])
                else:
                    oneFacTimeFig(abbrev, specialDict, opts)
                oneFacArrivalsFig(abbrev, specialDict, opts)

    if not opts.png:
        plt.show()

if __name__ == "__main__":
    main()
