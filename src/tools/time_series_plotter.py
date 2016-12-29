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

import os.path
import sys
from optparse import OptionParser
import types
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))

import pyrheautils
import schemautils
import pathogenbase as pth
from notes_plotter import readFacFiles, scanAllFacilities, checkInputFileSchema
from notes_plotter import importNotes
from notes_plotter import SCHEMA_DIR, INPUT_SCHEMA, CARE_TIERS, FAC_TYPE_TO_CATEGORY_MAP

def buildEnumNameStr(name, ind):
    return '%s_%d' % (name, ind)

def splitEnumNameStr(enumName):
    offset = enumName.rfind('_')
    patchName = enumName[:offset]
    ind = int(enumName[offset+1:])
    return patchName, ind


def occupancyTimeFig(specialDict, meanPopByCat=None):
    patchList = []
    indList = []
    for enumPatchName in specialDict:
        patchName, ind = splitEnumNameStr(enumPatchName)
        if patchName not in patchList:
            patchList.append(patchName)
        if ind not in indList:
            indList.append(ind)
        
    figs5, axes5 = plt.subplots(nrows=1, ncols=len(patchList))
    if len(patchList) == 1:
        axes5 = [axes5]
    for offset, patchName in enumerate(patchList):
        clrDict = defaultdict(lambda: None)
        for ind in indList:
            enumPatchName = buildEnumNameStr(patchName, ind)
            if enumPatchName in specialDict:
                data = specialDict[enumPatchName]
                try:
                    occDataList = data['occupancy']
                    assert isinstance(occDataList, types.ListType), \
                        'Special data %s is not a list' % enumPatchName
                    fields = defaultdict(list)
                    for d in occDataList:
                        for k, v in d.items():
                            fields[k].append(v)
                    assert 'day' in fields, 'Date field is missing for special data %s' % patchName
                    dayList = fields['day']
                    del fields['day']
                    keys = fields.keys()
                    keys.sort()
                    for k in keys:
                        l = fields[k]
                        assert len(l) == len(dayList), (('field %s is the wrong length in special data %s'
                                                         '(%d vs. %d)')
                                                        % (k, patchName, len(l), len(dayList)))
                        if not k.endswith('_all'):
                            # The '_all' curves count every arrival at the facility and
                            # are no longer of interest
                            if ind == indList[0]:
                                label = k
                            else:
                                label = None
                            baseLine, = axes5[offset].plot(dayList, l, label=label, c=clrDict[k])
                            clrDict[k] = baseLine.get_color()
                            if (meanPopByCat is not None and k in FAC_TYPE_TO_CATEGORY_MAP
                                    and FAC_TYPE_TO_CATEGORY_MAP[k] in meanPopByCat):
                                meanPop = meanPopByCat[FAC_TYPE_TO_CATEGORY_MAP[k]]
                                axes5[offset].plot(dayList, [meanPop] * len(dayList),
                                                   color=baseLine.get_color(),
                                                   linestyle='--')
                                    
                    axes5[offset].set_xlabel('Days')
                    axes5[offset].set_ylabel('Occupancy')
                    axes5[offset].legend()
                except Exception, e:
                    print e
        axes5[offset].set_title(patchName)
    figs5.tight_layout()
    figs5.canvas.set_window_title("Time History of Occupancy")


def pathogenTimeFig(specialDict):
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
                        if ind == indList[0]:
                            lbl = '%s' % pth.PthStatus.names[pthLvl]
                        else:
                            lbl = None
                        if np.count_nonzero(lVec):
                            with np.errstate(divide='ignore', invalid='ignore'):
                                scaleV = np.true_divide(np.asfarray(lVec), np.asfarray(totsThisCat))
                                scaleV[scaleV == np.inf] = 0.0
                                scaleV = np.nan_to_num(scaleV)
                                thisP = axes6[rowOff, colOff].plot(dayVec, scaleV, label=lbl, 
                                                                   c=clrDict[pthLvl])
                            if ind == indList[0]:
                                clrDict[pthLvl] = thisP[0].get_color()
                    if ind == indList[0]:
                        axes6[rowOff, colOff].set_xlabel('Days')
                        axes6[rowOff, colOff].set_ylabel('Pathogen Prevalence')
                        axes6[rowOff, colOff].legend()
                        axes6[rowOff, colOff].set_title(cat)
            except Exception, e:
                print 'Exception %s' % e
    figs6.tight_layout()
    figs6.canvas.set_window_title("Time History of Infection Status")


def mergeNotesFiles(notesPathList):
    """
    Given a list of path strings pointing to notes files, produce a dict containing all the
    time series info and all 'unusual' info in those notes.
    """
    specialDict = {}
    for ind, notesFName in enumerate(notesPathList):
        notesDict = importNotes(notesFName)
        for nm, dct in notesDict.items():
            if '_' not in nm or nm.startswith('Patch'):
                specialDict[buildEnumNameStr(nm, ind)] = dct
    return specialDict


def getTimeSeriesList(locKey, specialDict, specialDictKey):
    """
    key specifies the target entity, for example a loc abbrev or a facility category. These are
        typically in uppercase, but that is not enforced.
    specialDict is of the form output by mergeNotesFiles().
    specialDictKey is one of the second-level keys of specialDict, for example 'localpathogen'.
    Returns a list of tuples, of two possible forms.
    
    The first form appears if the time series contains
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
            
            intKeyL = []
            for keyKeyPair in fields.keys():
                try:
                    base, ind = splitEnumNameStr(keyKeyPair)
                    intKeyL.append(ind)
                except ValueError:
                    pass  # This field apparently lacks integer keys
            
            if intKeyL:
                curves = {}
                for intKey in intKeyL:
                    key = buildEnumNameStr(locKey, intKey)
                    if key in fields:
                        l = fields[key]
                        assert len(l) == len(dayV), (('field %s is the wrong length in special'
                                                         ' data %s (%d vs. %d)')
                                                        % (key, patchName, len(l), len(dayV)))
                        curves[intKey] = np.array(l)
                rsltL.append((dayV, curves))
            else:
                rsltL.append((dayV, fields[locKey]))

    return rsltL

def oneFacTimeFig(abbrev, specialDict, meanPop=None):
    fig, axes = plt.subplots(2,1)
    axes[0].set_xlabel('Days')
    axes[0].set_ylabel('Pathogen Prevalence')
    axes[0].set_title("%s History of Infection Status" % abbrev)
    axes[1].set_xlabel('Days')
    axes[1].set_ylabel('Occupancy')
    axes[1].set_title("%s History of Occupancy" % abbrev)

    tplList = getTimeSeriesList(abbrev, specialDict, 'localpathogen')
    clrDict = defaultdict(lambda: None)
    for idx, (dayVec, curves) in enumerate(tplList):
        totVec = sum(curves.values())
        scaledCurves = {}
        for pthLvl, lVec in curves.items():
            with np.errstate(divide='ignore', invalid='ignore'):
                scaleV = np.true_divide(np.asfarray(lVec), np.asfarray(totVec))
                scaleV[scaleV == np.inf] = 0.0
                scaleV = np.nan_to_num(scaleV)
            scaledCurves[pthLvl] = scaleV
        for pthLvl, lVec in scaledCurves.items():
            lbl = (None if idx else ('%s' % pth.PthStatus.names[pthLvl]))
            if np.count_nonzero(lVec):
                thisP, = axes[0].plot(dayVec, lVec, label=lbl, c=clrDict[pthLvl])
                clrDict[pthLvl] = thisP.get_color()
    axes[0].legend()

    popTplList = getTimeSeriesList(abbrev, specialDict, 'localoccupancy')
    popClr = None
    for dayVec, popVec in popTplList:
        baseLine, = axes[1].plot(dayVec, popVec, c=popClr)
        popClr = baseLine.get_color()
    if meanPop is not None:
        axes[1].plot(dayVec, [meanPop] * len(dayVec), c=popClr, linestyle='--')

    fig.tight_layout()
    fig.canvas.set_window_title(abbrev)


def main():
    """
    main
    """
    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl] run_descr.yaml
    """)
    parser.add_option('-n', '--notes', action='append', type='string',
                      help="Notes filename - may be repeated")
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    if not len(opts.notes):
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

    specialDict = mergeNotesFiles(opts.notes)

    if "ChicagoLand" in runDesc:
        occupancyTimeFig(specialDict, meanPopByCat=meanPopByCategory)
    else:
        occupancyTimeFig(specialDict) #, meanPopByCat=meanPopByCategory)
    pathogenTimeFig(specialDict)

    if 'trackedFacilities' in inputDict:
        for abbrev in inputDict['trackedFacilities']:
            if abbrev in facDict:
                oneFacTimeFig(abbrev, specialDict, meanPop=facDict[abbrev]['meanPop']['value'])

    plt.show()

if __name__ == "__main__":
    main()
