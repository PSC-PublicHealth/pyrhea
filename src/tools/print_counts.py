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


def printStuff(abbrevList, specialDict, facDict):
    for abbrev in abbrevList:
        tplList = getTimeSeriesList(abbrev, specialDict, 'localtierpathogen')
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
                    print '%s, %s, %s, %s' % (abbrev, CareTier.names[tier], np.sum(lVec), np.sum(totVecD[tier]))
                    found = True
            if not found:
                print 'No Data For %s' % abbrev

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
    parser.add_option('-n', '--notes', action='append', type='string',
                      help="Notes filename - may be repeated")
    
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    if not opts.notes or len(opts.notes)>1:
        parser.error('One --notes option is required')

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

    specialDict = mergeNotesFiles(opts.notes, False)

    assert 'trackedFacilities' in inputDict, 'No trackedFacilities?'
    printStuff(inputDict['trackedFacilities'], specialDict, facDict)

if __name__ == "__main__":
    main()
