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
from time_series_plotter import mergeNotesFiles
from time_series_plotter import buildEnumNameStr, splitEnumNameStr
from time_series_plotter import buildEnumEnumNameStr, splitEnumEnumNameStr
from time_series_plotter import getTimeSeriesList

logger = None

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
