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
import yaml
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))

import pyrheautils
import schemautils
import pathogenbase as pth
from notes_plotter import readFacFiles, checkInputFileSchema
from notes_plotter import SCHEMA_DIR, INPUT_SCHEMA
from time_series_plotter import mergeNotesFiles, getTimeSeriesList


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


def extractSamples(abbrev, time, specialDict):
    tplList = getTimeSeriesList(abbrev, specialDict, 'localpathogen')
    sampListDict = {pthLvl: [] for pthLvl in pth.PthStatus.names.values()}
    for dayVec, curves in tplList:
        totVec = sum(curves.values())
        scaledCurves = {}
        for pthLvl, lVec in curves.items():
            with np.errstate(divide='ignore', invalid='ignore'):
                scaleV = np.true_divide(np.asfarray(lVec), np.asfarray(totVec))
                scaleV[scaleV == np.inf] = 0.0
                scaleV = np.nan_to_num(scaleV)
            scaledCurves[pthLvl] = scaleV
            sampV = scaleV[dayVec == time]
            if len(sampV):
                sampListDict[pth.PthStatus.names[pthLvl]].append(float(sampV[0]))
    return sampListDict


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

    schemautils.setSchemaBasePath(SCHEMA_DIR)
    inputDict = checkInputFileSchema(args[0],
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    if 'trackedFacilities' not in inputDict or not len(inputDict['trackedFacilities']):
        raise RuntimeError('Run description file does not list any tracked facilities')
    
    modelDir = inputDict['modelDir']
    pyrheautils.PATH_STRING_MAP['MODELDIR'] = os.path.abspath(modelDir)
    implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])
    pyrheautils.PATH_STRING_MAP['IMPLDIR'] = implDir

    facDirList = [pyrheautils.pathTranslate(fPth) for fPth in inputDict['facilityDirs']]
    facDict = readFacFiles(facDirList)

    specialDict = mergeNotesFiles(opts.notes)
    
    sampTimes = [100,  # 2010 end of burn-in
                 815,  # 2012 begins
                 994,  # 2012 second half
                 1168, # 2013 begins
                 1350, # 2013 second half
                 1524, # 2014 begins
                 1880  # 2015 begins
                 ] 

    sampleL = []
    if 'trackedFacilities' in inputDict:
        for abbrev in inputDict['trackedFacilities']:
            if abbrev in facDict:
                for time in sampTimes:
                    sampleL.append({'abbrev': abbrev, 'time': time,
                                    'samples': extractSamples(abbrev, time, specialDict)})

    with open('time_pth_samples.yaml', 'w') as f:
        yaml.safe_dump(sampleL, f, indent=4,
                       encoding='utf-8', width=130, explicit_start=True)

if __name__ == "__main__":
    main()
