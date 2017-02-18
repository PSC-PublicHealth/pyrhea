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

import numpy as np

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))

import pyrheautils
import schemautils
import pathogenbase as pth
from notes_plotter import readFacFiles, checkInputFileSchema
from notes_plotter import SCHEMA_DIR, INPUT_SCHEMA
from time_series_plotter import mergeNotesFiles, getTimeSeriesList

DEFAULT_OUT_FILE = 'time_samples.yaml'
OUTPUT_SCHEMA = 'time_samples_schema.yaml'


def extractSamples(abbrev, time, specialDict):
    pthTplList = getTimeSeriesList(abbrev, specialDict, 'localpathogen')
    sampListDict = {pthLvl: [] for pthLvl in pth.PthStatus.names.values()}
    sampListDict['occupancy'] = []
    for dayVec, curves in pthTplList:
        totVec = sum(curves.values())
        for pthLvl, lVec in curves.items():
            with np.errstate(divide='ignore', invalid='ignore'):
                scaleV = np.true_divide(np.asfarray(lVec), np.asfarray(totVec))
                scaleV[scaleV == np.inf] = 0.0
                scaleV = np.nan_to_num(scaleV)
            sampV = scaleV[dayVec == time]
            if len(sampV):
                sampListDict[pth.PthStatus.names[pthLvl]].append(float(sampV[0]))
    popTplList = getTimeSeriesList(abbrev, specialDict, 'localoccupancy')
    for dayVec, popVec in popTplList:
        popV = np.asfarray(popVec)[dayVec == time]
        if len(popV): 
            sampListDict['occupancy'].append(float(popV[0]))
    popNcList = getTimeSeriesList(abbrev, specialDict, 'NEW_COLONIZED')
    return sampListDict


def main():
    """
    main
    """
    import sys
    print sys.argv
    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl] [--glob] [--proto prototype.yaml] [--out outname.yaml] run_descr.yaml
    """)
    parser.add_option('-n', '--notes', action='append', type='string',
                      help="Notes filename - may be repeated")
    
    parser.add_option('--proto', action='store', type='string',
                      help="Prototype YAML file to provide sampling times")
    
    parser.add_option('-o', '--out', action='store', type='string',
                      help="Output file (defaults to %s)" % DEFAULT_OUT_FILE)
    
    parser.add_option('--glob', action='store_true',
                      help=("Apply filename globbing for notes files."
                            "  (Remember to protect the filename string from the shell!)"))
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    if not opts.notes:
        parser.error('At least one --notes option is required')

    parser.destroy()

    schemautils.setSchemaBasePath(SCHEMA_DIR)
    
    if opts.out:
        outFileName = opts.out
    else:
        outFileName = DEFAULT_OUT_FILE
    
    sampTimeSet = set()
    if opts.proto:
        protoJSON = checkInputFileSchema(opts.proto,
                                         os.path.join(SCHEMA_DIR, OUTPUT_SCHEMA))
        for ent in protoJSON:
            sampTimeSet.add(ent['time'])
        sampTimes = list(sampTimeSet)
        sampTimes.sort()
    else:
        sampTimes = [x for x in range(100,200)]
        #sampTimes = [100,  # 2010 end of burn-in
        #             815,  # 2012 begins
        #             994,  # 2012 second half
        #             1168, # 2013 begins
        #             1350, # 2013 second half
        #             1524, # 2014 begins
        #             1880  # 2015 begins
        #             ] 
    
    print sampTimes
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

    specialDict = mergeNotesFiles(opts.notes, opts.glob)
    
    sampleL = []
    if 'trackedFacilities' in inputDict:
        for abbrev in inputDict['trackedFacilities']:
            if abbrev in facDict:
                if abbrev.find("RUSH") > -1:
                    for time in sampTimes:
                        print "abbrev: {0} time: {1}".format(abbrev,time)
                        sampleL.append({'abbrev': abbrev, 'time': time,
                                        'samples': extractSamples(abbrev, time, specialDict)})

    with open(outFileName, 'w') as f:
        yaml.safe_dump(sampleL, f, indent=4,
                       encoding='utf-8', width=130, explicit_start=True)

if __name__ == "__main__":
    main()
