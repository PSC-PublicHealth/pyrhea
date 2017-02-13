#! /usr/bin/env python

###################################################################################
# Copyright   2017, Pittsburgh Supercomputing Center (PSC).  All Rights Reserved. #
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

"""
This program runs somewhere on a cluster host, submitting batch jobs and reading
results to find the parameter point which minimizes a scalar measure.  The use
case is back-calculating run parameters to match a given set of measured values.
"""

import sys
import os.path
import signal
import optparse
import yaml
import numpy as np
import logging
import logging.config

import schemautils
import pyrheautils
from pyrhea import getLoggerConfig, checkInputFileSchema, loadPathogenImplementations
from facilitybase import CareTier, PthStatus
from map_transfer_matrix import parseFacilityData

import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools

SCHEMA_DIR = '../schemata'
INPUT_SCHEMA = 'parameter_optimizer_input_schema.yaml'
SAMPLE_SCHEMA = 'time_samples_schema.yaml'

logger = None

PATHOGEN_CLASS = None

def loadPthInfo(pthImplDir):
    global PATHOGEN_CLASS
    if not PATHOGEN_CLASS:
        logger.info('Loading pathogen implementation')
        modList = []
        implementationDir = pyrheautils.pathTranslate(pthImplDir)
        pyrheautils.PATH_STRING_MAP['POLICYDIR'] = implementationDir
        for newMod in pyrheautils.loadModulesFromDir(implementationDir,
                                                     requiredAttrList=['getPathogenClass']):
            modList.append(newMod)
            newPathogenClass = newMod.getPathogenClass()
            logger.info('Loaded %s implementing %s' % (newMod.__name__,
                                                       newPathogenClass.__name__))
        assert len(modList) == 1, 'Found more than one pathogen implementation in %s' % pthImplDir
        PATHOGEN_CLASS = modList[0].getPathogenClass()
        logger.info('Finished loading pathogen')
    p = PATHOGEN_CLASS.getEstimatedPrevalence(PthStatus.COLONIZED, 'MANO_900_S', 'SNF', CareTier.NURSING)
        

def getTargetVal(abbrev, tier, facDict):
    return PATHOGEN_CLASS.getEstimatedPrevalence(PthStatus.COLONIZED, abbrev, facDict['abbrev']['category'],
                                                 CareTier.HOSP)

def identifierToPath(sampleIdentifier):
    return '/home/welling/git/pyrhea/src/sim/time_samples_wide_pass9.yaml'

def getSampleValues(sampleIdentifier, time):
    """
    Return an N by M array of floats, where N is the number of values in a sample and
    M is the number of sample instances.
    """
    sampleTbl = checkInputFileSchema(identifierToPath(sampleIdentifier),
                                     os.path.join(SCHEMA_DIR, SAMPLE_SCHEMA))
    vecD = {}
    for elt in sampleTbl:
        if elt['time'] == time:
            abbrev = elt['abbrev']
            assert abbrev not in vecD, 'Redundant record for %s at time %s' % (abbrev, time)
            vecD[abbrev] = np.asarray(elt['samples']['COLONIZED'])
    if not vecD:
        raise RuntimeError('No samples in data file identified by %s for time %s'
                           % (sampleIdentifier, time))
    
    return vecD

def main():
    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)

    global logger
    logging.config.dictConfig(getLoggerConfig())
    logger = logging.getLogger(__name__)

    parser = optparse.OptionParser(usage="""
    %prog run_descr.yaml
    """)

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    parser.destroy()

    inputDict = checkInputFileSchema(args[0],
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    schemautils.setSchemaBasePath(SCHEMA_DIR)
    if 'modelDir' in inputDict:
        pyrheautils.PATH_STRING_MAP['MODELDIR'] = pyrheautils.pathTranslate(inputDict['modelDir'])
    if 'pathTranslations' in inputDict:
        for elt in inputDict['pathTranslations']:
            pyrheautils.PATH_STRING_MAP[elt['key']] = elt['value']

    loadPthInfo(inputDict['pathogenImplementationDir'], 'CRECore')
    facDirs = [pyrheautils.pathTranslate(dct) for dct in inputDict['facilityDirs']]
    facDict = parseFacilityData(facDirs)
    sampD = getSampleValues('someIdentifier', 465)
    for abbrev, valV in sampD.items():
        tier = 0
        tV = getTargetVal(abbrev, facDict[abbrev]['category'], tier)

if __name__ == "__main__":
    main()
