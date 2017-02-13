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
from collections import defaultdict
import numpy as np
from scipy.stats import norm
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
SAMPLE_SCHEMA = 'tier_time_samples_schema.yaml'
ALT_SAMPLE_SCHEMA = 'time_samples_schema.yaml'

EPSILON = 1.0e-6

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
        

def getTargetVal(abbrev, category, tier):
    return PATHOGEN_CLASS.getEstimatedPrevalence(PthStatus.COLONIZED,
                                                 abbrev, category, tier)

CAT_TIER_TBL = {'HOSPITAL': 'HOSP',
                'SNF': 'NURSING',
                'VSNF': 'SKILNRS',
                'LTACH': 'LTAC'}

def getSampleValues(samplePath, time, facDict):
    """
    Returns a dict structured as dct[abbrev][tierName] = array-of-samples
    """
    tierMode = True
    try:
        sampleTbl = checkInputFileSchema(pyrheautils.pathTranslate(samplePath),
                                         os.path.join(SCHEMA_DIR, SAMPLE_SCHEMA))
    except Exception, e:
        tierMode = False
        sampleTbl = checkInputFileSchema(pyrheautils.pathTranslate(samplePath),
                                         os.path.join(SCHEMA_DIR, ALT_SAMPLE_SCHEMA))

    vecDD = defaultdict(dict)
    for elt in sampleTbl:
        if elt['time'] == time:
            abbrev = elt['abbrev']
            if tierMode:
                for subElt in elt['tiers']:
                    vecDD[abbrev][subElt['tier']] = np.asarray(subElt['samples']['COLONIZED'])
            else:
                tier = CAT_TIER_TBL[facDict[abbrev]['category']]
                vecDD[abbrev][tier] = np.asarray(elt['samples']['COLONIZED'])
            
    if not vecDD:
        raise RuntimeError('No samples in data file identified by %s for time %s'
                           % (samplePath, time))
    return vecDD

def logitLink(vec):
    """Given a numpy vector of floats, convert to logit"""
    return np.log(vec/(1.0-vec))

def invLogitLink(vec):
    """Given a numpy vector of floats, invert the logit transformation"""
    expTerm = np.exp(vec)
    return expTerm/(expTerm + 1.0)

def printBadStuff(vec, label):
    print '%s: %s %s' % (label, vec[np.isnan(vec)], vec[np.isinf(vec)])

def main():
    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)

    global logger
    #logging.config.dictConfig(getLoggerConfig())
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

    loadPthInfo(inputDict['pathogenImplementationDir'])
    facDirs = [pyrheautils.pathTranslate(dct) for dct in inputDict['facilityDirs']]
    facDict = parseFacilityData(facDirs)
    sampDD = getSampleValues(inputDict['sampleFile'], inputDict['testTime'], facDict)
    logLik = 0.0
    testTuples = []
    nSamp = None
    for abbrev, subD in sampDD.items():
        for tierStr, sampV in subD.items():
            if nSamp:
                if sampV.shape[0] != nSamp:
                    raise RuntimeError('Discordant num samples (%s vs %s) for %s %s'
                                       % (sampV.shape[0], nSamp, abbrev, tierStr))
            else:
                nSamp = sampV.shape[0]
            tier = getattr(CareTier, tierStr)
            tVec = np.asarray(getTargetVal(abbrev, facDict[abbrev]['category'], tier))
            testTuples.append((sampV, tVec))
    sampMat = np.asarray([sampV for sampV, tVec in testTuples], dtype=np.double)
    fullTVec = np.asarray([tVec for sampV, tVec in testTuples], dtype=np.double)
    etaMat = logitLink(np.minimum(sampMat + EPSILON, 1.0 - EPSILON))
    etaTVec = logitLink(np.minimum(fullTVec + EPSILON, 1.0 - EPSILON))
    #printBadStuff(etaTVec, 'etaTVec')
    etaMeans = np.mean(etaMat, 1)
    #printBadStuff(etaMeans, 'etaMeans')
    etaStdv = np.std(etaMat, 1, ddof=1)
    #printBadStuff(etaMeans, 'etaStdvs')
    logLikV = norm.logpdf(etaTVec, loc=etaMeans, scale=etaStdv)
    # The following replaces any left-over nans with 0.0
    #printBadStuff(logLikV)
    logLikV[logLikV == np.inf] = 0.0
    logLikV = np.nan_to_num(logLikV)
    logLik = np.sum(logLikV)
    print logLik

if __name__ == "__main__":
    main()
