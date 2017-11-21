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

from __future__ import print_function
import future
import six
import os.path
import sys
from optparse import OptionParser
from collections import defaultdict
import yaml
import logging
import logging.config
import pika

import numpy as np

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))

import pyrheautils
import schemautils
import stats
from notes_plotter import readFacFiles, checkInputFileSchema, fullCRVFromLOSModel
from notes_plotter import SCHEMA_DIR, INPUT_SCHEMA
from pyrhea import getLoggerConfig

"""
Facilities to exclude, typically because they violate the schema (which isn't detected here)
"""
EXCLUDE_THESE_FACILITIES = ['CHOC', 'HBRI', 'HSUB', 'CMIS', 'AGBP']

DAYS_PER_YEAR = 365  # But might not, in HERMES for example
#DAYS_PER_YEAR = 10  # But might not, in HERMES for example

def configureLogging(cfgDict):
    global LOGGER
    logging.config.dictConfig(cfgDict)
    LOGGER = logging.getLogger(__name__)

class Mutator(object):
    @classmethod
    def select(cls, state, facIdx, facDict):
        #print(cls.__name__)
        return 1.0, None
    @classmethod
    def apply(s, state, descriptor):
        return state

def getPooledLOSCRV(facR):
    """
    Return the LOS CRV, pooled over CareTier and any other patient state, so that it
    depends only on facility abbrev
    """
    if 'losModel' in facR:
        return fullCRVFromLOSModel(facR['losModel'])
    else:
        raise RuntimeError('no losModel and meanLOS is not implemented yet')


class BinSwap(Mutator):
    """
    Move a fraction of the contents of one bin to another bin at the same facility
    """
    @classmethod
    def select(cls, state, facIdx, facDict):
        baseAR, baseDesc = super(BinSwap, cls).select(state, facIdx, facDict)
        facN = np.random.randint(0, state.shape[0])  # facility
        abbrev = facIdx[facN]
        dayN = np.random.randint(0, state.shape[1])  # time since arrival
        newDayN = np.random.randint(0, state.shape[1])  # new time since arrival
        fracMov = np.random.random()
        # fraction of weight to move
        facR = facDict[facIdx[facN]]
        if 'losCRV' not in facR:
            # losCRV represents the pooled LOS model over all conditions, including CareTier
            facR['losCRV'] = getPooledLOSCRV(facR)
        losCRV = facR['losCRV']
        # subState is meanPop * LOS PDF
        subState = np.sum(state[facN, :, :, :], axis=(1, 2))
        delta = fracMov * subState[dayN]
        #print("subState before: %s" % subState)
        subState[dayN] -= delta
        subState[newDayN] += delta
        #print("subState after: %s" % subState)
        #print("losCRV on dayN is %s" % losCRV.logpdf(float(dayN)))
        #print("losCRV on newDayN is %s" % losCRV.logpdf(float(newDayN)))
        beforeLnLik = (subState[dayN] * losCRV.logpdf(float(dayN) + 0.5)
                       + subState[newDayN] * losCRV.logpdf(float(newDayN) + 0.5))
        afterLnLik = ((subState[dayN] - delta) * losCRV.logpdf(float(dayN) + 0.5)
                       + (subState[newDayN] + delta) * losCRV.logpdf(float(newDayN) + 0.5))
        aR = min(afterLnLik/beforeLnLik, 1.0)
        #print('lnLiks: %s %s %s' % (beforeLnLik, afterLnLik, afterLnLik/beforeLnLik))
        return aR, (facN, dayN, newDayN, fracMov, baseDesc)
    @classmethod
    def apply(cls, state, descriptor):
        facN, dayN, newDayN, fracMov, baseDesc = descriptor
        state = super(BinSwap, cls).apply(state, baseDesc)
        delta = fracMov * state[facN, dayN, :, :]
        state[facN, dayN, :, :] -= delta
        state[facN, newDayN, :, :] += delta
        return state

class DirectTransfer(Mutator):
    pass

class IndirectTransfer(Mutator):
    pass

class NHIndirect(Mutator):
    pass

MUTATIONS = [BinSwap] #, DirectTransfer, IndirectTransfer, NHIndirect]

def main():
    """
    main
    """
    parser = OptionParser(usage="""
    %prog [--verbose] run_descr.yaml
    """)
    parser.add_option('-v', '--verbose', action='store_true',
                      help="request verbose output")

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

#     if not opts.notes:
#         parser.error('At least one --notes option is required')

    parser.destroy()

    configureLogging(getLoggerConfig())

    schemautils.setSchemaBasePath(SCHEMA_DIR)
    inputDict = checkInputFileSchema(args[0],
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    modelDir = inputDict['modelDir']
    pyrheautils.PATH_STRING_MAP['MODELDIR'] = os.path.abspath(modelDir)
    implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])
    pyrheautils.PATH_STRING_MAP['IMPLDIR'] = implDir

    facDirList = [pyrheautils.pathTranslate(fPth) for fPth in inputDict['facilityDirs']]
    facDict = readFacFiles(facDirList)
    for abbrev in EXCLUDE_THESE_FACILITIES:
        if abbrev in facDict:
            del facDict[abbrev]

    hospL = [abbrev for abbrev in facDict
             if facDict[abbrev]['category'] in ['HOSPITAL','LTAC', 'LTACH']]
    hospL.sort()
    hospIdx = {idx: abbrev for idx, abbrev in enumerate(hospL)}

    notHomeL = [abbrev for abbrev in facDict if facDict[abbrev]['category'] != 'COMMUNITY']
    notHomeL.sort()
    notHomeIdx = {idx: abbrev for idx, abbrev in enumerate(notHomeL)}

    notHospL = [abbrev for abbrev in facDict
                if facDict[abbrev]['category'] not in ['HOSPITAL','LTAC', 'LTACH']]
    notHospL.sort()
    notHospIdx = {idx: abbrev for idx, abbrev in enumerate(notHospL)}

    # Allocate and initialize the state
    nFac = len(notHomeIdx)
    nHosp = len(hospIdx)
    state = np.zeros((nFac, DAYS_PER_YEAR, DAYS_PER_YEAR, nHosp))
    cellsPerFac = 1
    for rng in state.shape[1:]:
        cellsPerFac *= rng
    for fN in xrange(nFac):
        abbrev = notHomeIdx[fN]
        if 'meanPop' in facDict[abbrev]:
            meanPop = facDict[abbrev]['meanPop']['value']
            state[fN,:,:,:] = float(meanPop)/cellsPerFac
        else:
            LOGGER.fatal('%s (%s) has no meanPop', abbrev, facDict[abbrev]['category'])
            raise RuntimeError('%s (%s) has no meanPop' % (abbrev, facDict[abbrev]['category']))

#    for loop in xrange(10):  # @UnusedVariable
    while True:
        mutator = np.random.choice(MUTATIONS)
        acceptRatio, descriptor = mutator.select(state, notHomeIdx, facDict)
        if np.random.random() < acceptRatio:
            print('aR %s -> update %s %s' % (acceptRatio, mutator.__name__, descriptor))
            state = mutator.apply(state, descriptor)
        else:
            print('aR %s -> reject %s' % (acceptRatio, mutator.__name__))

#############
# Input types:
# - lists by facility of fraction that return to some facility within a year (indirect readmit)
# - time to readmission line list and counts
# - and an infinite set of other combinations of stuff
#
# Consider a histogram, facilities x days at home, initially zero
# Consider N future hospital admissions on day 365 which are presumed to be readmissions
# 
# Every facility has a population, parameterized by
# (local LOS, time since last hosp, last hosp idx) (dim 365*365*31 = 4129975)
# Drop some large N of samples into this model and start doing random mutations:
# 1) A day passes, updating someone's local LOS and time since last hosp
# 2) A direct transfer occurs
# 3) An indirect transfer occurs
# 4) A nh-readmit occurs

#     with open(outFileName, 'w') as f:
#         yaml.safe_dump(sampleL, f, indent=4,
#                        encoding='utf-8', width=130, explicit_start=True)

if __name__ == "__main__":
    main()
