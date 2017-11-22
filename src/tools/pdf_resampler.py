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

"""
Model days per year
"""
DAYS_PER_YEAR = 365  # But might not, in HERMES for example
#DAYS_PER_YEAR = 10  # But might not, in HERMES for example

CENSUS_DAYS = 365 # The ground truth data really was gathered for a 1 year interval

"""
Numbers of facilities, initialized when the model is loaded
"""
N_FAC = None
N_HOSP = None
DIRECT_TRANSFER_TABLE = None

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
    if '_pooledLOSCRV' not in facR:
        if 'losModel' in facR:
            facR['_pooledLOSCRV'] = fullCRVFromLOSModel(facR['losModel'])
        else:
            raise RuntimeError('no losModel and meanLOS is not implemented yet')
    return facR['_pooledLOSCRV']

def getMeanPooledLOSCRV(facR):
    """
    Return the mean of the pooled LOS distribution, with caching.
    """
    if '_meanPooledLOSCRV' not in facR:
        facR['_meanPooledLOSCRV'] = getPooledLOSCRV(facR).mean()
        print('getMeanPooledLOSCRV %s' % facR['abbrev'])
    return facR['_meanPooledLOSCRV']


class BinSwap(Mutator):
    """
    Change LOS date, staying within a facility
    """
    @classmethod
    def select(cls, state, facIdx, facDict):
        facN, dayN, hospDay, hospIdx, isTransfer = state
        baseAR, baseDesc = super(BinSwap, cls).select(state, facIdx, facDict)
        facR = facDict[facIdx[facN]]
        newDayN = np.random.randint(0, DAYS_PER_YEAR)  # new time since arrival
        losCRV = getPooledLOSCRV(facR)
        aR = baseAR * min(losCRV.pdf(float(newDayN) + 0.5)/losCRV.pdf(float(dayN) + 0.5),
                          1.0)
        return aR, (newDayN, baseDesc)
    @classmethod
    def apply(cls, state, descriptor):
        newDayN, baseDesc = descriptor
        state = super(BinSwap, cls).apply(state, baseDesc)
        facN, dayN, hospDay, hospIdx, isTransfer = state
        return (facN, newDayN, hospDay, hospIdx, isTransfer)


class FacSwap(Mutator):
    """
    Change facilities. Acceptance ratio is based on relative mean population, so there is no
    correlation between before-and-after dates- that would require a conditional probability
    """
    @classmethod
    def select(cls, state, facIdx, facDict):
        facN, dayN, hospDay, hospIdx, isTransfer = state
        baseAR, baseDesc = super(FacSwap, cls).select(state, facIdx, facDict)
        facR = facDict[facIdx[facN]]
        if 'meanPop' in facR:
            oMP = facR['meanPop']['value']
        else:
            raise RuntimeError('%s has no meanPop' % facR['abbrev'])
        newFacN = np.random.randint(0, N_FAC)  # new facility
        if newFacN == facN:
            # reject self-transfers because the mutator would conflict with the LOS pdf
            return 0.0, baseDesc
        newFacR = facDict[facIdx[newFacN]]
        if 'meanPop' in newFacR:
            nMP = newFacR['meanPop']['value']
        else:
            raise RuntimeError('%s has no meanPop' % newFacR['abbrev'])
        aR = baseAR * min((nMP/oMP), 1.0)
        return aR, (newFacN, baseDesc)
    @classmethod
    def apply(cls, state, descriptor):
        newFacN, baseDesc = descriptor
        newDayN = np.random.randint(0, DAYS_PER_YEAR)
        state = super(FacSwap, cls).apply(state, baseDesc)
        facN, dayN, hospDay, hospIdx, isTransfer = state
        return (newFacN, newDayN, hospDay, hospIdx, isTransfer)


class IsTransferSwap(Mutator):
    """
    Within a facility, change whether or not the sample was a direct transfer in.  Acceptance
    ratio is based on the fraction of patients that are direct transfers.
    """
    @classmethod
    def select(cls, state, facIdx, facDict):
        facN, dayN, hospDay, hospIdx, isTransfer = state
        baseAR, baseDesc = super(IsTransferSwap, cls).select(state, facIdx, facDict)
        facR = facDict[facIdx[facN]]
        if 'meanPop' not in facR:
            raise RuntimeError('%s has no meanPop' % facR['abbrev'])
        meanPop = facR['meanPop']['value']
        if 'totalTransfersIn' not in facR:
            raise RuntimeError('%s has no totalTransfersIn' % facR['abbrev'])
        ttI = facR['totalTransfersIn']['value']
        meanPooledLOS = getMeanPooledLOSCRV(facR)
        transferPop = ttI*(meanPooledLOS/CENSUS_DAYS)
        transferFrac = min(1.0, transferPop/meanPop)
        if isTransfer:
            return baseAR * (1.0 - transferFrac), (transferFrac, baseDesc)
        else:
            return baseAR * transferFrac, (transferFrac, baseDesc)
    @classmethod
    def apply(cls, state, descriptor):
        transferFrac, baseDesc = descriptor
        state = super(IsTransferSwap, cls).apply(state, baseDesc)
        facN, dayN, hospDay, hospIdx, isTransfer = state
        isTransfer = (0 if isTransfer else 1)
        return (facN, dayN, hospDay, hospIdx, isTransfer)



class DirectTransfer(Mutator):
    pass

class IndirectTransfer(Mutator):
    pass

class NHIndirect(Mutator):
    pass

MUTATIONS = [BinSwap, FacSwap, IsTransferSwap] #, DirectTransfer, IndirectTransfer, NHIndirect]

def main():
    """
    main
    """

    global N_FAC
    global N_HOSP
    global DIRECT_TRANSFER_TABLE
    
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

#     directTransferTblPath = pyrheautils.pathTranslate('$(MODELDIR)/transfer_counts.yaml')
#     with open(directTransferTblPath, 'rU') as f:
#         transferDict = yaml.load(f)
#     invertDict = {}
#     for src, r in transferDict.items():
#         totTransfersOut = float(sum([v for k, v in r.items()]))
#         if totTransfersOut == 0.0:
#             print('%s has no outgoing transfers' % src)
#         for dst, v in r.items():
#             if dst not in invertDict:
#                 invertDict[dst] = {}
#             invertDict[dst][src] = v
#     DIRECT_TRANSFER_TABLE = transferDict
#     print(DIRECT_TRANSFER_TABLE.keys())
#     print(DIRECT_TRANSFER_TABLE['SJUD'])

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
    N_FAC = len(notHomeIdx)
    N_HOSP = len(hospIdx)
    samples = np.zeros((N_FAC, DAYS_PER_YEAR, DAYS_PER_YEAR, N_HOSP, 2), dtype=np.int32)
    facIdx = 1
    losDay = 0
    hospDay = 1
    hospIdx = 1
    isTransfer = 0  # or 1
    state = (facIdx, losDay, hospDay, hospIdx, 0)
    for loop in xrange(100000):  # @UnusedVariable
#    while True:
        facIdx, losDay, hospDay, hospIdx, isTransfer = state
        samples[facIdx, losDay, hospDay, hospIdx, isTransfer] += 1
        mutator = np.random.choice(MUTATIONS)
        acceptRatio, descriptor = mutator.select(state, notHomeIdx, facDict)
        if np.random.random() < acceptRatio:
            print('aR %s -> update %s %s' % (acceptRatio, mutator.__name__, descriptor))
            state = mutator.apply(state, descriptor)
        else:
            print('aR %s -> reject %s %s' % (acceptRatio, mutator.__name__, descriptor))
    facIdx, losDay, hospDay, hospIdx, isTransfer = state
    print(np.sum(samples, 0)[:, hospDay, hospIdx, isTransfer])
    print(np.sum(samples, 1)[:, hospDay, hospIdx, isTransfer])
    print(np.sum(samples, (0,1))[hospDay, hospIdx, :])

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
