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

###########################
# Notes-
# -What are we doing with self-transfers? What should we be doing?
# -What does it mean if a FacilitySwap causes srcN to be a place that is never a src
#  for the new facility?  Does it even matter?  Should all FacilitySwaps set src to N_FAC?
# -We need a burn-in, and maybe to save samples only every few mutations
###########################

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
N_FAC = None  # Total number of facilities
N_SRC = None  # Possible places to come from - N_FAC + 1 for 'home'
N_HOSP = None  # Number that count as hospitals for indirect-transfer purposes
N_PREV_HOSP = None  # Possible most recent hospitals - N_HOSP + 1 for 'never hospitalized'
DIRECT_TRANSFER_TABLE = None  # This is indexed like DTT[thisLoc][whereTransferCameFrom]

def configureLogging(cfgDict):
    global LOGGER
    logging.config.dictConfig(cfgDict)
    LOGGER = logging.getLogger(__name__)

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
    return facR['_meanPooledLOSCRV']


class WeightedListPicker(object):
    def __init__(self, wD):
        """
        wD is a dict of the form {key: weight, ...} .  All the weights are summed and a
        random number is drawn and scaled by the total weight.  The key corresponding to
        the resulting fraction of the total weight is returned.  This implementation assumes
        the keys are integers.
        """
        pairL = [(float(wt), key) for key, wt in wD.items()]
        pairL.sort()
        tot = sum([wt for wt, key in pairL])
        runTot = 0.0
        intervalL = []
        for wt, key in pairL:
            runTot += wt
            intervalL.append((runTot/tot, key))
        self.breakV = np.asarray([rwt for rwt, key in intervalL])
        self.idxV = np.asarray([key for rwt, key in intervalL])
    def draw(self):
        offset = np.searchsorted(self.breakV, np.random.random())
        return self.idxV[offset]

def testWeightedListPicker():
    sampD = {i:np.random.random() for i in xrange(10)}
    print(sampD)
    wLP = WeightedListPicker(sampD)
    counter = np.zeros([10], dtype=np.int32)
    for i in xrange(1000000):
        counter[wLP.draw()] += 1
    tot = sum([val for idx, val in sampD.items()])
    for i in xrange(10):
        print('%s: %s vs %s' % (i, sampD[i]/tot, float(counter[i])/1000000.))


def getRandomDirectTransferSrc(facN, facR):
    if '_directTransferWLP' not in facR:
        wlp = WeightedListPicker(DIRECT_TRANSFER_TABLE[facN])
        facR['_directTransferWLP'] = wlp
    return facR['_directTransferWLP'].draw()

class SampleCollection(object):
    """
    This is implemented as a class to allow more efficient representations.
    """
    def __init__(self):
        self.samples = np.zeros((N_FAC, DAYS_PER_YEAR, DAYS_PER_YEAR, N_SRC),
                                dtype=np.int32)

    @classmethod
    def state(cls, facN, dayN, hospDayN, srcN, phN):
        """
        Arguments are facility number, days at that facility, days since last hospitalization,
        previous facility number, most recent hospital number
        """
        return (facN, dayN, hospDayN, srcN, phN)

    @classmethod
    def split(cls, state):
        """
        The opposite of state()
        """
        # This works as long as state is just a tuple
        return state

    def sample(self, state):
        facN, dayN, hospDayN, srcN, phN = state
        self.samples[facN, dayN, hospDayN, srcN] += 1


class Mutator(object):
    @classmethod
    def select(cls, state, facIdx, facDict):
        #print(cls.__name__)
        return 1.0, None
    @classmethod
    def apply(s, state, descriptor):
        return state


class BinSwap(Mutator):
    """
    Change LOS date, staying within a facility
    """
    @classmethod
    def select(cls, state, facIdx, facDict):
        facN, dayN, hospDay, srcN, phN = SampleCollection.split(state)
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
        facN, dayN, hospDay, srcN, phN = SampleCollection.split(state)
        return SampleCollection.state(facN, newDayN, hospDay, srcN, phN)


class FacSwap(Mutator):
    """
    Change facilities. Acceptance ratio is based on relative mean population, so there is no
    correlation between before-and-after dates- that would require a conditional probability
    """
    @classmethod
    def select(cls, state, facIdx, facDict):
        facN, dayN, hospDay, srcN, phN = SampleCollection.split(state)
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
        facN, dayN, hospDay, srcN, phN = SampleCollection.split(state)
        return SampleCollection.state(newFacN, newDayN, hospDay, srcN, phN)


class IsTransferSwap(Mutator):
    """
    Within a facility, change whether or not the sample was a direct transfer in.  Acceptance
    ratio is based on the fraction of patients that are direct transfers.
    """
    @classmethod
    def select(cls, state, facIdx, facDict):
        facN, dayN, hospDay, srcN, phN = SampleCollection.split(state)
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
        if srcN == N_FAC:
            return baseAR * transferFrac, (transferFrac, getRandomDirectTransferSrc(facN, facR),
                                           baseDesc)
        else:
            return baseAR * (1.0 - transferFrac), (transferFrac, N_FAC, baseDesc)

    @classmethod
    def apply(cls, state, descriptor):
        transferFrac, newSrcN, baseDesc = descriptor
        state = super(IsTransferSwap, cls).apply(state, baseDesc)
        facN, dayN, hospDay, srcN, phN = SampleCollection.split(state)
        return SampleCollection.state(facN, dayN, hospDay, newSrcN, phN)


class DirectTransferSrcSwap(Mutator):
    """
    Within a facility, if the sample has srcN that points to another facility,
    swap that srcN for a different srcN.
    """
    @classmethod
    def select(cls, state, facIdx, facDict):
        facN, dayN, hospDay, srcN, phN = SampleCollection.split(state)
        baseAR, baseDesc = super(DirectTransferSrcSwap, cls).select(state, facIdx, facDict)
        if srcN == N_FAC:
            # No direct transfer src swaps if the src was the community
            return 0.0, (srcN, baseDesc)
        myDDT = DIRECT_TRANSFER_TABLE[facN]
        newSrcN = np.random.choice(myDDT.keys())
        if newSrcN == srcN:
            return 0.0, (srcN, baseDesc)  # why bother
        nWt = float(myDDT[newSrcN])
        if srcN in myDDT:
            #oWt = (float(myDDT[srcN]) if srcN in myDDT else 0.0)
            oWt = float(myDDT[srcN])
            return baseAR * min(1.0, nWt/oWt), (newSrcN, baseDesc)
        else:
            return baseAR, (newSrcN, baseDesc)

    @classmethod
    def apply(cls, state, descriptor):
        newSrcN, baseDesc = descriptor
        state = super(DirectTransferSrcSwap, cls).apply(state, baseDesc)
        facN, dayN, hospDay, srcN, phN = SampleCollection.split(state)
        return SampleCollection.state(facN, dayN, hospDay, newSrcN, phN)


class IsIndirectTransferSwap(Mutator):
    """
    To turn a direct transfer into an indirect transfer:
    -srcN must be a hospital
    -facN must be a hospital
    -srcN becomes phN, and srcN becomes 'home'
    -dayN has to get split between dayN and hospDay
    Factors contributing to relative odds:
     -current odds: LOS distro with dayN, odds of transfer from srcN to facN
     -new odds: home LOS distro for new hospDay, local LOS distro for remaining dayN,
      odds of transfer from new phN to home, odds of transfer home to facN
    Constraints on state:
     -srcN must be home
    """
    @classmethod
    def select(cls, state, facIdx, facDict):
        facN, dayN, hospDay, srcN, phN = SampleCollection.split(state)
        baseAR, baseDesc = super(DirectTransferSrcSwap, cls).select(state, facIdx, facDict)
        if srcN == N_FAC:
            # No direct transfer src swaps if the src was the community
            return 0.0, (srcN, baseDesc)
        myDDT = DIRECT_TRANSFER_TABLE[facN]
        newSrcN = np.random.choice(myDDT.keys())
        if newSrcN == srcN:
            return 0.0, (srcN, baseDesc)  # why bother
        nWt = float(myDDT[newSrcN])
        if srcN in myDDT:
            #oWt = (float(myDDT[srcN]) if srcN in myDDT else 0.0)
            oWt = float(myDDT[srcN])
            return baseAR * min(1.0, nWt/oWt), (newSrcN, baseDesc)
        else:
            return baseAR, (newSrcN, baseDesc)

    @classmethod
    def apply(cls, state, descriptor):
        newSrcN, baseDesc = descriptor
        state = super(DirectTransferSrcSwap, cls).apply(state, baseDesc)
        facN, dayN, hospDay, srcN, phN = SampleCollection.split(state)
        return SampleCollection.state(facN, dayN, hospDay, newSrcN, phN)


class NHIndirect(Mutator):
    pass

MUTATIONS = [BinSwap, FacSwap, IsTransferSwap, DirectTransferSrcSwap]

def main():
    """
    main
    """

    global N_FAC
    global N_SRC
    global N_HOSP
    global N_PREV_HOSP
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

    hospL = [abbrev for abbrev in facDict
             if facDict[abbrev]['category'] in ['HOSPITAL','LTAC', 'LTACH']]
    hospL.sort()
    hospIdx = {idx: abbrev for idx, abbrev in enumerate(hospL)}

    notHomeL = [abbrev for abbrev in facDict if facDict[abbrev]['category'] != 'COMMUNITY']
    notHomeL.sort()
    facIdx = {idx: abbrev for idx, abbrev in enumerate(notHomeL)}
    invFacIdx = {abbrev: idx for idx, abbrev in enumerate(notHomeL)}

    notHospL = [abbrev for abbrev in facDict
                if facDict[abbrev]['category'] not in ['HOSPITAL','LTAC', 'LTACH']]
    notHospL.sort()
    notHospIdx = {idx: abbrev for idx, abbrev in enumerate(notHospL)}

    # Load and invert the direct transfer table
    directTransferTblPath = pyrheautils.pathTranslate('$(MODELDIR)/transfer_counts.yaml')
    with open(directTransferTblPath, 'rU') as f:
        transferDict = yaml.load(f)
    invertDict = {}
    for src, rec in transferDict.items():
        totTransfersOut = float(sum([val for k, val in rec.items()]))
        if totTransfersOut == 0.0:
            print('%s has no outgoing transfers' % src)
        for dst, val in rec.items():
            dstN = invFacIdx[dst]
            if dstN not in invertDict:
                invertDict[dstN] = {}
            invertDict[dstN][invFacIdx[src]] = val
    DIRECT_TRANSFER_TABLE = invertDict

    # Allocate and initialize the state
    N_FAC = len(facIdx)
    N_SRC = N_FAC + 1  # the additional value is for 'came from home'
    N_HOSP = len(hospIdx)
    N_PREV_HOSP = N_HOSP + 1 # the additional value is 'no previous hospitalization'
    samples = SampleCollection()
    facN = 1
    losDay = 0
    hospDay = 1
    srcN = N_FAC  # which means this sample came from 'home'
    prevHospN = N_HOSP  # which means no prev hospitalization
    state = samples.state(facN, losDay, hospDay, srcN, prevHospN)
    for loop in xrange(100000):  # @UnusedVariable
#    while True:
        samples.sample(state)
        mutator = np.random.choice(MUTATIONS)
        acceptRatio, descriptor = mutator.select(state, facIdx, facDict)
        if np.random.random() < acceptRatio:
            print('aR %s -> update %s %s' % (acceptRatio, mutator.__name__, descriptor))
            state = mutator.apply(state, descriptor)
        else:
            print('aR %s -> reject %s %s' % (acceptRatio, mutator.__name__, descriptor))
    facN, losDay, hospDay, srcN, prevHospN = SampleCollection.split(state)
    print(np.sum(samples.samples, 0)[:, hospDay, srcN])
    print(np.sum(samples.samples, 1)[:, hospDay, srcN])

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
