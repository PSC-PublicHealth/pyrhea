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
This script parses a file of arrival times to produce synthetic transfer matrix 
statistics for comparison with 'ground truth' transfer matrix statistics collected
from actual patient data.
"""

import sys
import os.path
cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))
import signal
import optparse
import yaml
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import logging
import logging.config
import yaml
import random
import cPickle as pickle

import schemautils
import pyrheautils
import tools_util as tu
from pyrhea import getLoggerConfig, checkInputFileSchema, loadPathogenImplementations
from pathogenbase import PthStatus
from typebase import PatientStatus, DiagClassA, PatientOverallHealth
from facilitybase import CareTier

SCHEMA_DIR = os.path.join(os.path.dirname(__file__), '../schemata')
INPUT_SCHEMA = 'rhea_input_schema.yaml'

DIAG_TIER_MAP = {DiagClassA.WELL: CareTier.HOME, DiagClassA.NEEDSREHAB: CareTier.NURSING,
                 DiagClassA.NEEDSLTAC: CareTier.LTAC, DiagClassA.SICK: CareTier.HOSP,
                 DiagClassA.VERYSICK: CareTier.ICU, DiagClassA.NEEDSVENT: CareTier.VENT,
                 DiagClassA.NEEDSSKILNRS: CareTier.SKILNRS}
TIER_DIAG_MAP = {v: k for k, v in DIAG_TIER_MAP.items()}

LOGGER = None

DEFAULT_LOW_DATE = 0
DEFAULT_HIGH_DATE = -1  # meaning include all records
DEFAULT_INPUT = 'ofile.pkl'

_DIAG_CLASS_A_MAP = {v: k for k, v in DiagClassA.names.items()}
_PATIENT_OVERALL_HEALTH_MAP = {v: k for k, v in PatientOverallHealth.names.items()}
_PTH_STATUS_MAP = {v: k for k, v in PthStatus.names.items()}


def tierFromCat(cat):
    return {'COMMUNITY': CareTier.HOME, 'SNF': CareTier.NURSING, 'LTACH': CareTier.LTAC,
            'VSNF': None, 'HOSPITAL': CareTier.HOSP, 'LTAC': CareTier.LTAC,
            'NURSINGHOME': CareTier.NURSING}[cat]


def mergeSelfTransfers(tupleL):
    """
    If two adjacent tuples represent a direct transfer to self, merge them.
    """
    rslt = []
    if tupleL:
        prevLoc, prevArrival, prevDeparture = tupleL[0]
        for loc, arrival, departure in tupleL[1:]:
            if loc != prevLoc or arrival != prevDeparture:
                rslt.append((prevLoc, prevArrival, prevDeparture))
                prevLoc, prevArrival, prevDeparture = loc, arrival, departure
            else:
                prevDeparture = departure
        rslt.append((prevLoc, prevArrival, prevDeparture))
    return rslt


def filterPath(tupleL, testFun, facDict, patName):
    rslt = [tpl for tpl in tupleL if testFun(tpl[0], facDict, patName)]
    return rslt


def overlapsRange(arrive, depart, low, high):
    return (arrive <= high and depart >= low)

def countVisits(tplL, rangeLow, rangeHigh, facDict, accumCountsDct, accumStayDct):
    assert rangeLow <= rangeHigh, 'cutoff dates are not in order'
    if tplL:
        dct = defaultdict(lambda: 0)
        sDct = defaultdict(lambda: 0.0)
        foundOverlap = False
        for loc, arrive, depart in tplL:
            if overlapsRange(arrive, depart, rangeLow, rangeHigh):
                if loc.endswith('cache'):
                    loc = loc[:-5]
                ctg = facDict[loc]['category']
                dct[ctg] += 1
                sDct[ctg] += (depart - arrive)
                foundOverlap = True
            else:
                if foundOverlap:
                    # We're out of time
                    break
        effectiveStart = max(rangeLow, tplL[0][1])
        effectiveEnd = min(rangeHigh, tplL[-1][2])
        for ctg, val in dct.items():
            accumCountsDct[ctg] += val
        for ctg, val in sDct.items():
            accumStayDct[ctg] += val
        return dict(dct), effectiveStart, effectiveEnd
    else:
        return {}, rangeLow, rangeHigh


class State(object):
    def __init__(self):
        self.patientD = {}
    def copy(self):
        state = State()
        state.patientD = self.patientD.copy()
        return state
    def add(self, patNm, pthStatus, origin, date, pLPSC, alive):
        if patNm in self.patientD:
            oldPthStatus, oldOrigin, oldDate, oldPLPSC, oldAlive = self.patientD[patNm]
            if (pthStatus != oldPthStatus or oldOrigin != origin or pLPSC != oldPLPSC
                or oldAlive != alive):
                self.patientD[patNm] = (pthStatus, origin, date, pLPSC, alive)
            else:
                pass
        else:
            self.patientD[patNm] = (pthStatus, origin, date, pLPSC, alive)
    def rm(self, patNm):
        if patNm in self.patientD:
            del self.patientD[patNm]
    def getCounts(self, facDict, reqPthStatus = None):
        rslt = defaultdict(int)
        if reqPthStatus is None:
            for patNm, (pthStatus, origin, date, pLPSC, alive) in self.patientD.items():
                idx = 'None' if pLPSC is None else facDict[pLPSC]['category']
                rslt[idx] += 1
        else:
            for patNm, (pthStatus, origin, date, pLPSC, alive) in self.patientD.items():
                if pthStatus == reqPthStatus:
                    idx = 'None' if pLPSC is None else facDict[pLPSC]['category']
                    rslt[idx] += 1
        return {key: val for key, val in rslt.items()}  # turn it into a regular dict
    def getTuple(self, patNm):
        pthStatus, origin, date, pLPSC, alive = self.patientD[patNm]
        return pthStatus, origin, date, pLPSC, alive
    def __str__(self):
        terms = []
        for patNm, (pthS, orig, date, pLPSC, alive) in self.patientD.items():
            pthStatusName = 'None' if (pthS is None) else PthStatus.names[pthS]
            terms.append('%s: (%s, %s, %s, %s)' % (patNm, pthStatusName, pLPSC,
                                                   ('alive' if alive else 'died'),
                                                   date))
        dctS = ', '.join(terms)
        return 'State(%s)' % dctS
    def __add__(self, other):
        rslt = self.copy()
        for patNm, (pthStatus, origin, date, pLPSC, alive) in other.patientD.items():
            rslt.add(patNm, pthStatus, origin, date, pLPSC, alive)
        return rslt
    def __iadd__(self, other):
        for patNm, (pthStatus, origin, date, pLPSC, alive) in other.patientD.items():
            self.add(patNm, pthStatus, origin, date, pLPSC, alive)
        return self


DEAD_LIST = []
LAST_TOD = -1

def buildStateChain(placeEvtD, targetLoc, targetTier, facDict):
    global DEAD_LIST
    global LAST_TOD
    try:
        eventL = placeEvtD[(targetLoc, targetTier)] + placeEvtD[(targetLoc, None)]
        eventL.sort()
    except KeyError:
        eventL = []
    day = 0
    stateD = {day: State()}
    while eventL:
        nextEvent = eventL.pop(0)
        nextEDay = nextEvent[0]
        assert day <= nextEDay, 'Events out of order; day %d vs expected %d' % (day, nextEDay)
        while day < nextEDay:
            stateD[day + 1] = stateD[day].copy()
            day += 1
        stateNow = stateD[day]
        if len(nextEvent) == 6:
            # update, possible departure
            #print nextEvent
            patNm, loc, patStatus, pLPSC, alive = nextEvent[1:]  # pLPSC is 'place of last pth status change'
            if patNm in DEAD_LIST:
                raise RuntimeError('zombie! %s' % str(nextEvent))
            pthStatus = None if patStatus is None else patStatus.pthStatus
            if not alive:
                print 'ping %s' % patNm
                stateNow.rm(patNm)
                DEAD_LIST.append(patNm)
                LAST_TOD = max(LAST_TOD, day)
            elif (day != 0 and patStatus not in [None, targetTier]) or loc != targetLoc:
                # departure
                oldPS, oldOrigin, oldDate, oldPLPSC, oldAlive = stateNow.getTuple(patNm)
                assert oldAlive, 'found %s dead at day %s' % (patNm, oldAlive)
                if oldPS != pthStatus:
                    newPthStatus = pthStatus
                    backD = (day + oldDate)/2
                    for tD in range(backD, day):
                        stateD[tD].add(patNm, newPthStatus, -1, backD, pLPSC, alive)  # -1 being the code for local
                stateNow.rm(patNm)
            else:
                # time=0 creation events
                stateNow.add(patNm, pthStatus, targetTier, day, pLPSC, alive)
        elif len(nextEvent) == 7:
            # arrival
            #print 'point 2 ', nextEvent
            patNm, loc, patStatus, srcAddr, pLPSC, alive = nextEvent[1:]
            if patNm in DEAD_LIST:
                raise RuntimeError('zombie! (point 2) %s' % str(nextEvent))
            assert alive, '%s arrived at %s dead on day %s' % (patNm, loc, day)
            pthStatus = None if patStatus is None else patStatus.pthStatus
            srcLoc, srcTier = srcAddr
            if srcTier is None:
                srcTier = tierFromCat(facDict[srcLoc]['category'])
            stateNow.add(patNm, pthStatus, srcTier, day, pLPSC, alive)
        else:
            raise RuntimeError('Unexpected event format %s' % str(nextEvent))
#         print nextEvent
#         print day
#         print sum(stateD[day].getCounts(facDict).values()), stateD[day].getCounts(facDict)
    for nm in DEAD_LIST:
        if nm in stateD:
            raise RuntimeError('Vampire!')
    return stateD, day


def main():
    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)

    global DEAD_LIST

    global LOGGER
    logging.config.dictConfig(getLoggerConfig())
    LOGGER = logging.getLogger(__name__)

    parser = optparse.OptionParser(usage="""
    %prog [-L low_date] [-H high_date] -i infile.pkl run_descr.yaml
    """)
    parser.add_option('-L', '--low', action='store', type='int',
                      help='minimum date to include',
                      default=DEFAULT_LOW_DATE)
    parser.add_option('-H', '--high', action='store', type='int',
                      help='maximum date to include',
                      default=DEFAULT_HIGH_DATE)
    parser.add_option('-i', '--input', action='store', type='string',
                      help='file containing output of arrival_stream_parser',
                      default=DEFAULT_INPUT)
    parser.add_option('-f', '--facility', action='store', type='string',
                      help='restrict plot to this facility', default=None)

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('An input yaml file matching %s is required' % INPUT_SCHEMA)

    parser.destroy()
    lowDate = opts.low
    highDate = opts.high

    if opts.facility:
        targetFac = opts.facility

    inputDict = tu.readModelInputs(args[0])
    facDict = tu.getFacDict(inputDict)

    infile = opts.input 
    with open(infile, 'rU') as f:
        placeEvtD = pickle.load(f)

    if targetFac:
        cat = facDict[targetFac]['category']
        tier = tierFromCat(cat)
        if tier is None:
            sys.exit('There is no single tier of care for %s')
    else:
        #cat, tier = 'NURSINGHOME', CareTier.NURSING
        cat, tier = 'LTAC', CareTier.LTAC
        #cat, tier = 'HOSPITAL', CareTier.HOSP
        #targetFac = 'PLUM_24_S'
        #targetFac = 'ALDE_6120_S'
        #targetFac = 'EVAN_1300_S'
        #targetFac = 'RML_5601_L'    
        #targetFac = 'PRES_100_L'    
        #targetFac = 'THC_365_L'    
        

    allMergeD = defaultdict(lambda: defaultdict(int))
    pthMergeD = defaultdict(lambda: defaultdict(int))
    allKeySet = set()
    maxMaxDay = 0
    if targetFac is None:
        for fac, rec in facDict.items():
            if rec['category'] == cat:
                DEAD_LIST = []
                stateD, maxDay = buildStateChain(placeEvtD, fac, tier, facDict)
                #print fac,stateD[0]
                for day in xrange(maxDay+1):
                    allCounts = stateD[day].getCounts(facDict)
                    pthCounts = stateD[day].getCounts(facDict, PthStatus.COLONIZED)
                    for key, val in allCounts.items():
                        allMergeD[day][key] += val
                    for key, val in pthCounts.items():
                        pthMergeD[day][key] += val
                    allKeySet.update(allMergeD[day])
                maxMaxDay = max(maxMaxDay, maxDay)
            #print fac
            #break
    else:
        fac = targetFac
        rec = facDict[fac]
        assert rec['category'] == cat, 'Target facility %s is not in category %s' % (targetFac, cat)
        stateD, maxDay = buildStateChain(placeEvtD, fac, tier, facDict)
        #print fac,stateD[0]
        for day in xrange(maxDay+1):
            allCounts = stateD[day].getCounts(facDict)
            pthCounts = stateD[day].getCounts(facDict, PthStatus.COLONIZED)
            for key, val in allCounts.items():
                allMergeD[day][key] += val
            for key, val in pthCounts.items():
                pthMergeD[day][key] += val
            allKeySet.update(allMergeD[day])
        maxMaxDay = maxDay

#     for day in xrange(maxMaxDay+1):
#         allD = {key: val for key, val in allMergeD[day].items()}
#         pthD = {key: val for key, val in pthMergeD[day].items()}
#         print 'Day %s: %s %s: %s %s' % (day, sum(allD.values()),
#                                         float(sum(pthD.values()))/float(1+sum(allD.values())),
#                                         allD, pthD)
#         print '-----------'
    print 'Total dead: %d' % len(DEAD_LIST)
    print 'Last time of death: %d' % LAST_TOD

    keys = list(allKeySet)
    keys.sort()
    xV = np.linspace(0.0, maxMaxDay, num=maxMaxDay+1)
    yPthA = np.zeros((len(keys), maxMaxDay+1))
    yAllA = np.zeros((len(keys), maxMaxDay+1))
    for col in xrange(maxMaxDay+1):
        for row, key in enumerate(keys):
            yPthA[row, col] = pthMergeD[col][key]
            yAllA[row, col] = allMergeD[col][key]
    fig, axes = plt.subplots(2, 2)
    axes[0,0].stackplot(xV, yAllA, labels = keys)
    #axes[0,0].legend()
    axes[0,0].set_title('Sources of All Patients')
    axes[0,1].stackplot(xV, yPthA, labels = keys)
    handles, labels = axes[0, 1].get_legend_handles_labels()
    #axes[0,1].legend()
    axes[0,1].set_title('Sources of Colonization')
    ratioA = yPthA / yAllA
    for row, key in enumerate(keys):
        axes[1,0].plot(xV, ratioA[row, :], '-', label=key)
    axes[1,0].set_title('Relative Fractions')
    axes[1,1].remove()
    fig.legend(handles, labels, 'lower right')
    if targetFac is None:
        fig.suptitle('%s %s %s' % (infile, cat, CareTier.names[tier]))
    else:
        fig.suptitle('%s %s %s' % (infile, targetFac, CareTier.names[tier]))
    plt.show()

if __name__ == "__main__":
    main()
