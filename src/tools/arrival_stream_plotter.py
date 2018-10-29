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
    return {'COMMUNITY': CareTier.HOME, 'SNF': None, 'LTACH': CareTier.LTAC,
            'VSNF': None, 'HOSPITAL': CareTier.HOSP}[cat]


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
    def add(self, patNm, pthStatus, origin, date):
        if patNm in self.patientD:
            oldPthStatus, oldOrigin, oldDate = self.patientD[patNm]
            if pthStatus != oldPthStatus or oldOrigin != origin:
                self.patientD[patNm] = (pthStatus, origin, date)
            else:
                pass
        else:
            self.patientD[patNm] = (pthStatus, origin, date)
    def rm(self, patNm):
        if patNm in self.patientD:
            del self.patientD[patNm]
    def originStr(self, origin):
        if origin < 0:
            return 'local'
        else:
            return CareTier.names[origin]
    def getCounts(self, reqPthStatus = None):
        rslt = defaultdict(int)
        if reqPthStatus is None:
            for patNm, (pthStatus, origin, date) in self.patientD.items():
                rslt[self.originStr(origin)] += 1
        else:
            for patNm, (pthStatus, origin, date) in self.patientD.items():
                if pthStatus == reqPthStatus:
                    rslt[self.originStr(origin)] += 1
        return {key: val for key, val in rslt.items()}  # turn it into a regular dict
    def getTuple(self, patNm):
        pthStatus, origin, date = self.patientD[patNm]
        return pthStatus, origin, date
    def __str__(self):
        terms = []
        for patNm, (pthS, orig, date) in self.patientD.items():
            pthStatusName = 'None' if (pthS is None) else PthStatus.names[pthS]
            if orig is None:
                careTierName = 'tier_None'
            elif orig == -1:
                careTierName = 'local'
            else:
                careTierName = CareTier.names[orig]
            terms.append('%s: (%s, %s, %s)' % (patNm, pthStatusName, careTierName,
                                               date))
        dctS = ', '.join(terms)
        return 'State(%s)' % dctS
    def __add__(self, other):
        rslt = self.copy()
        for patNm, (pthStatus, origin, date) in other.patientD.items():
            rslt.add(patNm, pthStatus, origin, date)
        return rslt
    def __iadd__(self, other):
        for patNm, (pthStatus, origin, date) in other.patientD.items():
            self.add(patNm, pthStatus, origin, date)
        return self


def buildStateChain(placeEvtD, targetLoc, targetTier, facDict):
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
        #print 'nextEvent: %s' % str(nextEvent)
        assert day <= nextEDay, 'Events out of order; day %d vs expected %d' % (day, nextEDay)
        while day < nextEDay:
            stateD[day + 1] = stateD[day].copy()
            day += 1
        stateNow = stateD[day]
        if len(nextEvent) == 4:
            # update, possible departure
            patNm, loc, patStatus = nextEvent[1:]
            pthStatus = None if patStatus is None else patStatus.pthStatus
            if (day != 0 and patStatus not in [None, targetTier]) or loc != targetLoc:
                # departure
                oldPS, oldOrigin, oldDate = stateNow.getTuple(patNm)
                if oldPS != pthStatus:
                    newPthStatus = pthStatus
                    backD = (day + oldDate)/2
                    for tD in range(backD, day):
                        stateD[tD].add(patNm, newPthStatus, -1, backD)  # -1 being the code for local
                stateNow.rm(patNm)
            else:
                # time=0 creation events
                stateNow.add(patNm, pthStatus, targetTier, day)
        elif len(nextEvent) == 5:
            # arrival
            patNm, loc, patStatus, srcAddr = nextEvent[1:]
            pthStatus = None if patStatus is None else patStatus.pthStatus
            srcLoc, srcTier = srcAddr
            if srcTier is None:
                srcTier = tierFromCat(facDict[srcLoc]['category'])
            stateNow.add(patNm, pthStatus, srcTier, day)
        else:
            raise RuntimeError('Unexpected event format %s' % nextEvent)
        #print nextEvent
        #print day
        #print sum(stateD[day].getCounts().values()), stateD[day].getCounts()
    return stateD, day


def main():
    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)

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

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('An input yaml file matching %s is required' % INPUT_SCHEMA)

    parser.destroy()
    lowDate = opts.low
    highDate = opts.high

    inputDict = tu.readModelInputs(args[0])
    facDict = tu.getFacDict(inputDict)

    infile = opts.input 
    with open(infile, 'rU') as f:
        placeEvtD = pickle.load(f)

    cat, tier = 'SNF', CareTier.NURSING
    #targetFac = None
    #targetFac = 'PLUM_24_S'
    #targetFac = 'ALDE_6120_S'
    targetFac = 'EVAN_1300_S'    

    allMergeD = defaultdict(lambda: defaultdict(int))
    pthMergeD = defaultdict(lambda: defaultdict(int))
    maxMaxDay = 0
    if targetFac is None:
        for fac, rec in facDict.items():
            if rec['category'] == cat:
                stateD, maxDay = buildStateChain(placeEvtD, fac, tier, facDict)
                #print fac,stateD[0]
                for day in xrange(maxDay+1):
                    allCounts = stateD[day].getCounts()
                    pthCounts = stateD[day].getCounts(PthStatus.COLONIZED)
                    for key, val in allCounts.items():
                        allMergeD[day][key] += val
                    for key, val in pthCounts.items():
                        pthMergeD[day][key] += val
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
            allCounts = stateD[day].getCounts()
            pthCounts = stateD[day].getCounts(PthStatus.COLONIZED)
            for key, val in allCounts.items():
                allMergeD[day][key] += val
            for key, val in pthCounts.items():
                pthMergeD[day][key] += val
        maxMaxDay = maxDay

    for day in xrange(maxMaxDay+1):
        allD = {key: val for key, val in allMergeD[day].items()}
        pthD = {key: val for key, val in pthMergeD[day].items()}
        print 'Day %s: %s %s: %s %s' % (day, sum(allD.values()), float(sum(pthD.values()))/float(1+sum(allD.values())),
                                        allD, pthD)
        print '-----------'

    keys = CareTier.names.values() + ['local']
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
