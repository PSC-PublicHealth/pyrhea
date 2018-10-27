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

LOGGER = None

DEFAULT_LOW_DATE = 0
DEFAULT_HIGH_DATE = -1  # meaning include all records

_DIAG_CLASS_A_MAP = {v: k for k, v in DiagClassA.names.items()}
_PATIENT_OVERALL_HEALTH_MAP = {v: k for k, v in PatientOverallHealth.names.items()}
_PTH_STATUS_MAP = {v: k for k, v in PthStatus.names.items()}


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


class ParseLineError(RuntimeError):
    pass

def parseArrivalLine(line):
    if 'birth' in line: print line
    words =  line.split()
    while words[0].strip() == '(Pdb)':
        words = words[1:]
    if len(words) == 6:
        assert words[3] == 'arrives', 'Bad line format: %s' % line
        patName = words[2]
        dstName = words[4]
        date = int(words[5])
    elif len(words) == 7:
        assert words[4] == 'arrives', 'Bad line format: %s' % line
        assert words[2] == '-', 'Bad line format: %s' % line
        patName = words[3]
        dstName = words[5]
        date = int(words[6])
    else:
        raise ParseLineError('Bad line format: %s' % line)
    return patName, dstName, date

def parseStatusLine(line, patName):
    words = line.split()
    while words[0].strip() == '(Pdb)':
        words = words[1:]
    assert words[3] == patName, 'Bad line format for patient %s: %s' % (patName, line)
    assert words[6].startswith('PatientStatus('), 'Bad line format: %s' % line
    statS = ' '.join(words[6:])
    statS = statS.strip()[14:-1]
    terms = statS.split(',')
    kwargs = {'homeAddr': None}
    for term in terms:
        term = term.strip()
        if term.startswith('homeAddr'):
            continue  # we cannot reconstruct these
        try:
            key, val = term.split('=')
        except ValueError:
            continue  # fragments
        key = key.strip()
        if key == 'overall':
            val = _PATIENT_OVERALL_HEALTH_MAP[val]
        elif key == 'diagClassA':
            val = _DIAG_CLASS_A_MAP[val]
        elif key == 'pthStatus':
            val = _PTH_STATUS_MAP[val]
        elif key in ['startDateA', 'startDatePth']:
            val = int(val)
        elif key in ['relocateFlag', 'justArrived', 'canClear']:
            val = bool(val)
        else:
            continue  # fragments
        kwargs[key] = val
    return PatientStatus(**kwargs)

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
    cat arrivalTxt | %prog [-L low_date] [-H high_date] run_descr.yaml
    """)
    parser.add_option('-L', '--low', action='store', type='int',
                      help='minimum date to include',
                      default=DEFAULT_LOW_DATE)
    parser.add_option('-H', '--high', action='store', type='int',
                      help='maximum date to include',
                      default=DEFAULT_HIGH_DATE)

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('An input yaml file matching %s is required' % INPUT_SCHEMA)

    parser.destroy()
    lowDate = opts.low
    highDate = opts.high

    #inputDict = tu.readModelInputs(args[0])
    #facDict = tu.getFacDict(inputDict)

    patientLocs = {}

    lastArriving = None
    for line in sys.stdin.readlines():
        if 'DEBUG' in line and 'arrives' in line:
            try:
                patName, dstName, date = parseArrivalLine(line)
                bits = dstName.split('_')
                try:
                    int(bits[-1])
                    bits = bits[:-2]  # Some place names end in a ward tier and number
                except:
                    pass
                newLoc = '_'.join(bits[4:7])
                if patName not in patientLocs:
                    bits = patName.split('_')
                    if len(bits) == 8:
                        startLoc = bits[6]
                    else:
                        startLoc = '_'.join(bits[:-1])
                    patientLocs[patName] = [(startLoc, 0)]
                patientLocs[patName].append((newLoc, date))
                lastArriving = patName
            except ParseLineError, e:
                print e
        elif 'DEBUG' in line and 'status is' in line:
            try:
                patientStatus = parseStatusLine(line, lastArriving)
                oldEnt = patientLocs[lastArriving].pop()
                assert len(oldEnt) == 2, 'Two status lines for %s %s?' % (lastArriving, oldEnt)
                newLoc, date = oldEnt
                patientLocs[lastArriving].append((newLoc, date, patientStatus))
            except ParseLineError, e:
                print e
        else:
            pass

    targetLoc = 'NORT_660_H'
    targetDiagClA = DiagClassA.SICK
    relEventL = []
    placeEvtD = defaultdict(list)
    for patName, eventL in patientLocs.items():
        #print patName, 'eventL: ', eventL
        oldLoc, oldDate = eventL[0]
        if len(eventL[1]) >= 3:
            oldPS = eventL[1][2]
        else:
            oldPS = None
        eventL = eventL[1:]
        inFlag = (oldLoc == targetLoc and oldPS.diagClassA == targetDiagClA)
        if inFlag:
            relEventL.append((oldDate, patName, oldLoc, oldPS))
        oldTier = DIAG_TIER_MAP[oldPS.diagClassA] if oldPS is not None else None
        oldAddr = (oldLoc, oldTier)
        placeEvtD[oldAddr].append((oldDate, patName, oldLoc, oldPS))
        for evt in eventL:
            if len(evt) >= 3:
                newLoc, newDate, newPS = evt
            else:
                newLoc, newDate = evt
                newPS = None
            if inFlag or (newLoc == targetLoc and newPS.diagClassA == targetDiagClA):
                relEventL.append((newDate, patName, newLoc, newPS))
            newTier = DIAG_TIER_MAP[newPS.diagClassA] if newPS is not None else None
            newAddr = (newLoc, newTier)
            if newAddr != oldAddr:
                # capture departure event
                placeEvtD[oldAddr].append((newDate, patName, newLoc, newPS))
                placeEvtD[newAddr].append((newDate, patName, newLoc, newPS, oldAddr))
            else:
                placeEvtD[newAddr].append((newDate, patName, newLoc, newPS))
            oldLoc, oldDate, oldPS, oldTier, oldAddr = newLoc, newDate, newPS, newTier, newAddr
            inFlag = (oldLoc == targetLoc and oldPS.diagClassA == targetDiagClA)
    relEventL.sort()
    newPlaceEvtD = {}
    for key, lst in placeEvtD.items():
        lst.sort()
        newPlaceEvtD[key] = lst
    placeEvtD = newPlaceEvtD
    for a, b in zip(enumerate(relEventL),
                    enumerate(placeEvtD[(targetLoc, DIAG_TIER_MAP[targetDiagClA])])):
        aInd, (aDate, aPat, aLoc, aPS) = a
        bInd, (bDate, bPat, bLoc, bPS) = b
        aDCA = aPS.diagClassA if aPS is not None else None
        bDCA = bPS.diagClassA if bPS is not None else None
        print '%3d %3d %15s %10s %3d %3d %15s %10s' % (aInd, aDate, aPat, aDCA, bInd, bDate, bPat, bDCA)

    with open('ofile.pkl', 'w') as f:
        pickle.dump(placeEvtD, f, 2)

if __name__ == "__main__":
    main()
