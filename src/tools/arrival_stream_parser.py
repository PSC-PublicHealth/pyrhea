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


def parseDeathLine(line):
    words = line.split()
    while words[0].strip() == '(Pdb)':
        words = words[1:]
    assert words[4:6] == ['died', 'at'], 'Bad line format: %s' % line
    if len(words) == 8:
        patName = words[3]
        locName = None
        date = int(words[-1])
    elif len(words) == 10:
        patName = words[3]
        locName = words[6]
        date = int(words[-1])
    else:
        raise RuntimeError('death line format error on %s' % line)
#     words = patName.split('_')
#     while words:
#         lW = words.pop()
#         try:
#             idx = int(lW)
#             break
#         except:
#             pass
#     patName = '_'.join([birthLocFromPatientName(patName), str(idx)])
    return patName, locName, date


def tierFromPatientStatus(patientStatus):
    if patientStatus is not None:
        tier = DIAG_TIER_MAP[patientStatus.diagClassA]
        if patientStatus.overall == PatientOverallHealth.FRAIL and tier == CareTier.HOME:
            tier = CareTier.NURSING
    else:
        tier = None
    return tier


def birthLocFromPatientName(patName):
    words = patName.split('_')
    if 'Patch' in words:
        offset = words.index('Patch')
        if offset < 0:
            raise RuntimeError('name format parsing error 1 on %s' % patName)
        words = words[offset+3:]
    if 'birth' in words:
        words.remove('birth')
    try:
        int(words[-1])
    except Exception, e:
        raise RuntimeError('name format parsing error 2 on %s' % patName)
    words = words[:-1]
    try:
        int(words[-1])
        words = words[:-2]  # Some place names end in a ward tier and number
    except:
        pass
    rslt = '_'.join(words)
    #print '%s -> %s' % (patName, rslt)
    return rslt


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
                bits = patName.split('_')
                if len(bits) == 8:
                    startLoc = bits[6]
                else:
                    startLoc = '_'.join(bits[:-1])
            except ParseLineError, e:
                print e
        elif 'DEBUG' in line and 'status is' in line:
            # Always immediately follows the associated 'arrives' line
            try:
                patientStatus = parseStatusLine(line, patName)
                if patName not in patientLocs:
                    # add a creation event at time 0
                    patientLocs[patName] = [(startLoc, 0, None, True)]
                patientLocs[patName].append((newLoc, date, patientStatus, True))
            except ParseLineError, e:
                print e
        elif 'DEBUG' in line and 'died at' in line:
            try:
                patName, locName, date = parseDeathLine(line)
                if patName not in patientLocs:
                    birthLoc = birthLocFromPatientName(patName) if locName is None else locName
                    patientLocs[patName] = [(birthLoc, 0, None, True)]
                if locName is None:
                    oldLoc, oldDate, oldPS, oldAlive = patientLocs[patName][-1]
                    locName = oldLoc
                patientLocs[patName].append((locName, date, None, False))
            except ParseLineError, e:
                print e
        else:
            pass

    sampL = random.sample(patientLocs, 3)
    for patNm in sampL:
        print '%s: %s' % (patNm, patientLocs[patNm])

    placeEvtD = defaultdict(list)
    for patName, eventL in patientLocs.items():
        #print patName, 'eventL: ', eventL
        oldLoc, oldDate, oldPS, oldAlive = eventL[0]
        oldPthStatus = -1 if oldPS is None else oldPS.pthStatus
        eventL = eventL[1:]
        oldTier = tierFromPatientStatus(oldPS)
        oldAddr = (oldLoc, oldTier)
        pLPSC = None  # 'place of last pth status change'
        placeEvtD[oldAddr].append((oldDate, patName, oldLoc, oldPS, pLPSC, oldAlive))
        for evt in eventL:
            if len(evt) >= 4:
                newLoc, newDate, newPS, newAlive = evt
            else:
                raise RuntimeError('Bad event format')
            if newLoc is None and oldLoc is not None:
                newLoc = oldLoc
            newTier = tierFromPatientStatus(newPS)
            newAddr = (newLoc, newTier)
            newPthStatus = -1 if newPS is None else newPS.pthStatus
            if newPthStatus != oldPthStatus:
                pLPSC = newLoc
            if newAddr != oldAddr:
                # capture departure event
                placeEvtD[oldAddr].append((newDate, patName, newLoc, newPS, pLPSC, newAlive))
                if newAlive:
                    placeEvtD[newAddr].append((newDate, patName, newLoc, newPS,
                                               oldAddr, pLPSC, newAlive))
            else:
                placeEvtD[newAddr].append((newDate, patName, newLoc, newPS, pLPSC, newAlive))
            oldLoc, oldDate, oldPS, oldTier = newLoc, newDate, newPS, newTier
            oldAddr, oldPthStatus, oldAlive = newAddr, newPthStatus, newAlive
    newPlaceEvtD = {}
    for key, lst in placeEvtD.items():
        lst.sort()
        newPlaceEvtD[key] = lst
    placeEvtD = newPlaceEvtD

    with open('ofile.pkl', 'w') as f:
        pickle.dump(placeEvtD, f, 2)

if __name__ == "__main__":
    main()
