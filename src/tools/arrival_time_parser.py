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

import schemautils
import pyrheautils
import tools_util as tu
from pyrhea import getLoggerConfig, checkInputFileSchema, loadPathogenImplementations
from facilitybase import CareTier, PthStatus
from map_transfer_matrix import parseFacilityData
from plotting_utilities import pltMtxImageFig

import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools

SCHEMA_DIR = os.path.join(os.path.dirname(__file__), '../schemata')
INPUT_SCHEMA = 'rhea_input_schema.yaml'

LOGGER = None

HOSP_CAT_LIST = ['HOSPITAL', 'LTAC', 'LTACH']

DEFAULT_INPUT_TXT = 'arrivals_year_run_ChicagoLand.txt'
DEFAULT_LOW_DATE = 0
DEFAULT_HIGH_DATE = -1  # meaning include all records

DIRECT_TRANSFER_DATA = ['$(MODELDIR)/direct_transfer_counts.yaml']
INDIRECT_TRANSFER_DATA = ['$(MODELDIR)/hosp_indirect_transfer_counts.yaml',
                           '$(MODELDIR)/nh_readmit_transfer_counts.yaml']
#INDIRECT_TRANSFER_DATA = ['$(MODELDIR)/hosp_indirect_transfer_counts.yaml',
#                          '$(MODELDIR)/nh_readmit_fake_transfer_counts.yaml']

def allVisible(loc, facDict, patName):
    return True

def isVisibleDirect(loc, facDict, patName):
    """Return True if the given facility name is visible for direct transfers.  Loc may be None."""
    return (loc is not None and facDict[loc]['category'] != 'COMMUNITY')

def isVisibleIndirect(loc, facDict, patName):
    """Return True if the given facility name is visible for indirect transfers.  Loc may be None."""
    return isVisibleDirect(loc, facDict, patName)

INVISIBLE_HOSP_PATIENTS = set()
VISIBLE_HOSP_PATIENTS = set()

def isVisibleSomeHosp(loc, facDict, patName):
    global INVISIBLE_HOSP_PATIENTS
    global VISIBLE_HOSP_PATIENTS
    if loc is None:
        return False
    else:
        if loc in facDict:
            cat = facDict[loc]['category']
        elif loc.endswith('cache'):
            cat = facDict[loc[:-5]]['category']

        if cat == 'COMMUNITY':
            return False
        elif cat == 'HOSPITAL':
            if patName in INVISIBLE_HOSP_PATIENTS:
                return False
            elif patName in VISIBLE_HOSP_PATIENTS:
                return True
            else:
                rnd = random.random()
                if rnd <= 0.0:
                    INVISIBLE_HOSP_PATIENTS.add(patName)
                    return False
                else:
                    VISIBLE_HOSP_PATIENTS.add(patName)
                    return True
        else:
            return True


def addDepartureTime(tupleL):
    """
    Arriving tuples have the form (placeName, arrivalDate).  Return (placeName, arrivalDate, departureDate).
    """
    rslt = []
    prevLoc, prevArrival = tupleL[0]
    for newLoc, newArrival in tupleL[1:]:
        rslt.append((prevLoc, prevArrival, newArrival))
        prevLoc, prevArrival = newLoc, newArrival
    rslt.append((prevLoc, prevArrival, prevArrival+1))
    return rslt


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


def accumulate(mtx, tplL, lowThresh, highThresh, lowDate, highDate, facIdxTbl, facDict, commIdx):
    counts = 0
    if len(tplL) > 1:
        oldLoc, oldArrive, oldDepart = tplL[0]
        oldInHosp = (facDict[oldLoc]['category'] in HOSP_CAT_LIST)
        if oldInHosp:
            lastHospDate = oldDepart
        else:
            lastHospDate = None
        for newLoc, newArrive, newDepart in tplL[1:]:
            newInHosp = (facDict[newLoc]['category'] in HOSP_CAT_LIST)
            if newArrive >= lowDate and newArrive <= highDate:
                delta = newArrive - oldDepart
                if delta <= highThresh and delta >= lowThresh:
                    oldIdx = facIdxTbl[oldLoc] if oldLoc in facIdxTbl else commIdx
                    newIdx = facIdxTbl[newLoc] if newLoc in facIdxTbl else commIdx
                    mtx[oldIdx, newIdx] += 1
                    counts += 1
            oldLoc, oldArrive, oldDepart, oldInHosp = newLoc, newArrive, newDepart, newInHosp
    return counts

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

def mtxFromYaml(filePathL, protoMtx, facIdxTbl, commIdx):
    mtx = np.zeros_like(protoMtx)
    for dataPth in [pyrheautils.pathTranslate(pth) for pth in filePathL]:
        with open(dataPth, 'rU') as f:
            data = yaml.safe_load(f)
            for srcLoc, dstD in data.items():
                for dstLoc, ct in dstD.items():
                    srcIdx = facIdxTbl[srcLoc] if srcLoc in facIdxTbl else commIdx
                    dstIdx = facIdxTbl[dstLoc] if dstLoc in facIdxTbl else commIdx
                    mtx[srcIdx, dstIdx] += ct
    return mtx

class ParseLineError(RuntimeError):
    pass

def parseLine(line):
    if 'birth' in line: print line
    words =  line.split()
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
    %prog --input recs.txt [-L low_date] [-H high_date] run_descr.yaml
    """)
    
    parser.add_option('-i', '--input', action='store', type='string',
                      help=('The file of text records to be parsed (default %s)'
                            % DEFAULT_INPUT_TXT),
                      default=DEFAULT_INPUT_TXT)
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
    inputTxt = opts.input
    lowDate = opts.low
    highDate = opts.high

    inputDict = tu.readModelInputs(args[0])
    facDict = tu.getFacDict(inputDict)
    print 'IMPLEMENTING SPECIAL PATCH FOR WAUK_2615_H'
    facDict['WAUK_2615_H'] = {'category':'HOSPITAL'}
    print 'FIND OUT THE REAL ANSWER AND DELETE THIS!'

    # Sort by category, then alphabetically by name
    facL = [(facDict[fac]['category'], fac) for fac in facDict.keys()
            if facDict[fac]['category'] != 'COMMUNITY']
    facL.sort()
    facL = [fac for facCat, fac in facL]
    facIdxTbl = {key:idx for idx, key in enumerate(facL)}
    commIdx = len(facL)

    patientLocs = {}

    with open(os.path.join(os.path.dirname(__file__), inputTxt), 'rU') as f:
        for line in f.readlines():
            try:
                patName, dstName, date = parseLine(line)
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

            except ParseLineError, e:
                print e

    directCounts = 0
    indirectCounts = 0
    directMtx = np.zeros([commIdx+1, commIdx+1], dtype=np.int32)
    indirectMtx = np.zeros_like(directMtx)
    accumCountsD = defaultdict(lambda:0)
    accumStayD = defaultdict(lambda: 0.0)
    netEffStart = None
    netEffEnd = None
    for patName, tupleL in patientLocs.items():

        filtTupleL = filterPath(addDepartureTime(tupleL), isVisibleSomeHosp, facDict, patName)
        filtTupleL = mergeSelfTransfers(filtTupleL)
        countD, effStart, effEnd = countVisits(filtTupleL, lowDate, highDate, facDict,
                                               accumCountsD, accumStayD)
        if netEffStart is None or effStart < netEffStart:
            netEffStart = effStart
        if netEffEnd is None or effEnd > netEffEnd:
            netEffEnd = effEnd
        directCounts += accumulate(directMtx, filtTupleL, 0, 3, lowDate, highDate,
                                   facIdxTbl, facDict, commIdx)

        #filtTupleL = filterPath(addDepartureTime(tupleL), isVisibleIndirect, facDict, patName)
        filtTupleL = filterPath(addDepartureTime(tupleL), isVisibleSomeHosp, facDict, patName)
        filtTupleL = mergeSelfTransfers(filtTupleL)
        indirectCounts += accumulate(indirectMtx, filtTupleL, 4, 365, lowDate, highDate,
                                     facIdxTbl, facDict, commIdx)

    print 'distinct patients: %s' % len(patientLocs)
    print 'visits in date range %s to %s:' % (netEffStart, netEffEnd)
    for k, v in accumCountsD.items():
        totStay = accumStayD[k]
        ratio = (float(v) * 365.)/(float(len(patientLocs)) * (netEffEnd + 1 - netEffStart))
        print '   %8s: %d = %f per patient year, ave stay %s' % (k, v, ratio, totStay/v)
    print 'directCounts = %s' % directCounts
    print 'direct transfer matrix follows'
    print directMtx
    print 'indirectCounts = %s' % indirectCounts
    print 'indirect transfer matrix follows'
    print indirectMtx


    measDirectMtx = mtxFromYaml(DIRECT_TRANSFER_DATA, directMtx, facIdxTbl, commIdx)

    measIndirectMtx = mtxFromYaml(INDIRECT_TRANSFER_DATA, indirectMtx, facIdxTbl, commIdx)

    pltMtxImageFig(directMtx, measDirectMtx, indirectMtx, measIndirectMtx)

    plt.savefig('arrival_time_plots.svg', bbox_inches='tight')
   
    np.savez('arrival_time_arrays.npz',
            indirect_simulated=indirectMtx,
            direct_simulated=directMtx,
            indirect_measured=measIndirectMtx,
            direct_measured=measDirectMtx)


if __name__ == "__main__":
    main()
