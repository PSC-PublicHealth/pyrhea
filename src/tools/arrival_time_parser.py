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
import signal
import optparse
import yaml
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.stats import norm
import logging
import logging.config
import yaml
import random

import schemautils
import pyrheautils
from pyrhea import getLoggerConfig, checkInputFileSchema, loadPathogenImplementations
from facilitybase import CareTier, PthStatus
from map_transfer_matrix import parseFacilityData

import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools

SCHEMA_DIR = os.path.join(os.path.dirname(__file__), '../schemata')
INPUT_SCHEMA = 'rhea_input_schema.yaml'

LOGGER = None

hospCatList = ['HOSPITAL', 'LTAC', 'LTACH']

inputTxt = 'arrivals_year_run_ChicagoLand.txt'

directTransferData = ['$(MODELDIR)/direct_transfer_counts.yaml']
# indirectTransferData = ['$(MODELDIR)/constants/hosp_indirect_transfer_matrix.yaml',
#                         '$(MODELDIR)/constants/nh_readmit_transfer_matrix.yaml']
indirectTransferData = ['$(MODELDIR)/hosp_indirect_transfer_counts.yaml',
                        '$(MODELDIR)/nh_readmit_transfer_counts.yaml']

def allVisible(loc, facDict, patName):
    return True

def isVisibleDirect(loc, facDict, patName):
    """Return True if the given facility name is visible for direct transfers.  Loc may be None."""
    return (loc is not None and facDict[loc]['category'] != 'COMMUNITY')


def isVisibleIndirect(loc, facDict, patName):
    """Return True if the given facility name is visible for direct transfers.  Loc may be None."""
    return isVisibleDirect(loc, facDict, patName)

INVISIBLE_HOSP_PATIENTS = set()
VISIBLE_HOSP_PATIENTS = set()

def isVisibleSomeHosp(loc, facDict, patName):
    global INVISIBLE_HOSP_PATIENTS
    global VISIBLE_HOSP_PATIENTS
    if loc is None:
        return False
    else:
        cat = facDict[loc]['category']
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


def accumulate(mtx, tplL, lowThresh, highThresh, facIdxTbl, facDict, commIdx):
    counts = 0
    if len(tplL) > 1:
        oldLoc, oldArrive, oldDepart = tplL[0]
        oldInHosp = (facDict[oldLoc]['category'] in hospCatList)
        if oldInHosp:
            lastHospDate = oldDepart
        else:
            lastHospDate = None
        for newLoc, newArrive, newDepart in tplL[1:]:
            newInHosp = (facDict[newLoc]['category'] in hospCatList)
    #             if (newInHosp and lastHospDate is not None and (newDate - lastHospDate) <= 365):
    #                 print 'hosp-to-hosp transfer %s -> %s' % (oldLoc, newLoc)
            delta = newArrive - oldDepart
            if delta <= highThresh and delta >= lowThresh:
                oldIdx = facIdxTbl[oldLoc] if oldLoc in facIdxTbl else commIdx
                newIdx = facIdxTbl[newLoc] if newLoc in facIdxTbl else commIdx
                mtx[oldIdx, newIdx] += 1
                counts += 1
            oldLoc, oldArrive, oldDepart, oldInHosp = newLoc, newArrive, newDepart, newInHosp
    return counts


def plotLogImg(mtx, minVal, maxVal):
    fltIm = np.log10(mtx.astype(np.float32) + 1.0)
    myMin = np.log10(float(minVal) + 1.0)
    myMax = np.log10(float(maxVal) + 1.0)
#     fltIm = mtx.astype(np.float32)
#     myMin, myMax = minVal, maxVal
    return plt.imshow(fltIm,
#                       cmap=cm.RdYlGn,
                      cmap=cm.viridis,
                      vmin=myMin, vmax=myMax,
#                       interpolation='bilinear',
#                       origin='lower',
#                       extent=[-3, 3, -3, 3]
               )

def plotImg(mtx, minVal, maxVal):
    fltIm = mtx.astype(np.float32)
    myMin, myMax = minVal, maxVal
    return plt.imshow(fltIm,
                      cmap=cm.bwr,
                      vmin=myMin, vmax=myMax,
                      )

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


def main():
    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)

    global LOGGER
    #logging.config.dictConfig(getLoggerConfig())
    LOGGER = logging.getLogger(__name__)

    parser = optparse.OptionParser(usage="""
    %prog run_descr.yaml
    """)

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('An input yaml file matching %s is required' % INPUT_SCHEMA)

    parser.destroy()

    inputDict = checkInputFileSchema(args[0],
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    schemautils.setSchemaBasePath(SCHEMA_DIR)
    if 'modelDir' in inputDict:
        pyrheautils.PATH_STRING_MAP['MODELDIR'] = pyrheautils.pathTranslate(inputDict['modelDir'])
    if 'pathTranslations' in inputDict:
        for elt in inputDict['pathTranslations']:
            pyrheautils.PATH_STRING_MAP[elt['key']] = elt['value']

#    loadPthInfo(inputDict['pathogenImplementationDir'])
    facDirs = [pyrheautils.pathTranslate(dct) for dct in inputDict['facilityDirs']]
    facDict = parseFacilityData(facDirs)
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
    directMtx = np.zeros([commIdx+1, commIdx+1], dtype=np.int32)
    indirectMtx = np.zeros_like(directMtx)

    patientLocs = {}

    with open(os.path.join(os.path.dirname(__file__), inputTxt), 'rU') as f:
        for line in f.readlines():
            if 'birth' in line: print line
            words =  line.split()
            if len(words) != 6 or words[3] != 'arrives':
                print 'Bad line format: %s' % line
                continue
            patName = words[2]
            dstName = words[4]
            date = int(words[5])
            if patName in patientLocs:
                bits = dstName.split('_')
                newLoc = '_'.join(bits[4:7])
                patientLocs[patName].append((newLoc, date))
            else:
                bits = patName.split('_')
                if len(bits) == 8:
                    startLoc = bits[6]
                else:
                    startLoc = '_'.join(bits[:-1])
                patientLocs[patName] = [(startLoc, 0)]

    directCounts = 0
    indirectCounts = 0
    for patName, tupleL in patientLocs.items():

        #filtTupleL = filterPath(addDepartureTime(tupleL), isVisibleDirect, facDict, patName)
        filtTupleL = filterPath(addDepartureTime(tupleL), isVisibleSomeHosp, facDict, patName)
        filtTupleL = mergeSelfTransfers(filtTupleL)
        directCounts += accumulate(directMtx, filtTupleL, 0, 3, facIdxTbl, facDict, commIdx)

        #filtTupleL = filterPath(addDepartureTime(tupleL), isVisibleIndirect, facDict, patName)
        filtTupleL = filterPath(addDepartureTime(tupleL), isVisibleSomeHosp, facDict, patName)
        filtTupleL = mergeSelfTransfers(filtTupleL)
        indirectCounts += accumulate(indirectMtx, filtTupleL, 4, 365, facIdxTbl, facDict, commIdx)

    print 'directCounts = %s' % directCounts
    print 'direct transfer matrix follows'
    print directMtx
    print 'indirectCounts = %s' % indirectCounts
    print 'indirect transfer matrix follows'
    print indirectMtx
    maxVal = max(np.max(directMtx), np.max(indirectMtx))

    ax11 = plt.subplot(2, 2, 1)
    ax11.set_title('log direct')
    plotLogImg(directMtx, 0.0, float(maxVal))

    ax12 = plt.subplot(2, 2, 2)
    ax12.set_title('log indirect')
    plotLogImg(indirectMtx, 0.0, float(maxVal))

    plt.colorbar(ax=[ax11, ax12])

    measDirectMtx = mtxFromYaml(directTransferData, directMtx, facIdxTbl, commIdx)
    sclMtx = (float(np.sum(directMtx))/float(np.sum(measDirectMtx))) * measDirectMtx
    #deltaMtx = (2.0*(directMtx - sclMtx)/(directMtx + sclMtx))
    directDeltaMtx = directMtx - sclMtx
    #directDeltaMtx[0:100, 0:100] = 0.0

    measIndirectMtx = mtxFromYaml(indirectTransferData, indirectMtx, facIdxTbl, commIdx)
    sclMtx = (float(np.sum(indirectMtx))/float(np.sum(measIndirectMtx))) * measIndirectMtx
    #deltaMtx = (2.0*(indirectMtx - sclMtx)/(directMtx + sclMtx))
    indirectDeltaMtx = indirectMtx - sclMtx
    #indirectDeltaMtx[0:100, 0:100] = 0.0
    lim = max(np.max(np.fabs(directDeltaMtx)), np.max(np.fabs(indirectDeltaMtx)))

    ax21 = plt.subplot(2, 2, 3)
    ax21.set_title('normalized direct delta')
    plotImg(directDeltaMtx, -lim, lim)

    ax22 = plt.subplot(2, 2, 4)
    ax22.set_title('normalized indirect delta')
    plotImg(indirectDeltaMtx, -lim, lim)

    plt.colorbar(ax=[ax21, ax22])
    plt.show()

if __name__ == "__main__":
    main()
