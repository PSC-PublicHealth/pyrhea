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

_rhea_svn_id_ = "$Id: map_facilities.py 133 2016-03-07 20:12:44Z welling $"

import os.path
import optparse
import logging
import signal
import yaml
import re
from datetime import datetime
import numpy as np
from map_transfer_matrix import parseFacilityData
from notes_plotter import findFacImplCategory
from notes_plotter import importNotes
from pyrheautils import loadModulesFromDir
import schemautils
import geodatamanager
import geomap

SCHEMA_DIR = '../schemata'
INPUT_SCHEMA = 'rhea_input_schema.yaml'

logger = None


def buildLclWorkVec(orderedFacList, facDict, implementationDir, catToImplDict):
    implDict = {}
    for mod in loadModulesFromDir(implementationDir, ['category', 'estimateWork']):
        implDict[mod.category] = mod
    wtDict = {}
    for abbrev, descr in facDict.items():
        wtDict[abbrev] = implDict[catToImplDict[descr['category']]].estimateWork(descr)
    for key in wtDict.keys():
        if wtDict[key] == '':
            print 'No weight for %s' % key
            wtDict[key] = 0.0
    lclWorkVec = np.array([float(wtDict[k]) for k in orderedFacList])
    return lclWorkVec


def buildTransferMatrix(orderedFacList, facDict, notesPath):
    idx = {nm: i for i, nm in enumerate(orderedFacList)}
    notes = importNotes(notesPath)
    transferMatrix = np.zeros((len(orderedFacList), len(orderedFacList)))
    for fmAbbrev in orderedFacList:
        fmKey = '%s_%s' % (facDict[fmAbbrev]['category'].upper(), fmAbbrev)
        if fmKey in notes:
            for k, v in notes[fmKey].items():
                if k.endswith('_transfer'):
                    toAbbrev = k[:-len('_transfer')]
                    transferMatrix[idx[fmAbbrev], idx[toAbbrev]] += v
        else:
            print 'no notes for %s' % fmKey
    return transferMatrix


def checkInputFileSchema(fname, schemaFname):
    if logger is None:
        myLogger = logging.getLogger(__name__)
    else:
        myLogger = logger
    try:
        with open(fname, 'rU') as f:
            inputJSON = yaml.safe_load(f)
        validator = schemautils.getValidator(os.path.join(os.path.dirname(__file__), schemaFname))
        nErrors = sum([1 for e in validator.iter_errors(inputJSON)])  # @UnusedVariable
        if nErrors:
            myLogger.error('Input file violates schema:')
            for e in validator.iter_errors(inputJSON):
                myLogger.error('Schema violation: %s: %s' %
                               (' '.join([str(word) for word in e.path]), e.message))
            raise RuntimeError('Input file violates its schema')
        else:
            return inputJSON
    except Exception, e:
        myLogger.error('Error checking input against its schema: %s' % e)
        raise RuntimeError('Input file violates its schema')


class WorkPartition(object):
    """
    This class embodies the work associated with a single partition of the parallel job.
    """
    def __init__(self, orderedFacList, workVec, transferMatrix):
        self.idxFac = {idx: nm for idx, nm in enumerate(orderedFacList)}
        self.facIdx = {nm: idx for idx, nm in enumerate(orderedFacList)}
        self.wVec = workVec
        self.tMatrix = transferMatrix

    def swap(self, idx1, idx2):
        """
        Completely swap two indices.
        """
        self.wVec[[idx1, idx2]] = self.wVec[[idx2, idx1]]
        self.tMatrix[[idx1, idx2], :] = self.tMatrix[[idx2, idx1], :]
        self.tMatrix[:, [idx1, idx2]] = self.tMatrix[:, [idx2, idx1]]
        fac1 = self.idxFac[idx1]
        fac2 = self.idxFac[idx2]
        self.facIdx[fac1] = idx2
        self.facIdx[fac2] = idx1
        self.idxFac[idx1] = fac2
        self.idxFac[idx2] = fac1

    def dump(self):
        """
        Print a representation for diagnostic purposes
        """
        for i, abbrev in self.idxFac.items():
            wt = self.wVec[i]
            iBack = self.facIdx[abbrev]
            print '%d: %s: %f %d' % (i, abbrev, wt, iBack)
        print self.tMatrix

    def calcInternalWork(self):
        """
        Calculate the current total work for this partition
        """ 
        if len(self.wVec):
            return np.sum(self.wVec)
        else:
            return 0.0

    def split(self, idx):
        """
        Split into two WorkPartitions, with the first idx indices in the first partition
        and the remaining indices in the second.
        """
        sz = self.wVec.shape[0]
        wp1 = WorkPartition([self.idxFac[i] for i in xrange(idx)],
                            self.wVec[:idx],
                            self.tMatrix[:idx, :idx])
        wp2 = WorkPartition([self.idxFac[i] for i in xrange(idx, sz)],
                            self.wVec[idx:],
                            self.tMatrix[idx:, idx:])
        return wp1, wp2

    def calcTransferWork(self, splitIdx):
        """
        If this WorkPartition is split at splitIdx, calculate the work associated with transfers
        across the split.
        """
        wp1, wp2 = self.split(splitIdx)
        return self.tMatrix.sum() - (wp1.tMatrix.sum() + wp2.tMatrix.sum())

    def sortBy(self, func):
        """
        func(facAbbrev) is expected to return a scalar value.  The elements in this WorkPartition
        are re-ordered according to that scalar in ascending order
        """
        sortVec = []
        for abbrev in self.idxFac.values():
            sortVec.append((func(abbrev), abbrev))
        sortVec.sort()
        reorderedAbbrevList = [abbrev for val, abbrev in sortVec]  # @UnusedVariable
        for newIdx, abbrev in enumerate(reorderedAbbrevList):
            oldIdx = self.facIdx[abbrev]
            self.swap(newIdx, oldIdx)

    def facIter(self):
        """
        Returns an iterator over the ordered facility abbreviations
        """
        for idx in xrange(len(self.idxFac)):
            yield self.idxFac[idx]

    def __len__(self):
        return len(self.idxFac)


class PartitionSet(object):
    def __init__(self, workPartition):
        self.fullWP = workPartition

    def partition(self, nParts):
        pass

    def wPIter(self):
        yield self.fullWP


class BinaryPartitionSet(PartitionSet):
    def __init__(self, fullWP, lossFunc, sortFuncList, depth=0):
        super(BinaryPartitionSet, self).__init__(fullWP)
        self.lossFunc = lossFunc
        self.sortFuncList = sortFuncList
        self.depth = depth
        self.kid1 = None
        self.kid2 = None

    def partition(self, nParts):
        minMinLoss = None
        for sortFunc in self.sortFuncList:
            self.fullWP.sortBy(sortFunc)

            bestSplit = None
            minLoss = None
    #         f = open('/tmp/junk.curve', 'a')
            for splitHere in xrange(len(self.fullWP)):
                botWP, topWP = self.fullWP.split(splitHere)
                lossVal = self.lossFunc(botWP.calcInternalWork(),
                                        topWP.calcInternalWork(),
                                        self.fullWP.calcTransferWork(splitHere))
    #             f.write('%s %s\n' % (splitHere, lossVal))
                if minLoss is None or lossVal < minLoss:
                    bestSplit = splitHere
                    minLoss = lossVal
    #         f.write('\n')
    #         f.write('\n')
    #         f.close()
            if minMinLoss is None or minLoss < minMinLoss:
                bestFunc = sortFunc
                bestBestSplit = bestSplit
                minMinLoss = minLoss
        self.fullWP.sortBy(bestFunc)
        botWP, topWP = self.fullWP.split(bestBestSplit)
        print ('bestSplit %s of %s in direction %s -> minLoss: %s'
               % (bestBestSplit, len(self.fullWP),
                  self.sortFuncList.index(bestFunc),
                  minMinLoss))

        subParts = nParts / 2
        if subParts > 1:
            self.kid1 = BinaryPartitionSet(botWP, self.lossFunc, self.sortFuncList,
                                           depth=self.depth+1)
        else:
            self.kid1 = PartitionSet(botWP)
        self.kid1.partition(subParts)

        subParts = nParts - subParts
        if subParts > 1:
            self.kid2 = BinaryPartitionSet(topWP, self.lossFunc, self.sortFuncList,
                                           depth=self.depth+1)
        else:
            self.kid2 = PartitionSet(topWP)
        self.kid2.partition(subParts)

    def wPIter(self):
        if self.kid1 is None:
            raise RuntimeError('Call partition first!')
        for wP in self.kid1.wPIter():
            yield wP
        for wP in self.kid2.wPIter():
            yield wP


def lossFunc(work1, work2, crossWork):
    workDif = work1 - work2
    wt1 = 1.0
    #wt2 = 0.005
    wt2 = 0.002
    return wt1*(workDif*workDif) + wt2*(crossWork*crossWork)


def drawMap(partitionSet, abbrevTractDict, facDict, geoDataPathList,
            stateCodeRE, countyCodeRE, countySet, catToImplDict):
    ctrLon = sum([r['longitude'] for r in facDict.values()]) / len(facDict)
    ctrLat = sum([r['latitude'] for r in facDict.values()]) / len(facDict)

    myMap = geomap.Map(geoDataPathList, stateCodeRE, countyCodeRE, ctrLon, ctrLat,
                       annotate=False,
                       regexCodes=True)
    LTRED = '#cc6666'
    RED = '#cc0000'
    LTMAGENTA = '#cc66cc'
    MAGENTA = '#cc00cc'
    LTBLUE = '#6666cc'
    BLUE = '#0000cc'
    LTCYAN = '#66cccc'
    CYAN = '#00cccc'
    LTGREEN = '#66cc66'
    GREEN = '#00cc00'
    LTYELLOW = '#cccc66'
    YELLOW = '#cccc00'

    tractGeoIDDict = {}
    for geoID in myMap.tractIDList():
        propRec = myMap.getTractProperties(geoID)
        tractGeoIDDict[(propRec['STATE'] + propRec['COUNTY'], propRec['TRACT'])] = geoID

    clrTupleSeq = [(LTRED, RED), (LTMAGENTA, MAGENTA), (LTBLUE, BLUE),
                   (LTCYAN, CYAN), (LTGREEN, GREEN), (LTYELLOW, YELLOW)]
    for idx, wP in enumerate(partitionSet.wPIter()):
        clr1, clr2 = clrTupleSeq[idx % len(clrTupleSeq)]
        for abbrev in wP.facIter():
            if abbrev in abbrevTractDict:
                geoID = tractGeoIDDict[abbrevTractDict[abbrev]]
                myMap.plotTract(geoID, clr1)
            else:
                rec = facDict[abbrev]
                implStr = catToImplDict[rec['category']].upper()
                mrk = {'HOSPITAL': '*', 'LTAC': '+', 'NURSINGHOME': 'o'}[implStr]
                myMap.plotMarker(rec['longitude'], rec['latitude'], mrk, rec['abbrev'], clr2)

    myMap.draw()


def isPowerOfTwo(num):
    """
    Thanks to http://code.activestate.com/recipes/577514-chek-if-a-number-is-a-power-of-two/ which
    has a much nicer solution than I was contemplating.
    """
    return ((num & (num - 1)) == 0) and num != 0


def main():
    global logger

    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)

    logging.basicConfig()

    defaultOutputName = "partition.yaml"
    defaultNumParts = 16

    parser = optparse.OptionParser(usage="""
    %prog [-v][-d][-n notes1.pkl [-n notes2.pkl [...]] [-o out.yaml] input.yaml
    """)
    parser.add_option("-v", "--verbose", action="store_true",
                      help="verbose output")
    parser.add_option("-d", "--debug", action="store_true",
                      help="debugging output")
    parser.add_option("-n", "--notes", action="append",
                      help="pickled notes file to provide transfer data")
    parser.add_option("-o", "--out", action="store", default=defaultOutputName,
                      help="specifies output file name (default %s)" % defaultOutputName)
    parser.add_option("-p", "--parts", action="store", type="int", default=defaultNumParts,
                      help=(("Number of parts for the partition (default %d). Should be a"
                            " power of 2.") % defaultNumParts))

    opts, args = parser.parse_args()

    if opts.debug:
        logging.getLogger('root').setLevel(logging.DEBUG)
    elif opts.verbose:
        logging.getLogger('root').setLevel(logging.INFO)
    else:
        logging.getLogger('root').setLevel(logging.WARNING)
    logger = logging.getLogger(__name__)

    if opts.notes:
        notesFileList = opts.notes
    else:
        notesFileList = []

    if len(args) == 1:
        inputPath = args[0]
    else:
        parser.error("A YAML-format file specifying prototype model parameters must be specified.")

    if isPowerOfTwo(opts.parts):
        numParts = opts.parts
    else:
        parser.error("The number of parts should be a power of two.")

    outputPath = opts.out
    parser.destroy()

    inputDict = checkInputFileSchema(inputPath, os.path.join(SCHEMA_DIR, INPUT_SCHEMA))

    gDM = geodatamanager.GeoDataManager()
#     stateCodeRE = '06'
#     countyCodeRE = '059'
    countyCodeRE = r'.*'
    stateCodeRE = r'.*'

    facDirList = inputDict['facilityDirs']
    implementationDir = inputDict['facilityImplementationDir']

    facDict = parseFacilityData(facDirList)

    catNames = list(set([rec['category'] for rec in facDict.values()]))

    if 'facilitySelectors' in inputDict:
        facImplRules = [(re.compile(rule['category']), rule['implementation'])
                       for rule in inputDict['facilitySelectors']]
    else:
        facImplRules = [(re.compile(cat), cat)
                        for cat in catNames]  # an identity map

    catToImplDict = {cat: findFacImplCategory(implementationDir, facImplRules, cat)
                     for cat in catNames}

    schemautils.setSchemaBasePath(SCHEMA_DIR)

    orderedFacList = facDict.keys()[:]
    orderedFacList.sort()

    transferMatrix = np.zeros((len(orderedFacList), len(orderedFacList)))
    if notesFileList:
        for notesPath in notesFileList:
            transferMatrix += buildTransferMatrix(orderedFacList, facDict, notesPath)
    #         diagPairs = zip(orderedFacList, transferMatrix.diagonal())
    #         for abbrev, val in diagPairs:
    #             assert val == 0 or facDict[abbrev]['category'].upper() in ['HOSPITAL'], 'fail %s %s' % (abbrev, val)
        transferMatrix /= float(len(notesFileList))

    lclWorkVec = buildLclWorkVec(orderedFacList, facDict, implementationDir,
                                  catToImplDict)

    fullWP = WorkPartition(orderedFacList, lclWorkVec, transferMatrix)

    def latSortFun(abbrev):
        return facDict[abbrev]['latitude']

    def lonSortFun(abbrev):
        return facDict[abbrev]['longitude']

    partitionSet = BinaryPartitionSet(fullWP, lossFunc, [latSortFun, lonSortFun])
    partitionSet.partition(numParts)

    abbrevTractDict = {}
    countySet = set()
    for abbrev, rec in facDict.items():
        if 'FIPS' in rec:
            countySet.add(abbrev)
            gDM.addFIPS(rec['FIPS'])
        if 'censusTract' in rec:
            assert 'FIPS' in rec, '%s has a tract but no FIPS code?' % abbrev
            abbrevTractDict[abbrev] = (rec['FIPS'], rec['censusTract'])

    drawMap(partitionSet, abbrevTractDict, facDict,
            gDM.getGeoDataFileList(), stateCodeRE, countyCodeRE, countySet, catToImplDict)

    partitionByAbbrev = {}
    for idx, wP in enumerate(partitionSet.wPIter()):
        for abbrev in wP.facIter():
            partitionByAbbrev[abbrev] = idx
    allData = {'generatedBy': __file__,
               'generatedTimeUTC': datetime.utcnow().isoformat(),
               'notesInputs': notesFileList,
               'modelInput': inputPath,
               'partition': partitionByAbbrev}
    with open(outputPath, 'w') as f:
        yaml.dump(allData, f)

if __name__ == "__main__":
    main()
