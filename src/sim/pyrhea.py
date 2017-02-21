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

import sys
import os
import yaml
from imp import load_source
from random import seed
from collections import defaultdict
import optparse
import logging
import logging.config
import pika
import cPickle as pickle
import re
import signal

import quilt.patches as patches
import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.notes.noteholder as noteholder
import schemautils
import pyrheautils

SCHEMA_DIR = os.path.join(os.path.dirname(__file__), '../schemata')
INPUT_SCHEMA = 'rhea_input_schema.yaml'

DEFAULT_OUTPUT_NOTES_NAME = 'notes.pkl'

_TRACKED_FACILITIES = []
_TRACKED_FACILITIES_SET = None

logger = None

from pyrheabase import ScenarioPolicy

def buildFacOccupancyDict(patch, timeNow):
    facTypeDict = {'day': timeNow}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    for fac in patch.allFacilities:
        tpName = type(fac).__name__
        if tpName not in facTypeDict:
            facTypeDict[tpName] = 0
            facTypeDict[tpName + '_all'] = 0

        if 0:
            if hasattr(fac, 'patientDataDict'):
                patientCount = 0
                allCount = 0
                for rec in fac.patientDataDict.values():
                    rec = pickle.loads(rec)
                    allCount += 1 + rec.prevVisits
                    if rec.departureDate is None:
                        patientCount += 1
                facTypeDict[tpName] += patientCount
                facTypeDict[tpName + '_all'] += allCount
        else:
            if hasattr(fac, 'patientStats'):
                patientCount = fac.patientStats.currentOccupancy
                allCount = fac.patientStats.totalOccupancy
                facTypeDict[tpName] += patientCount
                facTypeDict[tpName + '_all'] += allCount

    return facTypeDict


def buildFacPthDict(patch, timeNow):
    facPthDict = {'day': timeNow}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    for fac in patch.allFacilities:
        facPPC = {}
        for ward in fac.getWards():
            pPC = ward.iA.getPatientPthCounts(timeNow)
            for k, v in pPC.items():
                if k in facPPC:
                    facPPC[k] += v
                else:
                    facPPC[k] = v
        tpName = type(fac).__name__
        for k, v in facPPC.items():
            key = '%s_%s' % (tpName, k)
            if key in facPthDict:
                facPthDict[key] += v
            else:
                facPthDict[key] = v
    return facPthDict


def buildLocalOccupancyDict(patch, timeNow):
    global _TRACKED_FACILITIES_SET
    if not _TRACKED_FACILITIES_SET:
        _TRACKED_FACILITIES_SET = frozenset(_TRACKED_FACILITIES)
    facDict = {'day': timeNow}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    for fac in patch.allFacilities:
        if 0:
            if fac.abbrev in _TRACKED_FACILITIES_SET and hasattr(fac, 'patientDataDict'):
                patientCount = 0
                for rec in fac.patientDataDict.values():
                    rec = pickle.loads(rec)
                    if rec.departureDate is None:
                        patientCount += 1
                facDict[fac.abbrev] = patientCount
        else:
            if fac.abbrev in _TRACKED_FACILITIES_SET and hasattr(fac, 'patientStats'):
                patientCount = fac.patientStats.currentOccupancy
                facDict[fac.abbrev] = patientCount
            
    return facDict


def buildLocalPthDict(patch, timeNow):
    global _TRACKED_FACILITIES_SET
    if not _TRACKED_FACILITIES_SET:
        _TRACKED_FACILITIES_SET = frozenset(_TRACKED_FACILITIES)
    facPthDict = {'day': timeNow}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    for fac in patch.allFacilities:
        if fac.abbrev not in _TRACKED_FACILITIES_SET:
            continue
        facPPC = {}
        for ward in fac.getWards():
            pPC = ward.iA.getPatientPthCounts(timeNow)
            for k, v in pPC.items():
                if k in facPPC:
                    facPPC[k] += v
                else:
                    facPPC[k] = v
        for k, v in facPPC.items():
            facPthDict['%s_%s' % (fac.abbrev, k)] = v
    return facPthDict


def buildLocalTierPthDict(patch, timeNow):
    global _TRACKED_FACILITIES_SET
    if not _TRACKED_FACILITIES_SET:
        _TRACKED_FACILITIES_SET = frozenset(_TRACKED_FACILITIES)
    facPthDict = {'day': timeNow}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    for fac in patch.allFacilities:
        if fac.abbrev not in _TRACKED_FACILITIES_SET:
            continue
        facPPC = defaultdict(lambda: 0)
        for ward in fac.getWards():
            pPC = ward.iA.getPatientPthCounts(timeNow)
            for k, v in pPC.items():
                key = '%s_%s' % (ward.tier, k)
                facPPC[key] += v
        for key, v in facPPC.items():
            facPthDict['%s_%s' % (fac.abbrev, key)] = v  # so key is now abbrev_tier_pthStatus
    return facPthDict

# def buildLocalNColDict(patch,timeNow):
#     global _TRACKED_FACILITIES_SET
#     if not _TRACKED_FACILITIES_SET:
#         _TRACKED_FACILITIES_SET = frozenset(_TRACKED_FACILITIES)
#     facPthDict = {'day': timeNow}
#     assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
#     for fac in patch.allFacilities:
#         if fac.abbrev not in _TRACKED_FACILITIES_SET:
#             continue
#         facPPC = defaultdict(lambda: 0)
#         for ward in fac.getWards():
#             key = "{0}".format(fac.abbrev)
#             facPPC[key] += ward.newColonizationsSinceLastChecked
#     ### Need to figure out better way than to make sure that zeroing happens
#     return facPPC

def buildLocalTierNColDict(patch,timeNow):
    global _TRACKED_FACILITIES_SET
    if not _TRACKED_FACILITIES_SET:
        _TRACKED_FACILITIES_SET = frozenset(_TRACKED_FACILITIES)
    facPthDict = {'day': timeNow}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    for fac in patch.allFacilities:
        if fac.abbrev not in _TRACKED_FACILITIES_SET:
            continue
        facPPC = defaultdict(lambda: 0)
        facSum= 0.0
        for ward in fac.getWards():
            pPC = ward.miscCounters['newColonizationsSinceLastChecked']
            ward.miscCounters['newColonizationsSinceLastChecked'] = 0.0
            key = "{0}".format(ward.tier)
            facPPC[key] += pPC
            facSum += pPC
        for key, v in facPPC.items():
            facPthDict['{0}_{1}'.format(fac.abbrev,key)] = v
        #facPthDict['{0}'.format(fac.abbrev)] = facSum
    
    return facPthDict

def buildLocalTierCPDict(patch,timeNow):
    global _TRACKED_FACILITIES_SET
    if not _TRACKED_FACILITIES_SET:
        _TRACKED_FACILITIES_SET = frozenset(_TRACKED_FACILITIES)
    facTrtDict = {'day': timeNow}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    for fac in patch.allFacilities:
        if fac.abbrev not in _TRACKED_FACILITIES_SET:
            continue
        facPPC = defaultdict(lambda: 0)
        for ward in fac.getWards():
            key = "{0}".format(ward.tier)
            facPPC[key] += ward.miscCounters['patientDaysOnCP']
            ward.miscCounters['patientDaysOnCP'] = 0.0
        for key,v in facPPC.items():
            facTrtDict['{0}_{1}'.format(fac.abbrev,key)] = v  
            
    return facTrtDict
        
        
PER_DAY_NOTES_GEN_DICT = {'occupancy': buildFacOccupancyDict,
                          'localoccupancy': buildLocalOccupancyDict,
                          'pathogen': buildFacPthDict,
                          'localpathogen': buildLocalPthDict,
                          'localtierpathogen': buildLocalTierPthDict,
                          #'localnewcolonized': buildLocalNColDict,
                          'localtiernewcolonized':buildLocalTierNColDict,
                          'localtierCP': buildLocalTierCPDict
                          }


def loadPolicyImplementations(implementationDir):
    logger.info('Loading policy implementations')
    implList = []
    implementationDir = pyrheautils.pathTranslate(implementationDir)
    pyrheautils.PATH_STRING_MAP['POLICYDIR'] = implementationDir
    for newMod in pyrheautils.loadModulesFromDir(implementationDir,
                                                 requiredAttrList=['getPolicyClasses']):
        newPolicyClasses = newMod.getPolicyClasses()
        logger.info('Loaded %s implementing %s' % (newMod.__name__,
                                                   [cls.__name__
                                                    for cls in newPolicyClasses]))
        implList.extend(newPolicyClasses)
    return implList


def loadPathogenImplementations(implementationDir):
    logger.info('Loading infectious agent implementations')
    implementationDir = pyrheautils.pathTranslate(implementationDir)
    pyrheautils.PATH_STRING_MAP['PATHOGENDIR'] = implementationDir
    implDict = {}
    for newMod in pyrheautils.loadModulesFromDir(implementationDir,
                                                 requiredAttrList=['pathogenName',
                                                                   'getPathogenClass']):
        assert newMod.pathogenName not in implDict, \
            "Redundant definitions for pathogen %s" % newMod.category
        implDict[newMod.pathogenName] = newMod
        logger.info('Loaded %s implementing %s' % (newMod.__name__,
                                                   newMod.pathogenName))
    return implDict


def loadFacilityImplementations(implementationDir):
    logger.info('Loading facility implementations')
    # provide some string mapping
    implementationDir = pyrheautils.pathTranslate(implementationDir)
    pyrheautils.PATH_STRING_MAP['IMPLDIR'] = implementationDir
    implDict = {}
    for newMod in pyrheautils.loadModulesFromDir(implementationDir,
                                                 requiredAttrList=['category',
                                                                   'generateFull']):
        assert newMod.category not in implDict, \
            "Redundant definitions for category %s" % newMod.category
        implDict[newMod.category] = newMod
        logger.info('Loaded %s implementing %s' % (newMod.__name__, newMod.category))
    return implDict


def loadFacilityDescriptions(dirList, facImplDict, facImplRules):
    """
    This loads facilities, including communities, checking them against a schema.
    """
    logger.info('Parsing facilities descriptions')
    allKeySet = set()
    rawRecs = []
    for d in dirList[:]:
        d = pyrheautils.pathTranslate(d)
        kS, recs = yaml_tools.parse_all(d)
        rawRecs.extend(recs)
        allKeySet.update(kS)
    facRecs = []
    for rec in rawRecs:
        assert 'abbrev' in rec, "Facility description has no 'abbrev' field: %s" % rec
        assert 'category' in rec, \
            "Facility description for %(abbrev)s has no 'category' field" % rec
        facImplCategory = findFacImplCategory(facImplDict, facImplRules, rec['category'])
        if facImplCategory:
            facImpl = facImplDict[facImplCategory]
            if os.name != 'nt':
                nErrors = facImpl.checkSchema(rec)
                if nErrors:
                    logger.warning('dropping %s; %d schema violations' % (rec['abbrev'], nErrors))
                else:
                    facRecs.append(rec)
            else:
                facRecs.append(rec)
        else:
            raise RuntimeError('cannot find a facility implementation for %s, category %s',
                               rec['abbrev'], rec['category'])
    return facRecs


def distributeFacilities(comm, facDirs, facImplDict, facImplRules, partitionFile):
    """
    Estimate work associated with each facility and attempt to share it out fairly
    """
    if comm.size == 1:
        facRecs = loadFacilityDescriptions(facDirs, facImplDict, facImplRules)
        myFacList = facRecs
    else:
        assert partitionFile is not None, 'Expected a partition file'
        if comm.rank == 0:
            with open(partitionFile, 'rU') as f:
                parentDict = yaml.safe_load(f)
            partitionDict = parentDict['partition']
            facRecs = loadFacilityDescriptions(facDirs, facImplDict, facImplRules)
            assignments = defaultdict(list)
            for r in facRecs:
                assignments[partitionDict[r['abbrev']]].append(r)
            for targetRank in xrange(comm.size):
                if targetRank == comm.rank:
                    myFacList = assignments[targetRank]
                else:
                    comm.send(assignments[targetRank], dest=targetRank)
        else:
            myFacList = comm.recv(source=0)
        if comm.rank == 0:
            logger.info('Finished distributing facilities')
    return myFacList


def collectNotes(nhGroup, comm):
    allNotesGroup = noteholder.NoteHolderGroup()
    nhList = [allNotesGroup.copyNoteHolder(nh.getDict()) for nh in nhGroup.getnotes()]
    if comm.rank == 0:
        for targetRank in xrange(comm.size):
            if targetRank != comm.rank:
                dictList = comm.recv(source=targetRank)
                nhList.extend([allNotesGroup.copyNoteHolder(d) for d in dictList])
        return (allNotesGroup, nhList)
    else:
        dictList = [nh.getDict() for nh in nhList]
        comm.send(dictList, dest=0)
        return (None, None)


class TweakedOptParser(optparse.OptionParser):
    def setComm(self, comm):
        self.comm = comm

    def exit(self, code=None, msg=None):
        # print 'Here is your message: %s' % msg
        if code is None:
            code = -1
        if msg is None:
            msg = ""
        print msg
        # logger.critical(msg)
        if hasattr(self, 'comm') and self.comm.size > 1:
            self.comm.Abort(code)
        else:
            sys.exit(code)


def checkInputFileSchema(fname, schemaFname, comm=None):
    if logger is None:
        myLogger = logging.getLogger(__name__)
    else:
        myLogger = logger
    try:
        print "HERE"
        with open(fname, 'rU') as f:
            inputJSON = yaml.safe_load(f)
            if os.name != 'nt':
                validator = schemautils.getValidator(schemaFname)
                nErrors = sum([1 for e in validator.iter_errors(inputJSON)])  # @UnusedVariable
                if nErrors:
                    myLogger.error('Input file violates schema:')
                    for e in validator.iter_errors(inputJSON):
                        myLogger.error('Schema violation: %s: %s' %
                                       (' '.join([str(word) for word in e.path]), e.message))
                    if comm:
                        comm.Abort(2)
                    else:
                        raise RuntimeError('Input file violates schema')
                else:
                    return inputJSON
            else:
                return inputJSON
    except Exception, e:
        myLogger.error('Error checking input against its schema: %s' % e)
        if comm:
            comm.Abort(2)
        else:
            raise RuntimeError('Error checking input against its schema')




def createPerDayCB(patch, noteHolder, runDurationDays, recGenDict):
    """
    recGenDict is of the form {key: genFun} where genFun has the signature:
    
      dictOfNotes = genFun(patch, timeNow)
      
    The result is a set of notes of the form:
    
      key: [dictOfNotesAtTime0, dictOfNotesAtTime1, ..., dictOfNotesAtRunDurationDaysPlusOne]
    """
    def perDayCB(loop, timeNow):
#         for p in patchGroup.patches:
#             p.loop.printCensus()
        noteHolder.addNote({key: [recFun(patch, timeNow)]
                            for key, recFun in recGenDict.items()})
        if timeNow > runDurationDays:
            patch.group.stop()
    return perDayCB


def getLoggerConfig():
    """
    This routine reads and parses logger_cfg.yaml in the current directory and returns the
    result as a dict.
    """
    cfgFname = 'log_cfg.yaml'
    try:
        with open(cfgFname, 'rU') as f:
            cfgDict = yaml.safe_load(f)
    except IOError:
        print 'Cannot find %s - using default logging' % cfgFname
        cfgDict = {'version': 1}
    except Exception, e:
        print 'Failed to read or parse %s: %s' % (cfgFname, e)
        cfgDict = {'version': 1}
    return cfgDict


def configureLogging(cfgDict, cfgExtra):
    global logger
    logging.config.dictConfig(cfgDict)
    if cfgExtra is not None:
        logging.getLogger('root').setLevel(cfgExtra)
    logger = logging.getLogger(__name__)


class RankLoggingFilter(logging.Filter):
    def __init__(self, **kwargs):
        logging.Filter.__init__(self)
        self.rank = patches.getCommWorld().rank

    def filter(self, record):
        record.rank = self.rank
        return True


class PikaLogHandler(logging.Handler):

    class Producer(object):
        def __init__(self):
            conn = pika.BlockingConnection(pika.ConnectionParameters(host='localhost'))
            self.ch = conn.channel()
            self.ch.queue_declare(queue='pyrhea_logger')

        def message(self, message):
            self.ch.basic_publish(exchange='', routing_key='pyrhea_logger',
                                  body=pickle.dumps(message))

    def __init__(self):
        logging.Handler.__init__(self)
        self.broadcaster = PikaLogHandler.Producer()
        self.machine = os.uname()[1]

    def emit(self, record):
        message = {'source': 'logger', 'machine': self.machine,
                   'message': self.format(record), 'level': record.levelname,
                   'pathname': record.pathname, 'lineno': record.lineno,
                   'exception': record.exc_info, 'rank': record.rank}
        self.broadcaster.message(message)


class BurnInAgent(patches.Agent):
    def __init__(self, name, patch, burnInDays, noteHolderGroup):
        patches.Agent.__init__(self, name, patch)
        self.burnInDays = burnInDays
        self.noteHolderGroup = noteHolderGroup

    def run(self, startTime):
        timeNow = self.sleep(self.burnInDays)  # @UnusedVariable
        logger.info('burnInAgent is running')
        self.noteHolderGroup.clearAll(keepRegex='.*(([Nn]ame)|([Tt]ype)|([Cc]ode)|(_vol)|(occupancy)|(pathogen))')
        # and now the agent exits


class ScenarioStartAgent(patches.Agent):
    def __init__(self, name, patch, totalWaitDays, scenarioPolicyList):
        patches.Agent.__init__(self, name, patch)
        self.totalWaitDays = totalWaitDays
        self.scenarioPolicyList = scenarioPolicyList

    def run(self, startTime):
        timeNow = self.sleep(self.totalWaitDays)  # @UnusedVariable
        logger.info('Scenario is beginning')
        for sP in self.scenarioPolicyList:
            sP.begin(self, timeNow)
        # and now the agent exits


def findPolicies(policyClassList,
                 policyRules,
                 category,
                 abbrev = None):
    l = []
    abbrevFlag = False
    for pCl in policyClassList:
        for categoryRegex, classRegex in policyRules:
            if abbrev is not None and categoryRegex.match(abbrev) and classRegex.match(pCl.__name__):
                l.append(pCl)
                abbrevFlag = True
            else:
                if not abbrevFlag:
                    if categoryRegex.match(category) and classRegex.match(pCl.__name__):
                        l.append(pCl)
         
    return l


def findFacImplCategory(facImplDict,
                        facImplRules,
                        category):
    """
    The 'category' parameter comes from a facility description 'category' attribute.
    Map it to its matching facility implementation 'category' attribute using the
    supplied rules.
    """
    for facRegex, implStr in facImplRules:
        if facRegex.match(category):
            return implStr
    return None


def createFacImplMap(facImplDict, facImplRules):
    def myMap(category):
        val = findFacImplCategory(facImplDict, facImplRules, category)
        if val:
            return val
        else:
            raise RuntimeError('cannot find a facility implementation for %s' % category)
    return myMap


def initializeFacilities(patchList, myFacList, facImplDict, facImplRules,
                         policyClassList, policyRules,
                         PthClass, noteHolderGroup, comm, totalRunDays):
    """Distribute facilities across patches and initialize them"""
    offset = 0
    tupleList = [(p, [], [], []) for p in patchList]
    for facDescription in myFacList:
        
        patch, allIter, allAgents, allFacilities = tupleList[offset]
        facImplCategory = findFacImplCategory(facImplDict, facImplRules,
                                              facDescription['category'])
        if facImplCategory:
            facImpl = facImplDict[facImplCategory]
            facImplMapFun = createFacImplMap(facImplDict, facImplRules)
            facilities, wards, patients = \
                facImpl.generateFull(facDescription, patch,
                                     policyClasses=findPolicies(policyClassList,
                                                                policyRules,
                                                                facImpl.category,
                                                                facDescription['abbrev']),
                                     categoryNameMapper=facImplMapFun)
            
            for w in wards:
                w.setInfectiousAgent(PthClass(w, useWardCategory=facImplCategory))
                w.initializePatientPthState()
                w.initializePatientTreatment()
            for fac in facilities:
                nh = noteHolderGroup.createNoteHolder()
                nh.addNote({'rank': comm.rank})
                fac.setNoteHolder(nh)
            for rQList in [fac.reqQueues for fac in facilities]:
                allIter.extend(rQList)
            allIter.extend([fac.holdQueue for fac in facilities])
            allIter.extend(wards)
            allAgents.extend([fac.manager for fac in facilities])
            allAgents.extend(patients)
            allFacilities.extend(facilities)
        else:
            raise RuntimeError('Facility %(abbrev)s category %(category)s has no implementation' %
                               facDescription)
        offset = (offset + 1) % len(tupleList)

    for patch, allIter, allAgents, allFacilities in tupleList:
        patch.addInteractants(allIter)
        patch.addAgents(allAgents)
        patch.allFacilities = allFacilities
        logger.info('Rank %d patch %s: %d interactants, %d agents, %d facilities' %
                    (comm.rank, patch.name, len(allIter), len(allAgents), len(allFacilities)))
        patchNH = noteHolderGroup.createNoteHolder()
        patch.loop.addPerDayCallback(createPerDayCB(patch, patchNH, totalRunDays,
                                                    recGenDict=PER_DAY_NOTES_GEN_DICT))
        initialNote = {'name': patch.name, 'rank': comm.rank}
        for key, genFun in PER_DAY_NOTES_GEN_DICT.items():
            initialNote[key] = [genFun(patch, 0)]
        patchNH.addNote(initialNote)


def main():

    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    if os.name != "nt":
        signal.signal(signal.SIGUSR1, handle_pdb)

    comm = patches.getCommWorld()
    logging.basicConfig(format="[%d]%%(levelname)s:%%(name)s:%%(message)s" % comm.rank)

    if comm.rank == 0:
        parser = TweakedOptParser(usage="""
        %prog [-v][-d][-t][-D][-p=npatch][-C][-L=loglevel][-P=partitionfile.yaml][--seed SEED] input.yaml
        """)
        parser.setComm(comm)
        parser.add_option("-v", "--verbose", action="store_true",
                          help="verbose output")
        parser.add_option("-d", "--debug", action="store_true",
                          help="debugging output")
        parser.add_option("-t", "--trace", action="store_true",
                          help="enable greenlet thread tracing")
        parser.add_option("-D", "--deterministic", action="store_true",
                          help="deterministic mode")
        parser.add_option("-p", "--patches", action="store", type="int",
                          help="Patches per MPI rank (default 1)", default=1)
        parser.add_option("-C", "--census", action="store_true",
                          help="print census for each rank every tick")
        parser.add_option("-L", "--log", action="store", type="string",
                          help=("Set logging level "
                                "('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL')"),
                          default=None)
        parser.add_option("-o", "--out", action="store", type="string",
                          help=("Set output file name (default '%s')" %
                                DEFAULT_OUTPUT_NOTES_NAME))
        parser.add_option("-P", "--partition", action="store", type="string",
                          help=("yaml file defining the partition of locations to ranks"
                                " (no default)"))
        parser.add_option("--seed", action="store", type="int",
                          help="Use this value as the random seed")

        opts, args = parser.parse_args()
        if opts.log is not None:
            numLogLevel = getattr(logging, opts.log.upper(), None)
            if not isinstance(numLogLevel, int):
                parser.error("Invalid log level: %s" % opts.log)
        else:
            numLogLevel = None
        if opts.partition is None and comm.size > 1:
            parser.error('A partition file is required for parallel runs')
        clData = {'verbose': opts.verbose,
                  'debug': opts.debug,
                  'trace': opts.trace,
                  'deterministic': opts.deterministic,
                  'patches': opts.patches,
                  'printCensus': opts.census,
                  'logCfgDict': getLoggerConfig(),
                  'loggingExtra': numLogLevel,
                  'partitionFile': opts.partition,
                  'randomSeed': opts.seed
                  }
        if len(args) == 1:
            clData['input'] = checkInputFileSchema(args[0],
                                                   os.path.join(SCHEMA_DIR, INPUT_SCHEMA),
                                                   comm)
        else:
            parser.error("A YAML-format file specifying run parameters must be specified.")

        # NOTE- outputNotesName is defined only for rank 0.
        if opts.out:
            outputNotesName = opts.out
        elif 'notesFileName' in clData['input']:
            outputNotesName = clData['input']['notesFileName']
        else:
            outputNotesName = DEFAULT_OUTPUT_NOTES_NAME
        parser.destroy()
    else:
        clData = None
        
    try:
        patchGroup = None  # for the benefit of diagnostics during init
        clData = comm.bcast(clData, root=0)
        if 'modelDir' in clData['input']:
            pyrheautils.PATH_STRING_MAP['MODELDIR'] = os.path.abspath(clData['input']['modelDir'])
        if 'pathTranslations' in clData['input']:
            for elt in clData['input']['pathTranslations']:
                pyrheautils.PATH_STRING_MAP[elt['key']] = elt['value']
        configureLogging(clData['logCfgDict'], clData['loggingExtra'])
    
        verbose = clData['verbose']  # @UnusedVariable
        debug = clData['debug']  # @UnusedVariable
        trace = clData['trace']
        deterministic = clData['deterministic']
        inputDict = clData['input']
    
        if deterministic:
            seed(1234 + comm.rank)  # Set the random number generator seed
        elif clData['randomSeed']:
            seed(clData['randomSeed'] + comm.rank)
        elif 'randomSeed' in inputDict:
            seed(inputDict['randomSeed'] + comm.rank)

        if 'trackedFacilities' in inputDict:
            global _TRACKED_FACILITIES
            if _TRACKED_FACILITIES:
                _TRACKED_FACILITIES.extend(inputDict['trackedFacilities'])
            else:
                _TRACKED_FACILITIES = inputDict['trackedFacilities'][:]

        schemautils.setSchemaBasePath(SCHEMA_DIR)
        pthImplDict = loadPathogenImplementations(inputDict['pathogenImplementationDir'])
        assert len(pthImplDict) == 1, 'Simulation currently supports exactly one pathogen'
        PthClass = pthImplDict.values()[0].getPathogenClass()
    
        facImplDict = loadFacilityImplementations(inputDict['facilityImplementationDir'])
        if 'facilitySelectors' in inputDict:
            facImplRules = [(re.compile(rule['category']), rule['implementation'])
                            for rule in inputDict['facilitySelectors']]
        else:
            facImplRules = [(re.compile(category), category)
                            for category in facImplDict.keys()]  # an identity map
    
        policyClassList = loadPolicyImplementations(inputDict['policyImplementationDir'])
        policyRules = [(re.compile(rule['category']), re.compile(rule['policyClass']))
                       for rule in inputDict['policySelectors'] if 'category' in rule.keys() ]
        policyRules += [(re.compile(rule['locationAbbrev']),re.compile(rule['policyClass']))
                           for rule in inputDict['policySelectors'] if 'locationAbbrev' in rule.keys() ]
        
        myFacList = distributeFacilities(comm, inputDict['facilityDirs'],
                                         facImplDict, facImplRules,
                                         clData['partitionFile'])
        logger.info('Rank %d has facilities %s' % (comm.rank, [rec['abbrev'] for rec in myFacList]))
    
        patchGroup = patches.PatchGroup(comm, trace=trace, deterministic=deterministic,
                                        printCensus=clData['printCensus'])
        noteHolderGroup = noteholder.NoteHolderGroup()
        patchList = [patchGroup.addPatch(patches.Patch(patchGroup))
                     for i in xrange(clData['patches'])]  # @UnusedVariable
        totalRunDays = inputDict['runDurationDays'] + inputDict['burnInDays']
    
        # Add some things for which we need only one instance
        patchList[0].addAgents([BurnInAgent('burnInAgent', patchList[0],
                                            inputDict['burnInDays'], noteHolderGroup)])
        if 'scenarioWaitDays' in inputDict:
            scenarioPolicyClasses = []
            policyClasses = (findPolicies(policyClassList, policyRules, 'scenario')
                             or [])
            scenarioPolicyClasses = [pC for pC in policyClasses
                                     if issubclass(pC, ScenarioPolicy)]
            if not scenarioPolicyClasses:
                scenarioPolicyClasses = [ScenarioPolicy]
            scenarioPolicies = [sPC('scenario_%d' % ind, patchList[0])
                                for ind, sPC in enumerate(scenarioPolicyClasses)]
            ssA = ScenarioStartAgent('scenarioStartAgent', patchList[0],
                                     inputDict['burnInDays'] + inputDict['scenarioWaitDays'],
                                     scenarioPolicies)
            patchList[0].addAgents([ssA])
    
        initializeFacilities(patchList, myFacList, facImplDict, facImplRules,
                             policyClassList, policyRules,
                             PthClass, noteHolderGroup, comm, totalRunDays)
    except Exception, e:
        if patchGroup:
            logger.error('%s exception during initialization: %s; traceback follows' %
                         (patchGroup.name, e))
            import traceback
            traceback.print_exc(file=sys.stderr)
            logger.info('%s exiting due to exception during initialization' % patchGroup.name)
            logging.shutdown()
            sys.exit('Exception during initialization')
        else:
            logger.error('%s exception during initialization: %s; traceback follows' %
                         ('<no patchGroup yet>', e))
            import traceback
            traceback.print_exc(file=sys.stderr)
            logger.info('%s exiting due to exception during initialization' % '<no patchGroup yet>')
            logging.shutdown()
            sys.exit('Exception during initialization')
    
    logger.info('%s #### Ready to Run #### (from main)' % patchGroup.name)
    try:
        exitMsg = patchGroup.start()
        logger.info('%s #### all done #### (from main); %s' % (patchGroup.name, exitMsg))
    except Exception, e:
        logger.error('%s exception %s; traceback follows' % (patchGroup.name, e))
        import traceback
        traceback.print_exc(file=sys.stderr)
    finally:
        logger.info('%s writing notes and exiting' % patchGroup.name)

        allNotesGroup, allNotesList = collectNotes(noteHolderGroup, comm)  # @UnusedVariable
        if comm.rank == 0:
            d = {}
            for nh in allNotesGroup.getnotes():
                d[nh['name']] = nh.getDict()
    #             if 'occupancy' in nh:
    #                 recs = nh['occupancy']
    #                 with open(('occupancy_%s.csv' % nh['name']), 'w') as f:
    #                     csv_tools.writeCSV(f, recs[1].keys(), recs)
            with open(outputNotesName, 'w') as f:
                pickle.dump(d, f)

        logging.shutdown()

############
# Main hook
############

if __name__ == "__main__":
    main()
