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
import sys
import os
import logging
import logging.config
from random import seed
from collections import defaultdict
import optparse
import re
import signal
import types

import yaml
import pika
import six.moves.cPickle as pickle
import ujson

import quilt.patches as patches
import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.notes.noteholder as noteholder
import schemautils
import pyrheautils
from typebase import PatientOverallHealth
from registry import Registry
from policybase import ScenarioPolicy
from tauadjuster import TauAdjuster
import checkpoint
import bcz_monitor
from closuretricks import ClosureFixer
from daydata import DayDataGroup

BASE_DIR = os.path.dirname(__file__)
SCHEMA_DIR = os.path.join(BASE_DIR, '../schemata')
INPUT_SCHEMA = 'rhea_input_schema.yaml'
PTH_IMPLEMENTATIONS_DIR = os.path.join(BASE_DIR, 'pathogenImplementations/$(PATHOGEN)')
DEFAULT_OUTPUT_NOTES_NAME = 'notes.json'

_TRACKED_FACILITIES = []

LOGGER = None

# make command line data accessible to anyone
CL_DATA = None


def getTrackedFacilitiesSet():
    if ClosureFixer.TRACKED_FACILITIES_SET is None:
        ClosureFixer.TRACKED_FACILITIES_SET = frozenset(_TRACKED_FACILITIES)
    return ClosureFixer.TRACKED_FACILITIES_SET


def buildFacOccupancyDict(patch, timeNow):
    facTypeDict = {'day': timeNow, 'format': ['factype']}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    for fac in patch.allFacilities:
        tpName = type(fac).__name__
        if tpName not in facTypeDict:
            facTypeDict[tpName] = 0
            facTypeDict[tpName + '_all'] = 0

        if hasattr(fac, 'patientStats'):
            patientCount = fac.patientStats.currentOccupancy
            allCount = fac.patientStats.totalOccupancy
            facTypeDict[tpName] += patientCount
            facTypeDict[tpName + '_all'] += allCount

    #print('buildFacOccupancyDict: %s') % facTypeDict
    return facTypeDict


def buildFacHeldBedDict(patch, timeNow):
    facTypeDict = defaultdict(int)
    facTypeDict['day'] = timeNow
    facTypeDict['format'] = ['factype', 'bedallockey']
    assert hasattr(patch, 'allFacilities'), ('patch %s has no list of facilities!'
                                             % patch.name)
    for fac in patch.allFacilities:
        tpName = type(fac).__name__
        if hasattr(fac, 'bedAllocDict'):
            for key, ct in fac.bedAllocDict.items():
                facTypeDict[(tpName, key)] += ct

    return facTypeDict


def buildFacOverallHealthDict(patch, timeNow):
    facOHDict = {'day': timeNow, 'format': ['factype', 'ohidx']}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    seenFacTypes = []
    for fac in patch.allFacilities:
        ctD = defaultdict(lambda: 0)
        for ward in fac.getWards():
            for pOH in PatientOverallHealth.names:
                ct = ward.cumStats.popByOH(pOH)
                if ct != 0:  # COMMUNITY pops can actually be negative due to initialization issues
                    ctD[pOH] += ct
        tpName = type(fac).__name__
        if tpName not in seenFacTypes:
            for ohIdx in PatientOverallHealth.names:
                key = (tpName, ohIdx)
                facOHDict[key] = 0
            seenFacTypes.append(tpName)
        for pOH, ct in ctD.items():
            key = (tpName, pOH)
            assert key in facOHDict, 'No key %s?' % key
            facOHDict[key] += ct
    return facOHDict


def buildFacPthDict(patch, timeNow):
    facPthDict = {'day': timeNow, 'format': ['factype', 'pthidx']}
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
            key = (tpName, k)
            if key in facPthDict:
                facPthDict[key] += v
            else:
                facPthDict[key] = v
    return facPthDict


def buildLocalOccupancyDict(patch, timeNow):
    facDict = {'day': timeNow, 'format': ['fac']}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    for fac in patch.allFacilities:
        if hasattr(fac, 'patientStats'):
            patientCount = fac.patientStats.currentOccupancy
            facDict[fac.abbrev] = patientCount
    return facDict


def buildLocalPthDict(patch, timeNow):
    facPthDict = {'day': timeNow, 'format': ['fac', 'pthidx']}
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
        for k, v in facPPC.items():
            facPthDict[(fac.abbrev, k)] = v
    return facPthDict


def buildLocalTierPthDict(patch, timeNow):
    facPthDict = {'day': timeNow, 'format': ['fac', 'tieridx, pthidx']}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    for fac in patch.allFacilities:
        facPPC = defaultdict(lambda: 0)
        for ward in fac.getWards():
            pPC = ward.iA.getPatientPthCounts(timeNow)
            for k, v in pPC.items():
                key = (ward.tier, k)
                facPPC[key] += v
        for key, v in facPPC.items():
            wT, k = key
            facPthDict[(fac.abbrev, wT, k)] = v  # so key is now abbrev_tier_pthStatus
    return facPthDict

def generateLocalTierDictBuilder(counterKey):
    def buildDict(patch, timeNow):
        facTrtDict = {'day': timeNow, 'format': ['fac', 'tieridx', 'wardidx']}
        assert hasattr(patch, 'allFacilities'), ('patch %s has no list of facilities!'
                                                 % patch.name)
        for fac in patch.allFacilities:
            facPPC = defaultdict(lambda: 0)
            for ward in fac.getWards():
                key = (ward.tier, ward.wardNum)
                facPPC[key] = ward.miscCounters[counterKey]
                ward.miscCounters[counterKey] = 0.0
            for (tN, wN), v in facPPC.items():
                facTrtDict[(fac.abbrev, tN, wN)] = v
        return facTrtDict
    return buildDict

def generateFacTypeDictBuilder(counterKey):
    def buildDict(patch, timeNow):
        facTypeDict = {'day': timeNow, 'format': ['factype']}
        assert hasattr(patch, 'allFacilities'), ('patch %s has no list of facilities!'
                                                 % patch.name)
        for fac in patch.allFacilities:
            tpName = type(fac).__name__
            if tpName not in facTypeDict:
                facTypeDict[tpName] = 0
            for ward in fac.getWards():
                facTypeDict[tpName] += ward.miscCounters[counterKey]
                ward.miscCounters[counterKey] = 0.0
        return facTypeDict
    return buildDict

def tupleToNoteKey(tpl):
    if isinstance(tpl, types.TupleType):
        return '_'.join([str(elt) for elt in tpl])
    else:
        return tpl  # which isn't a tuple

def dictToNote(dct):
    """
    This takes a dict in the format produced by the daily notes routines
    and converts it to the format expected of the Notes output.
    """
    if 'format' in dct and dct['format'][-1] == 'wardidx':
        innerD = defaultdict(int)
        for key, val in dct.items():
            if isinstance(key, types.TupleType):
                innerKey = tuple(key[:-1])
                innerD[innerKey] += val
            elif key == 'format':
                innerD['format'] = val[:-1]
            else:
                innerD[key] = val
        return dictToNote(innerD)
    else:
        return {tupleToNoteKey(key): val for key, val in dct.items() if key != 'format'}

PER_DAY_NOTES_GEN_DICT = {'occupancy': buildFacOccupancyDict,
                          'localoccupancy': buildLocalOccupancyDict,
                          'pathogen': buildFacPthDict,
                          'occupancyByOH': buildFacOverallHealthDict,
                          'localpathogen': buildLocalPthDict,
                          'localtierpathogen': buildLocalTierPthDict,
                          'localtiernewcolonized': generateLocalTierDictBuilder('newColonizationsSinceLastChecked'),
                          'localtierCP': generateLocalTierDictBuilder('patientDaysOnCP'),
                          'localtierCREBundle':  generateLocalTierDictBuilder('creBundlesHandedOut'),
                          'localtierCRESwabs': generateLocalTierDictBuilder('creSwabsPerformed'),
                          'localtierarrivals': generateLocalTierDictBuilder('arrivals'),
                          'localtierdepartures': generateLocalTierDictBuilder('departures'),
                          'localtiercrearrivals': generateLocalTierDictBuilder('creArrivals'),
                          'localtierpassiveCP': generateLocalTierDictBuilder('passiveDaysOnCP'),
                          'localtierswabCP': generateLocalTierDictBuilder('swabDaysOnCP'),
                          'localtierotherCP': generateLocalTierDictBuilder('otherDaysOnCP'),
                          'localtierxdroCP': generateLocalTierDictBuilder('xdroDaysOnCP'),
                          'localtierpatientsOnCP': generateLocalTierDictBuilder('newPatientsOnCP'),
                          'localtiernthawed': generateFacTypeDictBuilder('nThawed'),
                          'bedHoldStats': buildFacHeldBedDict
                          }


for key, genFun in PER_DAY_NOTES_GEN_DICT.items():
    DayDataGroup().add(key, genFun)


def defineTrackingGroup(gpName, keyTypeL):
    def extractFun(ward, timeNow):
        patch = ward.patch
        rslt = []
        for key, tp in keyTypeL:  # @UnusedVariable
            dct = DayDataGroup().get(patch, key, timeNow)
            assert dct['format'] == ['fac', 'tieridx', 'wardidx'], (('Cannot build tracking group'
                                                                     ' including %s; wrong format')
                                                                     % key)
            assert dct['day'] == timeNow, (('Tracked quantity entry generator returned wrong date'
                                            ' day (%s vs %s)')
                                           % (dct['day'], timeNow))
            tplKey = (ward.fac.abbrev, ward.tier, ward.wardNum)
            rslt.append(dct[tplKey] if tplKey in dct else 0)
        return rslt
    bcz_monitor.addTrackable(gpName, extractFun, keyTypeL)

defineTrackingGroup('arrivals', [('localtierarrivals', 'uint32'), ('localtierdepartures', 'uint32')])

defineTrackingGroup('creCounters', [('localtierCP', 'uint32'),
                          ('localtierCREBundle', 'uint32'),
                          ('localtierCRESwabs', 'uint32'),
                          ('localtiercrearrivals', 'uint32'),
                          ('localtierpassiveCP', 'uint32'),
                          ('localtierswabCP', 'uint32'),
                          ('localtierotherCP', 'uint32'),
                          ('localtierxdroCP', 'uint32'),
                          ('localtierpatientsOnCP', 'uint32'),
                          ])

defineTrackingGroup('newColonizations', [('localtiernewcolonized', 'uint32')])


def dumpFacilitiesMap(filename, patchList):
    from facilitybase import CareTier
    from collections import defaultdict
    tau = {}
    fracColonized = {}
    category = {}
    bedcounts = defaultdict(int)
    wards = defaultdict(int)
    startPop = defaultdict(int)

    for patch in patchList:
        for fac in patch.allFacilities:
            abbr = fac.abbrev
            for ward in fac.getWards():
                if ward.tier == CareTier.HOME:
                    continue
                tier = CareTier.names[ward.tier]
                key = (abbr,tier)

                category[key] = fac.category
                bedcounts[key] += ward._nLocks  # must I really go this low level to get this info?
                wards[key] += 1
                startPop[key] += ward.cumStats.popTotal()

                try: tau[key] = ward.iA.tau
                except: tau[key] = ""
                try: fracColonized[key] = ward.iA.initialFracColonized
                except: fracColonized[key] = ""


    with open(filename, "w") as f:
        f.write("fac,tier,category,wards,beds,startPop,tau,fracColonized\n")
        for k in category.keys():
            abbr,tier = k
            f.write("%s,%s,%s,%s,%s,%s,%s,%s\n"%(abbr,tier,category[k],wards[k],
                                                 bedcounts[k],startPop[k],
                                                 tau[k],fracColonized[k]))




def loadPolicyImplementations(implementationDir):
    LOGGER.info('Loading policy implementations')
    implList = []
    implementationDir = pyrheautils.pathTranslate(implementationDir)
    for newMod in pyrheautils.loadModulesFromDir(implementationDir,
                                                 requiredAttrList=['getPolicyClasses']):
        newPolicyClasses = newMod.getPolicyClasses()
        LOGGER.info('Loaded %s implementing %s' % (newMod.__name__,
                                                   [cls.__name__
                                                    for cls in newPolicyClasses]))
        implList.extend(newPolicyClasses)
    return implList


def loadPathogenImplementations(implementationDir):
    LOGGER.info('Loading infectious agent implementations')
    implementationDir = pyrheautils.pathTranslate(implementationDir)
    implDict = {}
    for newMod in pyrheautils.loadModulesFromDir(implementationDir,
                                                 requiredAttrList=['pathogenName',
                                                                   'getPathogenClass']):
        assert newMod.pathogenName not in implDict, \
            "Redundant definitions for pathogen %s" % newMod.category
        implDict[newMod.pathogenName] = newMod
        LOGGER.info('Loaded %s implementing %s' % (newMod.__name__,
                                                   newMod.pathogenName))
    return implDict


def loadFacilityImplementations(implementationDir):
    LOGGER.info('Loading facility implementations')
    # provide some string mapping
    implementationDir = pyrheautils.pathTranslate(implementationDir)
    implDict = {}
    for newMod in pyrheautils.loadModulesFromDir(implementationDir,
                                                 requiredAttrList=['category',
                                                                   'generateFull']):
        assert newMod.category not in implDict, \
            "Redundant definitions for category %s" % newMod.category
        implDict[newMod.category] = newMod
        LOGGER.info('Loaded %s implementing %s' % (newMod.__name__, newMod.category))
    return implDict


def loadFacilityDescriptions(dirList, facImplDict, facImplRules):
    """
    This loads facilities, including communities, checking them against a schema.
    """
    LOGGER.info('Parsing facilities descriptions')
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

        if rec['abbrev'] in pyrheautils.facilitiesReplacementData:
            rec = pyrheautils.replaceYamlData(pyrheautils.facilitiesReplacementData[rec['abbrev']], rec)
        
        assert 'category' in rec, \
            "Facility description for %(abbrev)s has no 'category' field" % rec
        facImplCategory = findFacImplCategory(facImplDict, facImplRules, rec['category'])
        if facImplCategory:
            facImpl = facImplDict[facImplCategory]
            if os.name != 'nt':
                nErrors = facImpl.checkSchema(rec)
                if nErrors:
                    LOGGER.warning('dropping %s; %d schema violations' % (rec['abbrev'], nErrors))
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
            LOGGER.info('Finished distributing facilities')
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
        print(msg)
        # LOGGER.critical(msg)
        if hasattr(self, 'comm') and self.comm.size > 1:
            self.comm.Abort(code)
        else:
            sys.exit(code)


def checkInputFileSchema(fname, schemaFname, comm=None):
    if LOGGER is None:
        myLogger = logging.getLogger(__name__)
    else:
        myLogger = LOGGER
    try:
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
    except Exception as e:
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
        """ This function is called once per day for maintenance """
#         for p in patchGroup.patches:
#             p.loop.printCensus()
        tFS = getTrackedFacilitiesSet()
        dctD = {}
        for key in recGenDict:
            dctD = {key: [dictToNote(DayDataGroup().get(patch, key, timeNow))]
                    for key in recGenDict}
            dayD = DayDataGroup().get(patch, key, timeNow)
            if dayD['format'][0] == 'fac':
                # Drop data for facilities that are not of interest
                grpD = {dKey: dVal for dKey, dVal in dayD.items() if dKey[0] in tFS}
            else:
                grpD = dayD.copy()
            dctD[key] = [dictToNote(grpD)]
        noteHolder.addNote(dctD)
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
        print('Cannot find %s - using default logging' % cfgFname)
        cfgDict = {'version': 1}
    except Exception as e:
        print('Failed to read or parse %s: %s' % (cfgFname, e))
        cfgDict = {'version': 1}
    return cfgDict


def configureLogging(cfgDict, cfgExtra):
    global LOGGER
    logging.config.dictConfig(cfgDict)
    if cfgExtra is not None:
        logging.getLogger('root').setLevel(cfgExtra)
    LOGGER = logging.getLogger(__name__)


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
        LOGGER.info('burnInAgent is running at time %s', timeNow)
        self.noteHolderGroup.clearAll(keepRegex='.*(([Nn]ame)|([Tt]ype)|([Cc]ode)|(_vol)|(occupancy)|(pathogen))')
        # and now the agent exits


class ScenarioStartAgent(patches.Agent):
    def __init__(self, name, patch, totalWaitDays, scenarioPolicyList):
        patches.Agent.__init__(self, name, patch)
        self.totalWaitDays = totalWaitDays
        self.scenarioPolicyList = scenarioPolicyList

    def run(self, startTime):
        timeNow = self.sleep(self.totalWaitDays)  # @UnusedVariable
        LOGGER.info('Scenario is beginning')
        for sP in self.scenarioPolicyList:
            sP.begin(self, timeNow)
        # and now the agent exits


def findPolicies(policyClassList,
                 policyRulesDict,
                 category,
                 abbrev = None):
    l = []
    for pCl in policyClassList:
        for ruleKey in policyRulesDict.keys():
            categoryRegex, classRegex = ruleKey[0:2]
            if ((abbrev is not None and categoryRegex.match(abbrev))
                    or categoryRegex.match(category)) and classRegex.match(pCl.__name__):
                l.append(pCl)
                policyRulesDict[ruleKey] = True  # rule has been used
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
                         policyClassList, policyRulesDict,
                         PthClass, noteHolderGroup, comm, totalRunDays):
    """Distribute facilities across patches and initialize them"""
    offset = 0
    tupleList = [(p, [], [], []) for p in patchList]

    # Every patch gets one Registry instance
    for patch, allIter, allAgents, allFacilities in tupleList:
        registry = Registry(patch.name + '_Registry', patch)
        registry.setNoteHolder(noteHolderGroup.createNoteHolder())
        allIter.extend(registry.reqQueues)
        allIter.append(registry.holdQueue)
        allAgents.append(registry.manager)

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
                                                                policyRulesDict,
                                                                facImpl.category,
                                                                facDescription['abbrev']),
                                     categoryNameMapper=facImplMapFun)

            for w in wards:
                w.setInfectiousAgent(PthClass(w, implCategory=facImplCategory))
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
            for fac in facilities:
                fac.finalizeBuild(facDescription)
            allFacilities.extend(facilities)
        else:
            raise RuntimeError('Facility %(abbrev)s category %(category)s has no implementation' %
                               facDescription)
        offset = (offset + 1) % len(tupleList)

    for patch, allIter, allAgents, allFacilities in tupleList:
        patch.addInteractants(allIter)
        patch.addAgents(allAgents)
        patch.allFacilities = allFacilities
        LOGGER.info('Rank %d patch %s: %d interactants, %d agents, %d facilities' %
                    (comm.rank, patch.name, len(allIter), len(allAgents), len(allFacilities)))
        patchNH = noteHolderGroup.createNoteHolder()
        patch.loop.addPerDayCallback(createPerDayCB(patch, patchNH, totalRunDays,
                                                    recGenDict=PER_DAY_NOTES_GEN_DICT))
        initialNote = {'name': patch.name, 'rank': comm.rank}
        for key, genFun in PER_DAY_NOTES_GEN_DICT.items():
            initialNote[key] = [dictToNote(genFun(patch, 0))]
        patchNH.addNote(initialNote)


def main():

    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)

    global LOGGER
    global CL_DATA

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
        parser.add_option("-k", "--checkpoint", action="store", type="int", default=-1,
                          help="save the state and stop after n ticks")
        parser.add_option("-c", "--constantsFile", action="store", type="string", default=None,
                          help="python file defining the dict constantsReplacementData and/or facilitiesReplacementData")
        parser.add_option("--saveNewConstants", action="store", type="string", default=None,
                          help="save any modified constants files (from -c) to the specified directory")
        parser.add_option("-b", "--bczmonitor", action="store", type="string", default=None,
                          help="save pathogen status as a pandas data structure in the file specified")
        parser.add_option("--taumod", action="store_true", default=False,
                          help="run pyrhea in the taumod mode")
        parser.add_option("-n", "--disableNotes", action="store_true",
                          help="disable noteholder functions to save memory (a minimal notes file will still be written)")
        parser.add_option("-m", "--dumpFacilitiesMap", action="store", type="string", default=None,
                          help="write a facililties map to the file specified to facilitate post processing")

        opts, args = parser.parse_args()
        if opts.log is not None:
            numLogLevel = getattr(logging, opts.log.upper(), None)
            if not isinstance(numLogLevel, int):
                parser.error("Invalid log level: %s" % opts.log)
        else:
            numLogLevel = None
        if opts.partition is None and comm.size > 1:
            parser.error('A partition file is required for parallel runs')
        CL_DATA = {'verbose': opts.verbose,
                   'debug': opts.debug,
                   'trace': opts.trace,
                   'deterministic': opts.deterministic,
                   'patches': opts.patches,
                   'printCensus': opts.census,
                   'logCfgDict': getLoggerConfig(),
                   'loggingExtra': numLogLevel,
                   'partitionFile': opts.partition,
                   'randomSeed': opts.seed,
                   'checkpoint': opts.checkpoint,
                   'constantsFile': opts.constantsFile,
                   'saveNewConstants': opts.saveNewConstants,
                   'bczmonitor': opts.bczmonitor,
                   'taumod': opts.taumod,
                   'dumpFacilitiesMap': opts.dumpFacilitiesMap,
                   'disableNotes' : opts.disableNotes,
        }
        if len(args) == 1:
            CL_DATA['input'] = checkInputFileSchema(args[0],
                                                   os.path.join(SCHEMA_DIR, INPUT_SCHEMA),
                                                   comm)
        else:
            parser.error("A YAML-format file specifying run parameters must be specified.")

        # NOTE- outputNotesName is defined only for rank 0.
        if opts.out:
            outputNotesName = opts.out
        elif 'notesFileName' in CL_DATA['input']:
            outputNotesName = CL_DATA['input']['notesFileName']
        else:
            outputNotesName = DEFAULT_OUTPUT_NOTES_NAME
        CL_DATA['outputNotesName'] = outputNotesName  # make this available to the other ranks anyways
        parser.destroy()
    else:
        CL_DATA = None

    try:
        patchGroup = None  # for the benefit of diagnostics during init
        CL_DATA = comm.bcast(CL_DATA, root=0)

        pyrheautils.prepPathTranslations(CL_DATA['input'])

        # make output notes name available to the constants file for magic reasons
        pyrheautils.outputNotesName = CL_DATA['outputNotesName']
        pyrheautils.saveNewConstants = CL_DATA['saveNewConstants']

        # loading up the constants replacements needs to happen fairly early
        if CL_DATA['constantsFile'] is not None:
            pyrheautils.readConstantsReplacementFile(CL_DATA['constantsFile'])

        configureLogging(CL_DATA['logCfgDict'], CL_DATA['loggingExtra'])

        verbose = CL_DATA['verbose']  # @UnusedVariable
        debug = CL_DATA['debug']  # @UnusedVariable
        trace = CL_DATA['trace']
        deterministic = CL_DATA['deterministic']
        inputDict = CL_DATA['input']

        if deterministic:
            seed(1234 + comm.rank)  # Set the random number generator seed
        elif CL_DATA['randomSeed']:
            seed(CL_DATA['randomSeed'] + comm.rank)
        elif 'randomSeed' in inputDict:
            seed(inputDict['randomSeed'] + comm.rank)

        if 'trackedFacilities' in inputDict:
            global _TRACKED_FACILITIES
            if _TRACKED_FACILITIES:
                _TRACKED_FACILITIES.extend(inputDict['trackedFacilities'])
            else:
                _TRACKED_FACILITIES = inputDict['trackedFacilities'][:]

        schemautils.setSchemaBasePath(SCHEMA_DIR)

        pthImplDict = loadPathogenImplementations(pyrheautils.pathTranslate(PTH_IMPLEMENTATIONS_DIR))
        assert len(pthImplDict) == 1, 'Simulation currently supports exactly one pathogen'
        PthClass = pthImplDict.values()[0].getPathogenClass()
        pthName = pthImplDict.values()[0].pathogenName
    
        facImplDict = loadFacilityImplementations(inputDict['facilityImplementationDir'])
        if 'facilitySelectors' in inputDict:
            facImplRules = [(re.compile(rule['category']), rule['implementation'])
                            for rule in inputDict['facilitySelectors']]
        else:
            facImplRules = [(re.compile(category), category)
                            for category in facImplDict.keys()]  # an identity map

        policyClassList = loadPolicyImplementations(inputDict['policyImplementationDir'])
        policyRules = [(re.compile(rule['category']), re.compile(rule['policyClass']),
                        rule['category'], rule['policyClass'])
                       for rule in inputDict['policySelectors'] if 'category' in rule.keys() ]
        policyRules += [(re.compile(rule['locationAbbrev']),re.compile(rule['policyClass']),
                         rule['locationAbbrev'], rule['policyClass'])
                           for rule in inputDict['policySelectors']
                           if 'locationAbbrev' in rule.keys() ]
        # Add a data structure to track whether rules get used
        policyRulesDict = {pR: False for pR in policyRules}

        myFacList = distributeFacilities(comm, inputDict['facilityDirs'],
                                         facImplDict, facImplRules,
                                         CL_DATA['partitionFile'])
        LOGGER.info('Rank %d has facilities %s' % (comm.rank, [rec['abbrev'] for rec in myFacList]))

        patchGroup = patches.PatchGroup(comm, trace=trace, deterministic=deterministic,
                                        printCensus=CL_DATA['printCensus'])
        noteHolderGroup = noteholder.NoteHolderGroup()

        if comm.rank == 0 and CL_DATA['checkpoint'] != -1:
            checkpointer = checkpoint.checkpoint(CL_DATA['checkpoint'])
            patchList = [patchGroup.addPatch(patches.Patch(patchGroup, checkpointer=checkpointer))]
            for i in xrange(CL_DATA['patches'] - 1):
                patchList.append(patchGroup.addPatch(patches.Patch(patchGroup)))
        else:
            patchList = [patchGroup.addPatch(patches.Patch(patchGroup))
                         for i in xrange(CL_DATA['patches'])]  # @UnusedVariable
        totalRunDays = inputDict['runDurationDays'] + inputDict['burnInDays']

        # Add some things for which we need only one instance
        patchList[0].addAgents([BurnInAgent('burnInAgent', patchList[0],
                                            inputDict['burnInDays'], noteHolderGroup)])
        if 'scenarioWaitDays' in inputDict:
            scenarioPolicyClasses = []
            policyClasses = (findPolicies(policyClassList, policyRulesDict, 'scenario')
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

        # Check that constants replacements happened as expected
        pyrheautils.checkReplacementsWereUsed()
        if CL_DATA['saveNewConstants']:
            logging.shutdown()
            sys.exit('shutting down after saving new constants')

        initializeFacilities(patchList, myFacList, facImplDict, facImplRules,
                             policyClassList, policyRulesDict,
                             PthClass, noteHolderGroup, comm, totalRunDays)

        if CL_DATA['disableNotes']:
            noteHolderGroup.disableAll()
            
        monitorList = []
        tauAdjusterList = []
        if CL_DATA['bczmonitor'] is not None:
            for patch in patchList:
                m = bcz_monitor.Monitor(CL_DATA['bczmonitor'], patch)
                monitorList.append(m)
                patch.loop.addPerDayCallback(m.createDailyCallbackFn())
                if 'trackedValues' in inputDict:
                    for tv in inputDict['trackedValues']:
                        m.registerTrackable(tv)
                m.solidifyFields()

                if CL_DATA['taumod']:
                    ta = TauAdjuster(m)
                    tauAdjusterList.append(ta)
                    m.setStopTimeFn(1, ta.createCallbackFn())

        if CL_DATA['dumpFacilitiesMap'] is not None:
            dumpFacilitiesMap(CL_DATA['dumpFacilitiesMap'], patchList)

        # Check that all policy rules have been used, to avoid a common user typo problem
        quitNow = False
        for ruleKey, usedFlag in policyRulesDict.items():
            if not usedFlag:
                LOGGER.error('The policy rule associating %s with %s was never used- typo?',
                             ruleKey[2], ruleKey[3])
                quitNow = True
        if quitNow:
            raise RuntimeError('Probable typos found in the policy section of the input file')


    except Exception as e:
        if patchGroup:
            if LOGGER is None:
                LOGGER = logging.getLogger(__name__)
            LOGGER.error('%s exception during initialization: %s; traceback follows' %
                         (patchGroup.name, e))
            import traceback
            traceback.print_exc(file=sys.stderr)
            LOGGER.info('%s exiting due to exception during initialization' % patchGroup.name)
            logging.shutdown()
            sys.exit('Exception during initialization')
        else:
            if LOGGER is None:
                LOGGER = logging.getLogger(__name__)
            LOGGER.error('%s exception during initialization: %s; traceback follows' %
                         ('<no patchGroup yet>', e))
            import traceback
            traceback.print_exc(file=sys.stderr)
            LOGGER.info('%s exiting due to exception during initialization' % '<no patchGroup yet>')
            logging.shutdown()
            sys.exit('Exception during initialization')

    LOGGER.info('%s #### Ready to Run #### (from main)' % patchGroup.name)
    runFailed = False
    try:
        exitMsg = patchGroup.start()
        LOGGER.info('%s #### all done #### (from main); %s' % (patchGroup.name, exitMsg))
    except Exception as e:
        LOGGER.error('%s exception %s; traceback follows' % (patchGroup.name, e))
        import traceback
        traceback.print_exc(file=sys.stderr)
        runFailed = True
    finally:
        try:
            LOGGER.info('%s writing notes and exiting' % patchGroup.name)
    
            allNotesGroup, allNotesList = collectNotes(noteHolderGroup, comm)  # @UnusedVariable
            if comm.rank == 0:
                d = {}
                for nh in allNotesGroup.getnotes():
                    d[nh['name']] = nh.getDict()
        #             if 'occupancy' in nh:
        #                 recs = nh['occupancy']
        #                 with open(('occupancy_%s.csv' % nh['name']), 'w') as f:
        #                     csv_tools.writeCSV(f, recs[1].keys(), recs)
                if runFailed:
                    # edit output file name to show failure
                    baseNm, extNm = os.path.splitext(outputNotesName)
                    outputNotesName = baseNm + '_FAILED' + extNm
                if outputNotesName.lower().endswith('.json'):
                    with open(outputNotesName, 'w') as f:
                        f.write('{\n')
                        for k,v in d.items():
                            f.write('"{0}":{1},\n'.format(k,ujson.dumps(v)))
                        f.write("}\n")
                else:
                    with open(outputNotesName, 'w') as f:
                        pickle.dump(d, f)

                if CL_DATA['bczmonitor'] is not None:
                    for m in monitorList:
                        m.writeData()

        except Exception as e:
            LOGGER.error('%s an exception occurred while writing notes: %s'
                         % (patchGroup.name, e))

    logging.shutdown()

############
# Main hook
############

if __name__ == "__main__":
    main()
