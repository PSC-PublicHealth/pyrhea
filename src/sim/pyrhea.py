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

_rhea_svn_id_ = "$Id$"

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
import pickle
import re

import patches
import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.notes.noteholder as noteholder
import schemautils
import pyrheautils

schemaDir = '../schemata'
inputSchema = 'rhea_input_schema.yaml'

logger = None


def loadPolicyImplementations(implementationDir):
    logger.info('Loading policy implementations')
    implList = []
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
    implDict = {}
    for newMod in pyrheautils.loadModulesFromDir(implementationDir,
                                                 requiredAttrList=['category',
                                                                   'generateFull']):
        assert newMod.category not in implDict, \
            "Redundant definitions for category %s" % newMod.category
        implDict[newMod.category] = newMod
        logger.info('Loaded %s implementing %s' % (newMod.__name__, newMod.category))
    return implDict


def loadFacilityDescriptions(dirList, facImplDict):
    """
    This loads facilities, including communities, checking them against a schema.
    """
    logger.info('Parsing facilities descriptions')
    allKeySet = set()
    rawRecs = []
    for d in dirList:
        kS, recs = yaml_tools.parse_all(d)
        rawRecs.extend(recs)
        allKeySet.update(kS)
    facRecs = []
    for rec in rawRecs:
        assert 'abbrev' in rec, "Facility description has no 'abbrev' field: %s" % rec
        assert 'category' in rec, \
            "Facility description for %(abbrev) has no 'category' field" % rec
        nErrors = facImplDict[rec['category']].checkSchema(rec)
        if nErrors:
            logger.warning('dropping %s; %d schema violations' % (rec['abbrev'], nErrors))
        else:
            facRecs.append(rec)
    return facRecs


def distributeFacilities(comm, facDirs, facImplDict):
    """
    Estimate work associated with each facility and attempt to share it out fairly
    """
    if comm.rank == 0:
        facRecs = loadFacilityDescriptions(facDirs, facImplDict)
        assignments = defaultdict(list)
        tots = {k: 0 for k in xrange(comm.size)}
        pairL = [(facImplDict[rec['category']].estimateWork(rec), rec) for rec in facRecs]
        pairL.sort(reverse=True)
        for newWork, rec in pairL:
            leastWork, leastBusy = min([(v, k) for k, v in tots.items()])  # @UnusedVariable
            assignments[leastBusy].append(rec)
            tots[leastBusy] += newWork
        logger.info('Estimated work by rank: %s' % tots)
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


def checkInputFileSchema(fname, schemaFname, comm):
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
            comm.Abort(2)
        else:
            return inputJSON
    except Exception, e:
        myLogger.error('Error checking input against its schema: %s' % e)
        comm.Abort(2)


def buildFacOccupancyDict(patch, timeNow):
    facTypeDict = {'day': timeNow}
    assert hasattr(patch, 'allFacilities'), 'patch %s has no list of facilities!' % patch.name
    for fac in patch.allFacilities:
        tpName = type(fac).__name__
        if tpName not in facTypeDict:
            facTypeDict[tpName] = 0
            facTypeDict[tpName + '_all'] = 0
        if hasattr(fac, 'patientDataDict'):
            patientCount = 0
            allCount = 0
            for rec in fac.patientDataDict.values():
                allCount += 1 + rec.prevVisits
                if rec.departureDate is None:
                    patientCount += 1
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


def createPerDayCB(patch, noteHolder, runDurationDays):
    def perDayCB(loop, timeNow):
#         for p in patchGroup.patches:
#             p.loop.printCensus()
        oD = buildFacOccupancyDict(patch, timeNow)
        pD = buildFacPthDict(patch, timeNow)
        noteHolder.addNote({'occupancy': [oD], 'pathogen': [pD]})
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


def findPolicies(policyClassList,
                 policyRules,
                 category):
    l = []
    for pCl in policyClassList:
        for categoryRegex, classRegex in policyRules:
            if categoryRegex.match(category) and classRegex.match(pCl.__name__):
                l.append(pCl)
    return l

def main():

    comm = patches.getCommWorld()
    logging.basicConfig(format="[%d]%%(levelname)s:%%(name)s:%%(message)s" % comm.rank)

    if comm.rank == 0:
        parser = TweakedOptParser(usage="""
        %prog [-v][-d][-t][-D][-p=npatch][-C][-L=loglevel] input.yaml
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

        opts, args = parser.parse_args()
        if opts.log is not None:
            numLogLevel = getattr(logging, opts.log.upper(), None)
            if not isinstance(numLogLevel, int):
                parser.error("Invalid log level: %s" % opts.log)
        else:
            numLogLevel = None
        clData = {'verbose': opts.verbose,
                  'debug': opts.debug,
                  'trace': opts.trace,
                  'deterministic': opts.deterministic,
                  'patches': opts.patches,
                  'printCensus': opts.census,
                  'logCfgDict': getLoggerConfig(),
                  'loggingExtra': numLogLevel
                  }
        if len(args) == 1:
            clData['input'] = checkInputFileSchema(args[0],
                                                   os.path.join(schemaDir, inputSchema),
                                                   comm)
        else:
            parser.error("A YAML-format file specifying run parameters must be specified.")
        parser.destroy()
    else:
        clData = None
    clData = comm.bcast(clData, root=0)

    configureLogging(clData['logCfgDict'], clData['loggingExtra'])

    verbose = clData['verbose']  # @UnusedVariable
    debug = clData['debug']  # @UnusedVariable
    trace = clData['trace']
    deterministic = clData['deterministic']
    inputDict = clData['input']

    if deterministic:
        seed(1234 + comm.rank)  # Set the random number generator seed

    schemautils.setSchemaBasePath(schemaDir)
    pthImplDict = loadPathogenImplementations(inputDict['pathogenImplementationDir'])
    assert len(pthImplDict) == 1, 'Simulation currently supports exactly one pathogen'
    PthClass = pthImplDict.values()[0].getPathogenClass()

    facImplDict = loadFacilityImplementations(inputDict['facilityImplementationDir'])

    policyClassList = loadPolicyImplementations(inputDict['policyImplementationDir'])
    policyRules = [(re.compile(rule['category']), re.compile(rule['policyClass']))
                   for rule in inputDict['policySelectors']]

    myFacList = distributeFacilities(comm, inputDict['facilityDirs'], facImplDict)
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

    # Share the facilities out over the patches
    offset = 0
    tupleList = [(p, [], [], []) for p in patchList]
    for facDescription in myFacList:
        patch, allIter, allAgents, allFacilities = tupleList[offset]
        if facDescription['category'] in facImplDict:
            facImpl = facImplDict[facDescription['category']]
            facilities, wards, patients = \
                facImpl.generateFull(facDescription, patch,
                                     policyClasses=findPolicies(policyClassList,
                                                                policyRules,
                                                                facImpl.category))
            for w in wards:
                w.setInfectiousAgent(PthClass(w))
                w.initializePatientPthState()
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
        patch.loop.addPerDayCallback(createPerDayCB(patch, patchNH, totalRunDays))
        patchNH.addNote({'name': patch.name,
                         'occupancy': [buildFacOccupancyDict(patch, 0)],
                         'rank': comm.rank})

    logger.info('%s #### Ready to Run #### (from main)' % patchGroup.name)
    exitMsg = patchGroup.start()
    logger.info('%s #### all done #### (from main); %s' % (patchGroup.name, exitMsg))

    allNotesGroup, allNotesList = collectNotes(noteHolderGroup, comm)  # @UnusedVariable
    if comm.rank == 0:
        d = {}
        for nh in allNotesGroup.getnotes():
            d[nh['name']] = nh.getDict()
#             if 'occupancy' in nh:
#                 recs = nh['occupancy']
#                 with open(('occupancy_%s.csv' % nh['name']), 'w') as f:
#                     csv_tools.writeCSV(f, recs[1].keys(), recs)
        with open('notes.pkl', 'w') as f:
            pickle.dump(d, f)

    logging.shutdown()

############
# Main hook
############

if __name__ == "__main__":
    main()
