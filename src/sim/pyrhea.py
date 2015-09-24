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
import jsonschema
from random import seed
from collections import defaultdict
import optparse
import logging
import logging.config

import patches
import yaml_tools

inputSchema = 'rhea_input_schema.yaml'

patchGroup = None
logger = None


def loadFacilityImplementations(implementationDir):
    logger.info('Loading facility implementations')
    implDict = {}
    sys.path.append(implementationDir)
    for fname in os.listdir(implementationDir):
        name, ext = os.path.splitext(fname)
        if ext == '.py':
            newMod = load_source(name, os.path.join(implementationDir, fname))
            for requiredAttr in ['category', 'generateFull']:
                if not hasattr(newMod, requiredAttr):
                    raise RuntimeError('The facility implemenation in %s has no %s' %
                                       (os.path.join(implementationDir, fname),
                                        requiredAttr))
            assert newMod.category not in implDict, \
                "Redundant definitions for category %s" % newMod.category
            implDict[newMod.category] = newMod
            logger.info('Loaded %s implementing %s' % (fname, newMod.category))
        elif ext.startswith('.py'):
            pass  # .pyc files for example
        else:
            logger.info('Skipping non-python file %s' % fname)
    sys.path.pop()  # drop implementaionDir
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


def localCollectResults():
    d = {}
    for cl in patches.Agent.__subclasses__():
        if hasattr(cl, 'generateReport'):
            d.update(cl.generateReport())
    for cl in patches.Interactant.__subclasses__():
        if hasattr(cl, 'generateReport'):
            d.update(cl.generateReport())

    return d


def collectResults(comm):
    if comm.rank == 0:
        resultDict = localCollectResults()
        for targetRank in xrange(comm.size):
            if targetRank != comm.rank:
                d = comm.recv(source=targetRank)
                for k, v in d.items():
                    if k in resultDict:
                        resultDict[k] += v
                    else:
                        resultDict[k] = v
        return resultDict
    else:
        d = localCollectResults()
        comm.send(d, dest=0)
        return None


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
    try:
        with open(fname, 'rU') as f:
            inputJSON = yaml.safe_load(f)
        with open(os.path.join(os.path.dirname(__file__), schemaFname), 'rU') as f:
            schemaJSON = yaml.safe_load(f)
        validator = jsonschema.validators.validator_for(schemaJSON)(schema=schemaJSON)
        nErrors = sum([1 for e in validator.iter_errors(inputJSON)])  # @UnusedVariable
        if nErrors:
            logger.error('Input file violates schema:')
            for e in validator.iter_errors(inputJSON):
                logger.error('Schema violation: %s: %s' %
                             (' '.join([str(word) for word in e.path]), e.message))
            comm.Abort(2)
        else:
            return inputJSON
    except Exception, e:
        logger.error('Error checking input against its schema: %s' % e)
        comm.Abort(2)


def createPerDayCB(patchGroup, runDurationDays):
    def perDayCB(loop, timeNow):
        if timeNow > runDurationDays:
            patchGroup.stop()
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

import pika
import pickle

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


def main():
    global patchGroup

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
            clData['input'] = checkInputFileSchema(args[0], inputSchema, comm)
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

    facImplDict = loadFacilityImplementations(inputDict['facilityImplementationDir'])

    myFacList = distributeFacilities(comm, inputDict['facilityDirs'], facImplDict)
    logger.info('Rank %d has facilities %s' % (comm.rank, [rec['abbrev'] for rec in myFacList]))

    patchGroup = patches.PatchGroup(comm, trace=trace, deterministic=deterministic,
                                    printCensus=clData['printCensus'])
    patchList = [patchGroup.addPatch(patches.Patch(patchGroup))
                 for i in xrange(clData['patches'])]  # @UnusedVariable
    for patch in patchList:
        patch.loop.addPerDayCallback(createPerDayCB(patchGroup, inputDict['runDurationDays']))
    tupleList = [(patch, [], []) for patch in patchList]

    offset = 0
    for facDescription in myFacList:
        patch, allIter, allAgents = tupleList[offset]
        if facDescription['category'] in facImplDict:
            facilities, wards, patients = \
                facImplDict[facDescription['category']].generateFull(facDescription, patch)
            allIter.extend([fac.reqQueue for fac in facilities])
            allIter.extend([fac.holdQueue for fac in facilities])
            allIter.extend(wards)
            allAgents.extend([fac.manager for fac in facilities])
            allAgents.extend(patients)
        else:
            raise RuntimeError('Facility %(abbrev)s category %(category)s has no implementation' %
                               facDescription)
        offset = (offset + 1) % len(tupleList)

    for patch, allIter, allAgents in tupleList:
        patch.addInteractants(allIter)
        patch.addAgents(allAgents)
        logger.info('Rank %d patch %s: %d interactants, %d agents' % (comm.rank, patch.name,
                                                                      len(allIter),
                                                                      len(allAgents)))

    exitMsg = patchGroup.start()
    logger.info('%s all done (from main); %s' % (patchGroup.name, exitMsg))

    resultDict = collectResults(comm)
    if comm.rank == 0:
        print '-' * 40
        print 'Collected Results:'
        for k, v in resultDict.items():
            print '%s: %s' % (k, v)

    logging.shutdown()

############
# Main hook
############

if __name__ == "__main__":
    main()
