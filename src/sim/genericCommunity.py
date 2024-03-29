#! /usr/bin/env python

###################################################################################
# Copyright 2015-16, Pittsburgh Supercomputing Center (PSC). All Rights Reserved. #
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

import os.path
import sys
import random
import logging
import types
from collections import defaultdict

from phacsl.utils.collections.phacollections import DefaultDict
import phacsl.utils.collections.interdict as interdict
import cPickle as pickle
import gzip
from scipy.stats import expon, binom

import pyrheabase
import pyrheautils
from typebase import DiagClassA, CareTier, PatientOverallHealth
from freezerbase import FreezerError
from facilitybase import TreatmentProtocol, BirthQueue, HOMEQueue  # @UnusedImport
from facilitybase import Facility, Ward, PatientAgent, PatientStatusSetter, PatientRecord
from facilitybase import ClassASetter, PatientStatus, PatientDiagnosis, FacilityManager
from facilitybase import MissingPatientRecordError, decodeHistoryEntry, findQueueForTier
from quilt.netinterface import GblAddr
from stats import CachedCDFGenerator, BayesTree
import schemautils
from pathogenbase import PthStatus

import time
import psutil

logger = logging.getLogger(__name__)
infectionLogger = logging.getLogger("infectionTracking")

category = 'COMMUNITY'
_schema = 'communityfacts_schema.yaml'
_validator = None
_constants_values = '$(MODELDIR)/constants/community_constants.yaml'
_constants_schema = 'community_constants_schema.yaml'
_constants = None

BYPASS_KEY = '_bypass_'  # used by some derived classes

USE_CUSTOM_ENCODING = False  # the custom encoding turns out to be only marginally beneficial

def importCommunity(moduleName):
    attrs = (
        'category',
        '_schema',
        '_validator',
        '_constants_values',
        '_constants_schema',
        'CommunityWard',
        'CommunityManager',
        'Community',
        'generateFull',
        'estimateWork',
        'checkSchema')

    global logger
    logger = logging.getLogger(moduleName)

    module = sys.modules[moduleName]

    curModule = sys.modules[__name__]
    for attr in attrs:
        if not hasattr(module, attr):
            setattr(module, attr, getattr(curModule, attr))
    ###########
    # Load constants
    ###########
    global _constants
    _constants = pyrheautils.importConstants(_constants_values,
                                             _constants_schema)
    if not hasattr(module, '_constants'):
        setattr(module, '_constants', _constants)



def lencode(thing):
    if isinstance(thing, types.DictType):
        resultL = []
        tpL = []
        valL = []
        kL = thing.keys()
        kL.sort()
        for k in kL:
            subResultL, subTpL, subValL = lencode(thing[k])
            resultL.append((k, subResultL))
            tpL.extend(subTpL)
            valL.extend(subValL)
        return ('dict', resultL), tpL, valL
    elif isinstance(thing, (types.TupleType, types.ListType)):
        resultL = []
        tpL = []
        valL = []
        for v in list(thing):
            subResultL, subTpL, subValL = lencode(v)
            resultL.append(subResultL)
            tpL.extend(subTpL)
            valL.extend(subValL)
        return (type(thing).__name__, resultL), tpL, valL
    else:
        return type(thing).__name__, [type(thing).__name__], [thing]


def ldecode(typeTpl, valL):
    if isinstance(typeTpl, types.TupleType):
        k, info = typeTpl
    else:
        k = typeTpl
        info = []
    if k == 'dict':
        d = {}
        for tpl in info:
            k, elt = tpl
            newVal, valL = ldecode(elt, valL)
            d[k] = newVal
        return d, valL
    elif k in ['tuple', 'list', 'PatientStatus', 'PatientDiagnosis', 'GblAddr']:
        outL = []
        for elt in info:
            newVal, valL = ldecode(elt, valL)
            outL.append(newVal)
        if k == 'tuple':
            result = tuple(outL)
        elif k == 'list':
            result = outL
        elif k == 'PatientStatus':
            result = PatientStatus(*outL)
        elif k == 'PatientDiagnosis':
            result = PatientDiagnosis(*outL)
        elif k == 'GblAddr':
            result = GblAddr(*outL)
        return result, valL
    elif k == 'list':
        outL = []
        for elt in info:
            newVal, valL = ldecode(elt, valL)
            outL.append(newVal)
        return outL, valL
    else:
        return valL[0], valL[1:]

cacheVer = 8
LastMemCheck = time.time()

class Freezer(object):
    def __init__(self, ward):
        self.ward = ward
        self.frozenAgentClass = PatientAgent
        self.frozenAgentLoggerName = None
        self.frozenAgentTypePattern = None
        self.frozenAgentList = []

    def freezeAndStore(self, agent):
        if agent in self.ward.lockingAgentSet:
            self.ward.lockingAgentSet.remove(agent)
        self.ward.suspend(agent)
        assert agent not in self.ward.lockingAgentSet, 'It is still there!'
        assert agent in self.ward._lockQueue, 'It is not there!'
        self.ward._lockQueue.remove(agent)
        d = agent.__getstate__()
        agent.kill()
        if d['newLocAddr'] is None or d['newLocAddr'] == self.ward.getGblAddr():
            del d['newLocAddr']
        else:
            raise FreezerError('%s cannot freezedry %s because it is still in motion'
                               % (self.ward._name, d['name']))
        if self.frozenAgentLoggerName is None:
            self.frozenAgentLoggerName = d['loggerName']
            del d['loggerName']
        elif d['loggerName'] == self.frozenAgentLoggerName:
            del d['loggerName']
        else:
            raise FreezerError('%s cannot freezedry %s because it has the wrong logger'
                               % (self.ward._name, d['name']))
        if USE_CUSTOM_ENCODING:
            typeTpl, linL, valL = lencode(d)  # @UnusedVariable
            #print 'dictionary: %s' % str(d)
            #print 'types: %s' % str(typeTpl)
            #print 'linearized: %s' % str(linL)
            #print 'values: %s' % str(valL)
        else:
            typeTpl = 1
            valL = pickle.dumps(d, 2)
        if self.frozenAgentTypePattern is None:
            self.frozenAgentTypePattern = typeTpl
        elif typeTpl != self.frozenAgentTypePattern:
            raise FreezerError('%s cannot freezedry %s because it has the wrong type pattern'
                               % (self.ward._name, d['name']))
        if USE_CUSTOM_ENCODING:
            self.frozenAgentList.append(tuple(valL))
        else:
            self.frozenAgentList.append(valL)

    def removeAndThaw(self, frozenAgent, timeNow):
        self.frozenAgentList.remove(frozenAgent)

        if USE_CUSTOM_ENCODING:
            valL = list(frozenAgent)
            d, leftovers = ldecode(self.frozenAgentTypePattern, valL)
            assert not leftovers, ('%s had %s left over unfreezing agent'
                                   % (self.ward._name, leftovers))
        else:
            d = pickle.loads(frozenAgent)
        d['locAddr'] = self.ward.getGblAddr()
        d['newLocAddr'] = self.ward.getGblAddr()
        d['loggerName'] = self.frozenAgentLoggerName
        agent = self.frozenAgentClass.__new__(self.frozenAgentClass)
        agent.__setstate__(d)
        agent.reHome(self.ward.patch)
        # NEED to check locking agent count
        self.ward._lockQueue.append(agent)
        self.ward.awaken(agent)
        self.ward.lockingAgentSet.add(agent)
        return agent


class InvalidCacheVer(Exception):
    pass

def mapLMDBAbbrev(abbrev):
    return 'all'

def mapLMDBSegments(abbrev, iDict):
    global cacheVer
    segmentInfoId = 0
    idsPerFac = 1000000000
    startId = 100

    try:
        segmentInfo = iDict[segmentInfoId]
    except:
        segmentInfo = [cacheVer, startId, {}]
        # this will necessarily get written later so don't bother yet.

    if segmentInfo[0] != cacheVer:
        raise InvalidCacheVer()

    cv, nextId, segmentDict = segmentInfo
    if abbrev in segmentDict:
        return segmentDict[abbrev]

    ret = segmentDict[abbrev] = nextId
    nextId += idsPerFac
    iDict[segmentInfoId] = [cv, nextId, segmentDict]
    return ret


class CopyOnWriteLMDBFreezer(Freezer):
    _AgentListId = 0
    _SavedInfoListId = 1
    _SavedWardDataId = 2

    def __init__(self, ward, orig, changed, infoList, pType):
        """
        This is a freezer for community agents that is stored in an LMDB InterDict rather than
        in memory. It is stored in two parts so that the community agents can be cached and then
        changed community members are written to a second LMDB file.
        
        A given ward can have multiple freezers to store multiple types of community agents but
        it's useful to have all freezers in one pair of interdicts.  To this end there is an
        infoPointer which is a list of values that help the multiple freezers interact.

        ward: this is the ward that is associated with the freezer
        orig: this is the InterDict where the cached community agents are stored.  By default this
              is not modified.
        changed: any newly stored agents are saved here
        infoList: this is a list shared between all freezers using the same pair of interdicts with
                  the following values:
            origRO: if True (normal usage) write new agents to changed.  For building the cache set
                    to False and all new members will be written to orig
            maxOrig: when looking up an agent, any ID less than this will be found in orig, otherwise
                     found in changed
            nextID: next ID value to write
        pType: this is the key used for the dict of freezers in the ward
        """
        self.ward = ward
        self.frozenAgentClass = PatientAgent
        self.frozenAgentLoggerName = None
        self.frozenAgentTypePattern = None
        self.orig = orig
        self.changed = changed
        self.infoList = infoList
        self.pType = pType
        self.iDictOffset = self.ward.iDictOffset
        try:
            self.frozenAgentList = self.orig[self._AgentListId+self.iDictOffset][pType]
        except:
            self.frozenAgentList = set()

    def freezeAndStore(self, agent):
        if agent in self.ward.lockingAgentSet:
            self.ward.lockingAgentSet.remove(agent)
        self.ward.suspend(agent)
        assert agent not in self.ward.lockingAgentSet, 'It is still there!'
        assert agent in self.ward._lockQueue, 'It is not there!'
        self.ward._lockQueue.remove(agent)
        d = agent.__getstate__()
        agent.kill()
        if d['newLocAddr'] is None or d['newLocAddr'] == self.ward.getGblAddr():
            del d['newLocAddr']
        else:
            raise FreezerError('%s cannot freezedry %s because it is still in motion'
                               % (self.ward._name, d['name']))
        if self.frozenAgentLoggerName is None:
            self.frozenAgentLoggerName = d['loggerName']
            del d['loggerName']
        elif d['loggerName'] == self.frozenAgentLoggerName:
            del d['loggerName']
        else:
            raise FreezerError('%s cannot freezedry %s because it has the wrong logger'
                               % (self.ward._name, d['name']))

        origRO, maxOrig, nextId = self.infoList
        if origRO:
            self.changed[nextId] = d
        else:
            self.orig[nextId] = d
            self.infoList[1] = nextId
        self.frozenAgentList.add(nextId)
        self.infoList[2] += 1


    def removeAndThaw(self, frozenAgent, timeNow):
        self.frozenAgentList.remove(frozenAgent)
        origRO, maxOrig, nextId = self.infoList
        if frozenAgent <= maxOrig:
            d = self.orig[frozenAgent]
        else:
            d = self.changed[frozenAgent]
            del self.changed[frozenAgent]

        d['locAddr'] = self.ward.getGblAddr()
        d['newLocAddr'] = self.ward.getGblAddr()
        d['loggerName'] = self.frozenAgentLoggerName
        agent = self.frozenAgentClass.__new__(self.frozenAgentClass)
        agent.__setstate__(d)
        agent.reHome(self.ward.patch)
        # NEED to check locking agent count
        self.ward._lockQueue.append(agent)
        self.ward.awaken(agent)
        self.ward.lockingAgentSet.add(agent)
        global LastMemCheck
        if time.time() > 60.0 + LastMemCheck:
            p = psutil.Process()
            print p.memory_full_info()
            LastMemCheck = time.time()
        return agent

    def saveFreezerData(self):
        origRO, maxOrig, nextId = self.infoList
        if origRO:
            lDict = self.changed
        else:
            lDict = self.orig
        try:
            agentListDict = lDict[self._AgentListId+self.iDictOffset]
        except:
            agentListDict = {}
        agentListDict[self.pType] = self.frozenAgentList
        lDict[self._AgentListId+self.iDictOffset] = agentListDict


    def saveInfoList(self):
        # save the infoList
        # however we should save it with origRO flag set to True
        origRO, maxOrig, nextId = self.infoList
        if origRO:
            lDict = self.changed
        else:
            lDict = self.orig
        tmpIL = self.infoList[:]
        tmpIL[0] = True
        lDict[self._SavedInfoListId+self.iDictOffset] = tmpIL

    def getInfoList(self, fromOrig=True):
        if fromOrig:
            return self.orig[self._SavedInfoListId+self.iDictOffset]
        else:
            return self.changed[self._SavedInfoListId+self.iDictOffset]
            
    def saveWardData(self, data):
        origRO, maxOrig, nextId = self.infoList
        if origRO:
            lDict = self.changed
        else:
            lDict = self.orig

        lDict[self._SavedWardDataId+self.iDictOffset] = data

    def getWardData(self, fromOrig=True):
        if fromOrig:
            return self.orig[self._SavedWardDataId+self.iDictOffset]
        else:
            return self.changed[self._SavedWardDataId+self.iDictOffset]



def makeLMDBDirs():
    origDir = pyrheautils.pathTranslate('$(COMMUNITYCACHEDIR)')
    if not os.path.exists(origDir):
        os.makedirs(origDir)

    changedDir = pyrheautils.pathTranslate('$(AGENTDIR)')
    if not os.path.exists(changedDir):
        os.makedirs(changedDir)

def origFname(abbrev):
    abbrev = mapLMDBAbbrev(abbrev)
    return pyrheautils.pathTranslate('$(COMMUNITYCACHEDIR)/freezer_%s'%abbrev)

def changedFname(abbrev):
    abbrev = mapLMDBAbbrev(abbrev)
    return pyrheautils.pathTranslate('$(AGENTDIR)/freezer_%s'%abbrev)

def patientDataFname(abbrev):
    abbrev = mapLMDBAbbrev(abbrev)
    return pyrheautils.pathTranslate('$(AGENTDIR)/patientData_%s'%abbrev)

def cachePatientDataFname(abbrev):
    abbrev = mapLMDBAbbrev(abbrev)
    return pyrheautils.pathTranslate('$(COMMUNITYCACHEDIR)/patientData_%s'%abbrev)


INTERDICT_MAPPING = {}


def newFreezer(dd, ward, orig, changed, infoList, key):
    freezer = CopyOnWriteLMDBFreezer(ward, orig, changed, infoList, key)
    dd[key] = freezer
    return freezer

class CommunityWard(Ward):
    """This 'ward' type represents being out in the community"""
    def __init__(self, abbrev, name, patch, nBeds):
        Ward.__init__(self, name, patch, CareTier.HOME, nBeds=nBeds)
        self.checkInterval = 1  # so they get freeze dried promptly

        oName = origFname(abbrev)
        if oName in INTERDICT_MAPPING:
            self.orig = INTERDICT_MAPPING[oName]
        else:
            self.orig = interdict.InterDict(oName, convert_int=True, val_serialization='pickle')
            INTERDICT_MAPPING[oName] = self.orig

        self.changed = {}

        try:
            self.iDictOffset = mapLMDBSegments(abbrev, self.orig)
        except:
            self.orig.close()
            self.orig = interdict.InterDict(oName, overwrite_existing=True,
                                            convert_int=True, val_serialization='pickle')
            self.iDictOffset = mapLMDBSegments(abbrev, self.orig)

        self.infoList = [True,1,1]  # we'll overwrite this shortly- deals with a chicken/egg problm
        self.freezers = DefaultDict(lambda dd,key: newFreezer(dd, self,
                                                              self.orig, self.changed,
                                                              self.infoList, key))
        self.newArrivals = []

    def classify(self, agent, timeNow):
        """Return the PatientCategory appropriate for this agent for freeze drying"""
        return "base"

    def flushNewArrivals(self):
        """
        Forget any new arrivals- this prevents them from being freezedried (possibly redundantly).
        """
        self.newArrivals = []

    def lock(self, lockingAgent):
        # print 'new lock %s next wake %s' % (lockingAgent.name, lockingAgent.nextWakeTime())
        self.newArrivals.append(lockingAgent)
        return super(CommunityWard, self).lock(lockingAgent)

    def getPatientList(self):
        """Because most of the patients are freeze-dried, we cannot supply a patient list"""
        unfrozenL = super(CommunityWard, self).getPatientList()
        if len(unfrozenL) == self.cumStats.popTotal():
            return unfrozenL
        else:
            raise FreezerError('Requested list cannot be generated without thawing everyone')

    def handlePatientArrival(self, patientAgent, timeNow):
        super(CommunityWard, self).handlePatientArrival(patientAgent, timeNow)
        # If a patient lands here, make this fac its home unless it's FRAIL
        # (and thus should not be landing here at all...)
        if timeNow is not None:
            infectionLogger.info("%s arriving in community at time %d with colonization status %s"%(
                patientAgent.name, timeNow, PthStatus.names[patientAgent.getStatus().pthStatus]))
        if patientAgent.getStatus().overall != PatientOverallHealth.FRAIL:
            patientAgent.setStatus(homeAddr=findQueueForTier(CareTier.HOME,
                                                             self.fac.reqQueues).getGblAddr())


class CommunityManager(FacilityManager):

    def getProbThaw(self, patientCategory, dT):
        return self.fac.cachedCDFs[patientCategory].intervalProb(0, dT)

    def perTickActions(self, timeNow):
        dT = (timeNow - self.fac.collectiveStatusStartDate if timeNow is not None else 0)
        countD = defaultdict(lambda: 0)
        for ward in self.fac.getWards():
            if dT != 0:
                for patCat, freezer in ward.freezers.items():
                    assert patCat in self.fac.cachedCDFs, ('%s has no CDF for patient category %s' %
                                                           (self.fac.name, patCat))
                    pThaw = self.getProbThaw(patCat, dT)
                    nFroz = len(freezer.frozenAgentList)
                    try:
                        nThawed = binom.rvs(nFroz, pThaw)
                    except ValueError:
                        logger.error('ValueError for %s: nFroz=%s, pThaw=%s',
                                     self.fac.name, nFroz, pThaw)
                        raise
                    self.fac.getNoteHolder().addNote({('thawed_%s' % timeNow) : nThawed })
                    changedList = random.sample(freezer.frozenAgentList, nThawed)
                    for a in changedList:
                        thawedAgent = freezer.removeAndThaw(a, timeNow)
                        if thawedAgent.debug:
                            thawedAgent.logger.debug('%s unfreezedried %s at %s'
                                                     % (ward._name, thawedAgent.name, timeNow))
                        countD[thawedAgent.getStatus().overall] += 1
                    ward.miscCounters['nThawed'] += nThawed

            for agent in ward.newArrivals:
                if agent.debug:
                    agent.logger.debug('%s freezedrying %s at %s'
                                       % (ward._name, agent.name, timeNow))
                ward.freezers[ward.classify(agent, timeNow)].freezeAndStore(agent)
            ward.newArrivals = []
        if dT != 0:
            logger.debug('%s (%s) at time %s thawed by category %s',
                         self.fac.name, ('%s_%s' % ward.getGblAddr().getLclAddr()), timeNow,
                         {PatientOverallHealth.names[k]: v for k, v in countD.items()})
            self.fac.collectiveStatusStartDate = timeNow

    def allocateAvailableBed(self, tier):
        """
        We are going to pretend that we never run out of beds, since the agents just get
        frozen on arrival.
        """
        if tier in self.fac.wardTierDict:
            availableWards = self.fac.wardTierDict[tier]
            w = availableWards[0]
            bestW, bestN = self.fac.wardAddrDict[w.getGblAddr().getLclAddr()]
            for w in availableWards[1:]:
                ww, nn = self.fac.wardAddrDict[w.getGblAddr().getLclAddr()]
                if nn > bestN:
                    bestW = ww
                    bestN = nn
            self.fac.wardAddrDict[bestW.getGblAddr().getLclAddr()] = (bestW, bestN-1)
            return bestW
        else:
            return None


class Community(Facility):
    def __init__(self, descr, patch, policyClasses=None, categoryNameMapper=None,
                 managerClass=None, wardClass=None):
        if managerClass is None:
            managerClass = CommunityManager
        if wardClass is None:
            wardClass = CommunityWard
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr,
                          descr, patch,
                          reqQueueClasses=[pyrheabase.FacRequestQueue, BirthQueue, HOMEQueue],
                          policyClasses=policyClasses,
                          managerClass=managerClass,
                          categoryNameMapper=categoryNameMapper)
        descr = self.mapDescrFields(descr)
        meanPop = descr['meanPop']['value']
        popScale = _constants['populationScale']['value']
        meanPop *= popScale
        nBeds = int(round(3.0*meanPop/popScale))
        self.setCDFs(_constants['losModelMap'])
        self.addWard(wardClass(descr['abbrev'], '%s_%s_%s' % (category, patch.name,
                                                              descr['abbrev']),
                               patch, nBeds))
        self.collectiveStatusStartDate = 0
        self.treeCache = {}
        self.trappedPatientFlowDct = {}

        pdName = cachePatientDataFname(descr['abbrev'])
        self.patientDataDict = {}
        if pdName in INTERDICT_MAPPING:
            self.cachePatientDataDict = INTERDICT_MAPPING[pdName]
        else:
            self.cachePatientDataDict = interdict.InterDict(pdName, overwrite_existing=False,
                                                            integer_keys=False,
                                                            key_serialization='msgpack',
                                                            val_serialization='pickle')
            global cacheVer
            if "cacheVer" not in self.cachePatientDataDict:
                self.cachePatientDataDict["cacheVer"] = cacheVer
            elif cacheVer != self.cachePatientDataDict["cacheVer"]:
                self.cachePatientDataDict.close()
                self.cachePatientDataDict = interdict.InterDict(pdName, overwrite_existing=True,
                                                                integer_keys=False,
                                                                key_serialization='msgpack',
                                                                val_serialization='pickle')
                self.cachePatientDataDict["cacheVer"] = cacheVer
            INTERDICT_MAPPING[pdName] = self.cachePatientDataDict
        self.patientCacheIsBeingRegenerated = False

    def finalizeBuild(self, facDescr):
        if self.patientCacheIsBeingRegenerated:

            ward = self.getWards()[0]  # There can be only one

            # Let's write all of our data at once to the interdict, so for now use a normal dict:
            orig = ward.orig
            ward.orig = {}
            # do the same with the patientDataDict
            realCachePatientDataDict = self.cachePatientDataDict
            self.cachePatientDataDict = {}  # make sure no one is found while writing these

            # the patient data dict is currently an empty dict where things will be initially
            # written.
            # just leave this and we'll move things around after everyone is generated.

            ward.infoList = [False, ward.iDictOffset+5, ward.iDictOffset+5]

            # Use perTickActions to force newly created patients to get freeze-dried
            self.manager.perTickActions(None)  # timeNow == None means before sim start

            freezerList = []
            for pType,freezer in ward.freezers.items():
                logger.debug('Ward %s has %d %s', ward._name, len(freezer.frozenAgentList), pType)
                freezerList.append((pType,freezer.frozenAgentLoggerName))
                freezer.saveFreezerData()

            # since it's possible that there's no freezers because there are no agents, let's make
            # another temp freezer
            tFreezer = CopyOnWriteLMDBFreezer(ward, orig, ward.changed, ward.infoList, False)
            tFreezer.saveWardData((cacheVer, facDescr, freezerList))
            tFreezer.saveInfoList()

            orig.mset(ward.orig.items())
            orig.flush()

            for k in ward.orig.keys():
                del ward.orig[k]

            ward.orig = orig
            for f in ward.freezers.values():
                f.orig = orig

            realCachePatientDataDict.mset(self.patientDataDict.items())
            self.cachePatientDataDict = realCachePatientDataDict
            for k in self.patientDataDict.keys():
                del self.patientDataDict[k]
            self.patientDataDict = {}

            self.patientCacheIsBeingRegenerated = False


    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        self.treeCache = {}

    def setCDFs(self, losModelMap):
        self.cachedCDFs = {}
        for classKey, losModel in losModelMap.items():
            assert losModel['pdf'] == 'expon(lambda=$0)', \
                "Unexpected losModel form %s - only expon is supported" % losModel['pdf']
            rate = losModel['parms'][0]
            #rate *= (779./675.)  # This is the one that scales the get-sick rate up to the go-home rate
            #rate *= (675./779.)
            self.cachedCDFs[classKey] = CachedCDFGenerator(expon(scale=1.0/rate))

    def calcTierRateConstants(self, flowKey):
        """
        These constants represent the fraction of the newly-sick that die, that
        need rehab as a (v)SNF, and that need care at the HOSP, ICU, LTAC, SKILNRS,
        and VENT tiers respectively.
        """
        deathRate = _constants['communityDeathRate']['value']
        needsICURate = _constants['communityVerySickRate']['value']
        needsLTACRate = _constants['communityNeedsLTACRate']['value']
        needsRehabRate = _constants['communityNeedsRehabRate']['value']
        needsSkilNrsRate = _constants['communityNeedsSkilNrsRate']['value']
        needsVentRate = (_constants['communityNeedsVentRate']['value']
                         if 'communityNeedsVentRate' in _constants else 0.0)
        needsHospRate = 1.0 - (deathRate + needsICURate + needsLTACRate + needsRehabRate
                               + needsSkilNrsRate + needsVentRate)
        return (deathRate, needsRehabRate, needsHospRate, needsICURate, needsLTACRate,
                needsSkilNrsRate, needsVentRate)

    def updateModifiers(self, patientAgent, modifierDct):
        for ent in reversed(patientAgent.agentHistory):
            if decodeHistoryEntry(ent)['category'] != 'COMMUNITY':
                prevFacAbbrev = decodeHistoryEntry(ent)['abbrev']
                break
        else:
            prevFacAbbrev = None
        modifierDct[pyrheabase.TierUpdateModKey.FLOW_KEY] = prevFacAbbrev

    def getStatusChangeTree(self, patientAgent, modifierDct, startTime, timeNow):  # @UnusedVariable
        patientStatus = patientAgent.getStatus()
        ward = patientAgent.ward
        careTier = ward.tier
        self.updateModifiers(patientAgent, modifierDct)
        flowKey = modifierDct[pyrheabase.TierUpdateModKey.FLOW_KEY]
        assert careTier == CareTier.HOME, \
            "The community only offers CareTier 'HOME'; found %s" % careTier
        if patientAgent.getDiagnosis().diagClassA == DiagClassA.WELL:
            self.trappedPatientFlowDct[patientAgent.id] = flowKey # in case patient gets trapped
            key = (flowKey, startTime - patientStatus.startDateA,
                   timeNow - patientStatus.startDateA)
            infectionLogger.info("%s leaving community at time %d with colonization status %s"%(
                patientAgent.name, timeNow, PthStatus.names[patientStatus.pthStatus]))

            if key in self.treeCache:
                return self.treeCache[key]
            else:
                changeProb = 1.0  # since the agent was already randomly chosen to be un-frozen
                try:
                    (deathRate, needsRehabRate, needsHospRate, needsICURate, needsLTACRate,
                     needsSkilNrsRate, needsVentRate) = self.calcTierRateConstants(flowKey)
                except:
                    logger.error('exception on historyL %s',
                                 str(reversed(patientAgent.agentHistory)))
                    raise
                tree = BayesTree(BayesTree.fromLinearCDF([(deathRate,
                                                           ClassASetter(DiagClassA.DEATH)),
                                                          (needsHospRate,
                                                           ClassASetter(DiagClassA.SICK)),
                                                          (needsICURate,
                                                           ClassASetter(DiagClassA.VERYSICK)),
                                                          (needsRehabRate,
                                                           ClassASetter(DiagClassA.NEEDSREHAB)),
                                                          (needsLTACRate,
                                                           ClassASetter(DiagClassA.NEEDSLTAC)),
                                                          (needsSkilNrsRate,
                                                           ClassASetter(DiagClassA.NEEDSSKILNRS)),
                                                          (needsVentRate,
                                                           ClassASetter(DiagClassA.NEEDSVENT)),
                                                          ], tag='FATE'),
                                 PatientStatusSetter(),
                                 changeProb, tag='LOS')
                self.treeCache[key] = tree
                return tree
        else:
            # This patient doesn't belong in this ward
            logger.warning('fac %s patient: %s careTier %s with status %s startTime: %s: '
                           'this patient should be gone by now',
                            self.name, patientAgent.name, CareTier.names[careTier],
                            DiagClassA.names[patientAgent.getStatus().diagClassA], startTime)

            # This patient must use the same flow key as when it first tried to depart, or the
            # outgoing flow may not provide the CareTier it needs.
            flowKey = self.trappedPatientFlowDct[patientAgent.id]
            modifierDct[pyrheabase.TierUpdateModKey.FLOW_KEY] = flowKey

            return BayesTree(PatientStatusSetter())

    def prescribe(self, ward, patientId, patientDiagnosis, patientTreatment, # @UnusedVariable
                  modifierDct, timeNow=None):  # @UnusedVariable
        """
        This returns a tuple (careTier, patientTreatment)

        modifierDct is a back channel for downstream communication.  It may be modified by
        routines for which it is a parameter.
        """
        if patientDiagnosis.diagClassA == DiagClassA.WELL:
            if (patientDiagnosis.overall == PatientOverallHealth.HEALTHY
                or patientDiagnosis.overall == PatientOverallHealth.UNHEALTHY):
                return (CareTier.HOME, patientTreatment._replace(rehab=False))
            elif patientDiagnosis.overall == PatientOverallHealth.FRAIL:
                return (CareTier.NURSING, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
            return (CareTier.NURSING, patientTreatment._replace(rehab=True))
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSLTAC:
            return (CareTier.LTAC, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.SICK:
            return (CareTier.HOSP, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.VERYSICK:
            return (CareTier.ICU, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSSKILNRS:
            return (CareTier.SKILNRS, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.DEATH:
            return (None, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSVENT:
            return (CareTier.VENT, patientTreatment._replace(rehab=False))
        else:
            raise RuntimeError('Unknown DiagClassA %s' % str(patientDiagnosis.diagClassA))

    def getPatientRecord(self, patientId, timeNow=None):
        k = (self.abbrev, patientId)
        # get your pickled patient record
        try:
            ppr = self.patientDataDict[k]
        except KeyError:
            try:
                ppr = self.cachePatientDataDict[k]
            except: # (interdict doesn't give keyerrors) KeyError:
                if timeNow is None:
                    raise MissingPatientRecordError('Record not found and cannot create a new one')
                pR = PatientRecord(patientId, timeNow, isFrail=False)
                self.patientDataDict[k] = pickle.dumps(pR, 2)  # keep a copy
                pR._owningFac = self
                return pR

        pR = pickle.loads(ppr)
        pR._owningFac = self
        return pR

    def forgetPatientRecord(self, patientId):
        """
        The facility forgets it ever saw this patient.  Used to implement record-keeping errors.
        """
        if patientId in self.patientDataDict or patientId in self.cachePatientDataDict:
            with self.getPatientRecord(patientId) as pRec:
                pRec.forgetPathogenInfo()

    def mergePatientRecord(self, patientId, newPatientRec, timeNow):
        patientRec = self.getPatientRecord(patientId, timeNow)
        patientRec.merge(newPatientRec)
        delattr(patientRec, '_owningFac')
        self.patientDataDict[(self.abbrev, patientId)] = pickle.dumps(patientRec, 2)

    def patientRecordExists(self, patientId):
        k = (self.abbrev, patientId)
        if k in self.patientDataDict:
            return True
        return k in self.cachePatientDataDict

    def getInitialOverallHealth(self, ward, timeNow):  # @UnusedVariable
        fracUnhealthy = _constants['initialUnhealthyFrac']['value']
        if random.random() <= fracUnhealthy:
            return PatientOverallHealth.UNHEALTHY
        else:
            return PatientOverallHealth.HEALTHY


def _populate(fac, descr, patch):
    assert 'meanPop' in descr, \
        "Community description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = float(descr['meanPop']['value'])
    popScale = _constants['populationScale']['value']
    meanPop *= popScale
    logger.info('Generating the population for %s (%s freeze-dried people)'
                % (descr['abbrev'], meanPop))
    for i in xrange(int(round(meanPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.HOME)
        assert ward is not None, 'Ran out of beds populating %(abbrev)s!' % descr
        agent = PatientAgent('PatientAgent_HOME_%s_%d' % (ward._name, i), patch, ward)
        agent.reHome(fac.manager.patch)
        agent.setStatus(homeAddr=findQueueForTier(CareTier.HOME, fac.reqQueues).getGblAddr())
        ward.handlePatientArrival(agent, None)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, agent),
                              None)
    return []

TOTAL_COMMUNITY = 0

def generateFull(facilityDescr, patch, policyClasses=None, categoryNameMapper=None,
                 communityClass=None):
    global cacheVer
    makeLMDBDirs()
    global TOTAL_COMMUNITY

    if communityClass is None:
        communityClass = Community
    fac = communityClass(facilityDescr, patch, policyClasses=policyClasses,
                         categoryNameMapper=categoryNameMapper)


    # create cache dir    facilityDescr['abbrev']
    # create facility dir for this instance
    # integrate everything!!!

    wards = fac.getWards()
    if len(wards) != 1:
        raise RuntimeError("caching wards assumes 1 ward per community")
    ward = wards[0]

    orig = ward.orig
    changed = ward.changed
    # get myself a freezer to access the metadata
    tFreezer = CopyOnWriteLMDBFreezer(ward, orig, changed, ward.infoList, False)

    try:
        ward.infoList = tFreezer.getInfoList()
        wardInfo = tFreezer.getWardData()

        if wardInfo[0] != cacheVer:
            raise InvalidCacheVer()

        ver, fDesc,freezerList = wardInfo
        if fDesc != facilityDescr:
            raise Exception()

        agentCount = 0
        for freezerData in freezerList:
            pType, loggerName = freezerData
            freezer = ward.freezers[pType]
            freezer.frozenAgentLoggerName = loggerName
            agentCount += len(freezer.frozenAgentList)

        if agentCount < 1:
            raise Exception()

        TOTAL_COMMUNITY += agentCount
        PatientAgent.allocateIds(fac, agentCount)
        logger.info('read population for %s from cache (%s freeze-dried people, %s)'
                    %(fDesc['abbrev'], agentCount, TOTAL_COMMUNITY))

        return [fac], fac.getWards(), []

    except InvalidCacheVer:
        raise RuntimeError("community cache versions are out of sync within the same file!")
        #ward.orig.close()
        #oName = origFname(facilityDescr['abbrev'])
        #ward.orig = interdict.InterDict(oName, convert_int=True,
        #                                val_serialization='pickle', overwrite_existing=True)
        #InterdictMapping[oName] = ward.orig
    except:
        #import traceback
        #traceback.print_exc()
        # no cache file or whatever was there wasn't satisfactory for any reason.
        # since if we found a valid cache file we'd have returned already, just
        # fall out of the exception.
        pass

    # at this point we _should_ wipe our section of the interdict
    # it's not absolutely necessary but it's a big todo
    #  *** TODO wipe our section of the interdict using a new method keyRange()

    fac.patientCacheIsBeingRegenerated = True

    pop = _populate(fac, facilityDescr, patch)

    return [fac], fac.getWards(), pop


def estimateWork(facRec):
    return 1  # Because everyone is freeze-dried


def checkSchema(facilityDescr):
    global _validator
    if _validator is None:
        _validator = schemautils.getValidator(_schema)
    nErrors = sum([1 for e in _validator.iter_errors(facilityDescr)])  # @UnusedVariable
    return nErrors

################
# Constants are loaded when the client module calls importCommunity()
################
