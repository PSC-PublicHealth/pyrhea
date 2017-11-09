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
from scipy.stats import expon, binom
import logging
import types
from collections import defaultdict

from phacsl.utils.collections.phacollections import DefaultDict
import phacsl.utils.collections.interdict as interdict

import cPickle as pickle
import gzip

import pyrheabase
import pyrheautils
from facilitybase import DiagClassA, CareTier, TreatmentProtocol, BirthQueue, HOMEQueue
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent, PatientStatusSetter
from facilitybase import ClassASetter, PatientStatus, PatientDiagnosis, FacilityManager
from quilt.netinterface import GblAddr
from stats import CachedCDFGenerator, BayesTree
import schemautils

logger = logging.getLogger(__name__)

category = 'COMMUNITY'
_schema = 'communityfacts_schema.yaml'
_validator = None
_constants_values = '$(MODELDIR)/constants/community_constants.yaml'
_constants_schema = 'community_constants_schema.yaml'
_constants = None




def importCommunity(moduleName):
    attrs = (
        'category',
        '_schema',
        '_validator',
        '_constants_values',
        '_constants_schema',
        '_constants',
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


class FreezerError(RuntimeError):
    pass

class Freezer(object):
    def __init__(self, ward):
        self.ward = ward
        self.frozenAgentClass = PatientAgent
        self.frozenAgentLoggerName = None
        self.frozenAgentTypePattern = None
        self.frozenAgentList = []

    def freezeAndStore(self, agent):
        if agent in self.ward._lockingAgentList:
            self.ward._lockingAgentList.remove(agent)
        self.ward.suspend(agent)
        assert agent not in self.ward._lockingAgentList, 'It is still there!'
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
        if 0:
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
        if 0:
            self.frozenAgentList.append(tuple(valL))
        else:
            self.frozenAgentList.append(valL)

    def removeAndThaw(self, frozenAgent):
        self.frozenAgentList.remove(frozenAgent)
        if 0:
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
        self.ward._lockingAgentList.append(agent)
        return agent

class InvalidCacheVer(Exception):
    pass

def mapLMDBAbbrev(abbrev):
    return 'all'

def mapLMDBSegments(abbrev, iDict):
    cacheVer = 1
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
        This is a freezer for community agents that is stored in an LMDB InterDict rather than in memory.
        It is stored in two parts so that the community agents can be cached and then changed community
        members are written to a second LMDB file.
        
        A given ward can have multiple freezers to store multiple types of community agents but it's useful
        to have all freezers in one pair of interdicts.  To this end there is an infoPointer which is a list
        of values that help the multiple freezers interact.

        ward: this is the ward that is associated with the freezer
        orig: this is the InterDict where the cached community agents are stored.  By default this is not modified.
        changed: any newly stored agents are saved here
        infoList: this is a list shared between all freezers using the same pair of interdicts with the following values:
            origRO: if True (normal usage) write new agents to changed.  For building the cache set to False
                    and all new members will be written to orig
            maxOrig: when looking up an agent, any ID less than this will be found in orig, otherwise found in changed
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
        if agent in self.ward._lockingAgentList:
            self.ward._lockingAgentList.remove(agent)
        self.ward.suspend(agent)
        assert agent not in self.ward._lockingAgentList, 'It is still there!'
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

        
    def removeAndThaw(self, frozenAgent):
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
        self.ward._lockingAgentList.append(agent)
        return agent

    def saveFreezerData(self):
        origRO, maxOrig, nextId = self.infoList
        if origRO:
            lDict = self.changed
        else:
            lDict = self.orig
        try:
            agentListDict = self.changed[self._AgentListId+self.iDictOffset]
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
    return pyrheautils.pathTranslate('$(COMMUNITYCACHEDIR)/%s'%abbrev)

def changedFname(abbrev):
    abbrev = mapLMDBAbbrev(abbrev)
    return pyrheautils.pathTranslate('$(AGENTDIR)/%s'%abbrev)

InterdictMapping = {}


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
        if oName in InterdictMapping:
            self.orig = InterdictMapping[oName]
        else:
            self.orig = interdict.InterDict(oName, convert_int=True, val_serialization='pickle')
            InterdictMapping[oName] = self.orig

        cName = changedFname(abbrev)
        if cName in InterdictMapping:
            self.changed = InterdictMapping[cName]
        else:
            self.changed = interdict.InterDict(cName, overwrite_existing=True, convert_int=True, val_serialization='pickle')
            InterdictMapping[cName] = self.changed

        try:
            self.iDictOffset = mapLMDBSegments(abbrev, self.orig)
        except InvalidCacheVer:
            self.orig.close()
            self.orig = interdict.InterDict(oName, overwrite_existing=True, convert_int=True, val_serialization='pickle')
            self.iDictOffset = mapLMDBSegments(abbrev, self.orig)

        self.infoList = [True,1,1]  # we'll overwrite this shortly.  This deals with a chicken/egg problem
        self.freezers = DefaultDict(lambda dd,key: newFreezer(dd, self, self.orig, self.changed, self.infoList, key))
        self.newArrivals = []

    def classify(self, agent):
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


class CommunityManager(FacilityManager):

    def getProbThaw(self, patientCategory, dT):
        return self.fac.cachedCDFs[patientCategory].intervalProb(0, dT)

    def perTickActions(self, timeNow):
        dT = timeNow - self.fac.collectiveStatusStartDate
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
                        print 'ValueError for %s: nFroz=%s, pThaw=%s' % (self.fac.name, nFroz, pThaw)
                        raise
                    self.fac.getNoteHolder().addNote({('thawed_%s' % timeNow) : nThawed })
                    changedList = random.sample(freezer.frozenAgentList, nThawed)
                    for a in changedList:
                        thawedAgent = freezer.removeAndThaw(a)
                        if thawedAgent.debug:
                            thawedAgent.logger.debug('%s unfreezedried %s at %s'
                                                     % (ward._name, thawedAgent.name, timeNow))
            for agent in ward.newArrivals:
                if agent.debug:
                    agent.logger.debug('%s freezedrying %s at %s'
                                       % (ward._name, agent.name, timeNow))
                ward.freezers[ward.classify(agent)].freezeAndStore(agent)
            ward.newArrivals = []
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
                 managerClass=None):
        if managerClass is None:
            managerClass = CommunityManager
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr,
                          descr, patch,
                          reqQueueClasses=[pyrheabase.FacRequestQueue, BirthQueue, HOMEQueue],
                          policyClasses=policyClasses,
                          managerClass=managerClass,
                          categoryNameMapper=categoryNameMapper)
        descr = self.mapDescrFields(descr)
        meanPop = descr['meanPop']['value']
        nBeds = int(round(3.0*meanPop))
        losModel = _constants['communityLOSModel']
        assert losModel['pdf'] == 'expon(lambda=$0)', \
            "Unexpected losModel form %s for %s!" % (losModel['pdf'], descr['abbrev'])
        self.addWard(CommunityWard(descr['abbrev'], '%s_%s_%s' % (category, patch.name, descr['abbrev']),
                                   patch, nBeds))
        self.setCDFs(losModel)
        self.collectiveStatusStartDate = 0
        self.treeCache = {}

    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        self.treeCache = {}        

    def setCDFs(self, losModel):
        baseRate = losModel['parms'][0]
        self.cachedCDFs = {"base": CachedCDFGenerator(expon(scale=1.0/baseRate))}
        

    def getStatusChangeTree(self, patientStatus, ward, treatment, startTime, timeNow):
        careTier = ward.tier
        assert careTier == CareTier.HOME, \
            "The community only offers CareTier 'HOME'; found %s" % careTier
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        if key in self.treeCache:
            return self.treeCache[key]
        else:
            changeProb = 1.0  # since the agent was already randomly chosen to be un-frozen
            deathRate = _constants['communityDeathRate']['value']
            verySickRate = _constants['communityVerySickRate']['value']
            needsLTACRate = _constants['communityNeedsLTACRate']['value']
            needsRehabRate = _constants['communityNeedsRehabRate']['value']
            needsSkilNrsRate = _constants['communityNeedsSkilNrsRate']['value']
            needsVentRate = (_constants['communityNeedsVentRate']['value']
                             if 'communityNeedsVentRate' in _constants else 0.0)
            sickRate = 1.0 - (deathRate + verySickRate + needsLTACRate + needsRehabRate
                              + needsSkilNrsRate + needsVentRate)
            tree = BayesTree(BayesTree.fromLinearCDF([(deathRate,
                                                       ClassASetter(DiagClassA.DEATH)),
                                                      (sickRate,
                                                       ClassASetter(DiagClassA.SICK)),
                                                      (verySickRate,
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

    def prescribe(self, ward, patientId, patientDiagnosis, patientTreatment,
                  timeNow=None):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
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


def _populate(fac, descr, patch):
    assert 'meanPop' in descr, \
        "Community description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = float(descr['meanPop']['value'])
    logger.info('Generating the population for %s (%s freeze-dried people)'
                % (descr['abbrev'], meanPop))
    for i in xrange(int(round(meanPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.HOME)
        assert ward is not None, 'Ran out of beds populating %(abbrev)s!' % descr
        a = PatientAgent('PatientAgent_HOME_%s_%d' % (ward._name, i), patch, ward)
        a.reHome(fac.manager.patch)
        a._status = a._status._replace(homeAddr=ward.getGblAddr())
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              0)
        ward.flushNewArrivals()  # since we are about to manually freeze-dry.
        ward.freezers[ward.classify(a)].freezeAndStore(a)
    return []

def generateFull(facilityDescr, patch, policyClasses=None, categoryNameMapper=None,
                 communityClass=None):
    cacheVer = 1
    makeLMDBDirs()

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
            raise invalidCacheVer()

        ver, fDesc,freezerList,patientDataDict,patientStats = wardInfo
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
        
        fac.patientDataDict = patientDataDict
        fac.patientStats = patientStats
        logger.info('read population for %s from cache (%s freeze-dried people)'%(fDesc['abbrev'], agentCount))

        return [fac], fac.getWards(), []

    except InvalidCacheVer as e:
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
    
    # Let's write all of our data at once to the interdict, so for now use a normal dict:
    orig = ward.orig
    ward.orig = {}
    ward.infoList = [False, ward.iDictOffset+5, ward.iDictOffset+5]
    
    pop = _populate(fac, facilityDescr, patch)

    freezerList = []
    for pType,freezer in ward.freezers.items():
        freezerList.append((pType,freezer.frozenAgentLoggerName))
        freezer.saveFreezerData()
    
    patientDataDict = fac.patientDataDict
    patientStats = fac.patientStats

    # since it's possible that there's no freezers because there are no agents, let's make another temp freezer
    tFreezer = CopyOnWriteLMDBFreezer(ward, orig, changed, ward.infoList, False)
    tFreezer.saveWardData((cacheVer,facilityDescr, freezerList, patientDataDict,patientStats))
    tFreezer.saveInfoList()
    
    orig.mset(ward.orig.items())
    orig.flush()

    for k in ward.orig.keys():
        del(ward.orig[k])
        
    ward.orig = orig
    
    return [fac], fac.getWards(), pop


def estimateWork(facRec):
    return 1  # Because everyone is freeze-dried


def checkSchema(facilityDescr):
    global _validator
    if _validator is None:
        _validator = schemautils.getValidator(_schema)
    nErrors = sum([1 for e in _validator.iter_errors(facilityDescr)])  # @UnusedVariable
    return nErrors


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(_constants_values,
                                         _constants_schema)

