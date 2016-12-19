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

import random
from scipy.stats import expon, binom
import logging
import types

import pyrheabase
import pyrheautils
from phacsl.utils.collections.phacollections import enum
from facilitybase import DiagClassA, CareTier, TreatmentProtocol, BirthQueue, HOMEQueue
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent, PatientStatusSetter
from facilitybase import PatientStatus, PatientDiagnosis, FacilityManager, PthStatus
from quilt.netinterface import GblAddr
from stats import CachedCDFGenerator, BayesTree
from hospital import ClassASetter
import schemautils

logger = logging.getLogger(__name__)

category = 'COMMUNITY'
_schema = 'communityfacts_schema.yaml'
_validator = None
_constants_values = '$(MODELDIR)/constants/community_constants.yaml'
_constants_schema = 'community_constants_schema.yaml'
_constants = None

PatientCategory = enum('BASE', 'COLONIZED')

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


class CommunityWard(Ward):
    """This 'ward' type represents being out in the community"""

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
                raise CommunityWard.FreezerError('%s cannot freezedry %s because it is still in motion'
                                                 % (self.ward._name, d['name']))
            if self.frozenAgentLoggerName is None:
                self.frozenAgentLoggerName = d['loggerName']
                del d['loggerName']
            elif d['loggerName'] == self.frozenAgentLoggerName:
                del d['loggerName']
            else:
                raise CommunityWard.FreezerError('%s cannot freezedry %s because it has the wrong logger'
                                                 % (self.ward._name, d['name']))
            typeTpl, linL, valL = lencode(d)  # @UnusedVariable
    #         print 'dictionary: %s' % str(d)
    #         print 'types: %s' % str(typeTpl)
    #         print 'linearized: %s' % str(linL)
    #         print 'values: %s' % str(valL)
            if self.frozenAgentTypePattern is None:
                self.frozenAgentTypePattern = typeTpl
            elif typeTpl != self.frozenAgentTypePattern:
                raise CommunityWard.FreezerError('%s cannot freezedry %s because it has the wrong type pattern'
                                                 % (self.ward._name, d['name']))
            self.frozenAgentList.append(tuple(valL))

        def removeAndThaw(self, frozenAgent):
            self.frozenAgentList.remove(frozenAgent)
            valL = list(frozenAgent)
            d, leftovers = ldecode(self.frozenAgentTypePattern, valL)
            assert not leftovers, ('%s had %s left over unfreezing agent'
                                   % (self.ward._name, leftovers))
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


    def __init__(self, name, patch, nBeds):
        Ward.__init__(self, name, patch, CareTier.HOME, nBeds=nBeds)
        self.checkInterval = 1  # so they get freeze dried promptly
        self.freezers = {}
        for patientCategory in PatientCategory.names.keys():
            self.freezers[patientCategory] = CommunityWard.Freezer(self)
        self.newArrivals = []

    def classify(self, agent):
        """Return the PatientCategory appropriate for this agent"""
        if agent._status.pthStatus == PthStatus.COLONIZED:
            return PatientCategory.COLONIZED
        elif agent._status.pthStatus == PthStatus.CLEAR:
            return PatientCategory.BASE
        else:
            raise CommunityWard.FreezerError('%s has unexpected PthStatus %s' %
                                             (agent.name, PthStatus.names[agent._status.pthStatus]))

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

    def perTickActions(self, timeNow):
        dT = timeNow - self.fac.collectiveStatusStartDate
        for ward in self.fac.getWards():
            if dT != 0:
                for patCat, freezer in ward.freezers.items():
                    assert patCat in self.fac.cachedCDFs, ('%s has no CDF for %s' %
                                                           (self.fac.name, PatientCategory.names[patCat]))
                    pThaw = self.fac.cachedCDFs[patCat].intervalProb(0, dT)
                    nFroz = len(freezer.frozenAgentList)
                    nThawed = binom.rvs(nFroz, pThaw)
                    changedList = random.sample(freezer.frozenAgentList, nThawed)[:]
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
    def __init__(self, descr, patch, policyClasses=None, categoryNameMapper=None):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr,
                          descr, patch,
                          reqQueueClasses=[pyrheabase.FacRequestQueue, BirthQueue, HOMEQueue],
                          policyClasses=policyClasses,
                          managerClass=CommunityManager,
                          categoryNameMapper=categoryNameMapper)
        descr = self.mapDescrFields(descr)
        meanPop = descr['meanPop']['value']
        nBeds = int(round(3.0*meanPop))
        losModel = _constants['communityLOSModel']
        assert losModel['pdf'] == 'expon(lambda=$0)', \
            "Unexpected losModel form %s for %s!" % (losModel['pdf'], descr['abbrev'])
        self.addWard(CommunityWard('%s_%s_%s' % (category, patch.name, descr['abbrev']),
                                   patch, nBeds))
        baseRate = losModel['parms'][0]
        colonizedRate = _constants['communityColonizedReadmissionRateScale']['value'] * baseRate
        self.cachedCDFs = {PatientCategory.BASE: CachedCDFGenerator(expon(scale=1.0/baseRate)),
                           PatientCategory.COLONIZED: CachedCDFGenerator(expon(scale=1.0/colonizedRate))}
        self.collectiveStatusStartDate = 0
        self.treeCache = {}

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
            sickRate = 1.0 - (deathRate + verySickRate + needsLTACRate + needsRehabRate)
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
                                                      ], tag='FATE'),
                             PatientStatusSetter(),
                             changeProb, tag='LOS')
            self.treeCache[key] = tree
            return tree

    def prescribe(self, patientDiagnosis, patientTreatment):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
                return (CareTier.HOME, TreatmentProtocol.NORMAL)
            elif patientDiagnosis.overall == PatientOverallHealth.FRAIL:
                return (CareTier.NURSING, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
            return (CareTier.NURSING, TreatmentProtocol.REHAB)
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSLTAC:
            return (CareTier.LTAC, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.SICK:
            return (CareTier.HOSP, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.VERYSICK:
            return (CareTier.ICU, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.DEATH:
            return (None, TreatmentProtocol.NORMAL)
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
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              0)
        ward.flushNewArrivals()  # since we are about to manually freeze-dry.
        ward.freezers[ward.classify(a)].freezeAndStore(a)
    return []


def generateFull(facilityDescr, patch, policyClasses=None, categoryNameMapper=None):
    fac = Community(facilityDescr, patch, policyClasses=policyClasses,
                    categoryNameMapper=categoryNameMapper)
    return [fac], fac.getWards(), _populate(fac, facilityDescr, patch)


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
