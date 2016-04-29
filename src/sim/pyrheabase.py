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

import logging
from random import shuffle
import patches

logger = logging.getLogger(__name__)


class Location(patches.MultiInteractant):
    def __init__(self, name, patch, nCapacity, checkInterval=1):
        patches.MultiInteractant.__init__(self, name, nCapacity, patch)
        self.checkInterval = checkInterval

    def getDepartureMsgPayload(self, person):
        return None

    def getArrivalMsgPayload(self, person):
        return None

    def getReqQueueAddr(self):
        raise RuntimeError('Derived classes must implement this method')


class HoldQueue(patches.Interactant):
    def __init__(self, name, patch, debug=False):
        patches.Interactant.__init__(self, name, patch, debug=debug)
        self.keyCounter = 0
        self.heldDict = {}

    def getUniqueKey(self):
        key = self.keyCounter
        self.keyCounter += 1
        return key

    def lock(self, lockingAgent, key=None, debug=False):
        if key is not None:
            self.heldDict[key] = lockingAgent
        return patches.Interactant.lock(self, lockingAgent)

    def awaken(self, agentOrKey):
        if isinstance(agentOrKey, patches.Agent):
            agent = agentOrKey
        elif agentOrKey in self.heldDict:
            agent = self.heldDict[agentOrKey]
            del self.heldDict[agentOrKey]
        else:
            raise RuntimeError("%s does not know the agent or key %s" % (self._name, agentOrKey))
        return patches.Interactant.awaken(self, agent)


class Manager(patches.Agent):
    def __init__(self, name, patch, managementBase):
        patches.Agent.__init__(self, name, patch)
        self.timeless = True
        self.toManage = managementBase
        self.logger = logger.getChild('Manager')

    def handleRequest(self, req, logDebug, timeNow):
        if isinstance(req, SimpleMsg):
            return self.fac.handleIncomingMsg(req.__class__, req.payload, timeNow)
        else:
            raise RuntimeError("%s unexpectedly got the message %s" % (self.name, req.name))

    def run(self, startTime):
        timeNow = startTime  # @UnusedVariable
        logDebug = self.logger.isEnabledFor(logging.DEBUG)
        while True:
            foundAny = True
            while foundAny:
                foundAny = False
                for rQ in self.fac.reqQueues:
                    if rQ._lockQueue:
                        foundAny = True
                        req = rQ._lockQueue[0]
                        timeNow = self.handleRequest(req, logDebug, timeNow)
                        rQ.awaken(req)
            timeNow = self.sleep(0)  # @UnusedVariable


class SimpleMsg(patches.Agent):
    STATE_MOVING = 0
    STATE_ARRIVED = 1

    def __init__(self, name, patch, payload, destAddr, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.payload = payload
        self.destAddr = destAddr
        self.fsmstate = self.STATE_MOVING

    def run(self, startTime):
        timeNow = startTime  # @UnusedVariable
        while True:
            if self.fsmstate == self.STATE_MOVING:
                addr, final = self.patch.getPathTo(self.destAddr)
                if final:
                    self.fsmstate = self.STATE_ARRIVED
                timeNow = addr.lock(self)  # @UnusedVariable
            elif self.fsmstate == self.STATE_ARRIVED:
                break  # we are done

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['payload'] = self.payload
        d['fsmstate'] = self.fsmstate
        d['destAddr'] = self.destAddr
        return d

    def __setstate__(self, stateDict):
        patches.Agent.__setstate__(self, stateDict)
        self.payload = stateDict['payload']
        self.fsmstate = stateDict['fsmstate']
        self.destAddr = stateDict['destAddr']


class ArrivalMsg(SimpleMsg):
    pass


class DepartureMsg(SimpleMsg):
    pass


class RequestQueue(patches.Interactant):
    pass


class ManagementBase(object):
    def __init__(self, name, patch, managerClass=None, reqQueueClasses=None):
        if managerClass is None:
            managerClass = Manager
        if reqQueueClasses is None:
            reqQueueClasses = [RequestQueue]
        self.name = name
        self.manager = managerClass(name + '_Mgr', patch, self)
        self.reqQueues = [rQC(name+'_rQ', patch) for rQC in reqQueueClasses]
        for rQ in self.reqQueues:
            rQ.lock(self.manager)
        self.holdQueue = HoldQueue(name+'_hQ', patch)
        self.holdQueue.lock(self.manager)

    def getMsgPayload(self, msgType, person):
        if issubclass(msgType, ArrivalMsg):
            return person.loc.getArrivalMsgPayload(person)
        elif issubclass(msgType, DepartureMsg):
            return person.loc.getDepartureMsgPayload(person)
        else:
            raise RuntimeError('%s: payload request for unknown message type %s'
                               % (self.name, msgType.__name__))

    def handleIncomingMsg(self, msgType, payload, timeNow):
        return timeNow


class Person(patches.Agent):
    STATE_ATLOC = 0
    STATE_MOVING = 1
    STATE_JUSTARRIVED = 2

    def __init__(self, name, patch, loc, debug=False):
        """
        loc is the current location, an interactant to which the Person instance is locked.
        """
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.fsmstate = PatientAgent.STATE_ATLOC
        self.loc = loc
        self.locAddr = loc.getGblAddr()
        self.newLocAddr = None
        self.logger = logging.getLogger(__name__ + '.Person')

    def getPostArrivalPauseTime(self, timeNow):
        """
        This allows the insertion of an extra pause on arrival at a new location.  It is useful
        for desynchronizing the activity cycles of agents at a location which does not sample
        all agents daily.
        """
        return 0

    def getNewLocAddr(self, timeNow):
        """
        This method is called once each time the Person agent is active and returns a tuple
        of the form (newLocGlobalAddr, updatedTimeNow).  Returning a newLocGlobalAddr of
        self.locAddr (that is, the current value) indicates that the Person stays attached
        to the same location for this time slice.  Returning a newLocGlobalAddr
        of None signals 'death' and will result in self.handleDeath being called and the agent's
        thread exiting.  This method may not return for a long time as the agent waits for
        a new location to become available.
        """
        return self.locAddr, timeNow

    def handleArrival(self, timeNow):
        """
        An opportunity to do bookkeeping on arrival at self.loc.  The Person agent has already
        locked the interactant self.loc.
        """
        pass

    def handleDeath(self, timeNow):
        """
        The name says it all.  Do any bookkeeping needed to deal with the death of this Person.
        After this method returns, a DepartureMessage will be sent to inform its current location
        of departure and the agent's run method will exit.
        """
        pass

    def run(self, startTime):
        timeNow = startTime
        timeNow = self.sleep(self.getPostArrivalPauseTime(timeNow))
        while True:
            if self.debug:
                self.logger.debug('%s point 0: state %s day %s' %
                                  (self.name, self.fsmstate, timeNow))
            if self.fsmstate == Person.STATE_ATLOC:
                newLocAddr, timeNow = self.getNewLocAddr(timeNow)
                if newLocAddr is None:
                    self.handleDeath(timeNow)
                    self.patch.launch(DepartureMsg(self.name + '_depMsg',
                                                   self.patch,
                                                   self.loc.getDepartureMsgPayload(self),
                                                   self.loc.getReqQueueAddr()),
                                      timeNow)
                    timeNow = self.loc.unlock(self)
                    break
                elif newLocAddr == self.locAddr:
                    if self.debug:
                        self.logger.debug('%s point 9: sleeping for %s' %
                                          (self.name, self.loc.checkInterval))
                    timeNow = self.sleep(self.loc.checkInterval)
                else:
                    if self.debug:
                        self.logger.debug('%s point 1: new addr %s vs current %s' %
                                          (self.name, self.newLocAddr, self.locAddr))
                    self.patch.launch(DepartureMsg(self.name + '_depMsg',
                                                   self.patch,
                                                   self.loc.getDepartureMsgPayload(self),
                                                   self.loc.getReqQueueAddr()),
                                      timeNow)
                    self.newLocAddr = newLocAddr
                    timeNow = self.loc.unlock(self)
                    self.fsmstate = Person.STATE_MOVING

            elif self.fsmstate == Person.STATE_MOVING:
                self.loc = None
                addr, final = self.patch.getPathTo(self.newLocAddr)
                if final:
                    self.fsmstate = Person.STATE_JUSTARRIVED
                    self.loc = addr
                    if self.debug:
                        self.logger.debug('%s point 6: arrive tier %s day %s'
                                          % (self.name, self.ward.tier, timeNow))
                timeNow = addr.lock(self)
            elif self.fsmstate == Person.STATE_JUSTARRIVED:
                self.handleArrival(timeNow)
                self.patch.launch(ArrivalMsg(self.name + '_arvMsg',
                                             self.patch,
                                             self.loc.getArrivalMessagePayload(),
                                             self.loc.getReqQueueAddr()),
                                  timeNow)
                self.fsmstate = PatientAgent.STATE_ATLOC
                if self.debug:
                    self.logger.debug('%s point 7: day %s'
                                      % (self.name, timeNow))
                timeNow = self.sleep(self.getPostArrivalPauseTime(timeNow))
            if self.debug:
                self.logger.debug('%s point 2' % self.name)

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['newLocAddr'] = self.newLocAddr
        d['fsmstate'] = self.fsmstate
        return d

    def __setstate__(self, d):
        patches.Agent.__setstate__(self, d)
        self.newLocAddr = d['newLocAddr']
        self.fsmstate = d['fsmstate']


class Ward(Location):
    def __init__(self, name, patch, tier, nBeds):
        Location.__init__(self, name, patch, nBeds)
        self.tier = tier
        self.fac = None

    def getDepartureMsgPayload(self, person):
        return self.fac.getMsgPayload(DepartureMsg, person)

    def getArrivalMsgPayload(self, person):
        return self.fac.getMsgPayload(ArrivalMsg, self)

    def getReqQueueAddr(self):
        return self.fac.reqQueues[0].getGblAddr()


class FacRequestQueue(RequestQueue):
    pass


class FacilityManager(Manager):
    def __init__(self, name, patch, facility):
        super(FacilityManager, self).__init__(name, patch, facility)
        self.logger = logger.getChild('FacilityManager')

    @property
    def fac(self):
        return self.toManage

    def allocateAvailableBed(self, tier):
        if tier in self.fac.wardTierDict:
            for w in self.fac.wardTierDict[tier]:
                addr = w.getGblAddr().getLclAddr()
                ww, n = self.fac.wardAddrDict[addr]  # @UnusedVariable
                if n > 0:
                    self.fac.wardAddrDict[addr] = (w, n-1)
                    # print 'alloc: ward %s now has %d free beds of %d' % (w._name, n-1, w._nLocks)
                    return w
            return None
        else:
            return None

    def handleBedRequest(self, req, logDebug, timeNow):
        payload = req.payload
        tier = req.tier
        ward = self.allocateAvailableBed(tier)
        req.payload = self.fac.handleBedRequestResponse(ward, payload, timeNow)
        if ward is None:
            if logDebug:
                if tier in self.fac.wardTierDict:
                    self.logger.debug('%s: no beds available for %s' % (self.name, req.name))
                else:
                    self.logger.debug('%s: no such ward for %s' % (self.name, req.name))
            req.fsmstate = BedRequest.STATE_DENIEDWARD
        else:
            req.fsmstate = BedRequest.STATE_GOTWARD
            req.bedWard = ward.getGblAddr()
            if logDebug:
                self.logger.debug('%s: found a bed for %s' % (self.name, req.name))
        return timeNow

    def handleRequest(self, req, logDebug, timeNow):
        if isinstance(req, BedRequest):
            return self.handleBedRequest(req, logDebug, timeNow)
        else:
            return super(FacilityManager, self).handleRequest(req, logDebug, timeNow)


class Facility(ManagementBase):
    def __init__(self, name, patch, managerClass=None, reqQueueClasses=None):
        if managerClass is None:
            managerClass = FacilityManager
        if reqQueueClasses is None:
            reqQueueClasses = [FacRequestQueue]
        ManagementBase.__init__(self, name, patch, managerClass, reqQueueClasses)
        self.wardTierDict = {}  # lists of wards by tier
        self.wardAddrDict = {}  # dict of (ward, nFree) by lclAddr

    def addWard(self, ward):
        if ward.tier not in self.wardTierDict:
            self.wardTierDict[ward.tier] = []
        self.wardTierDict[ward.tier].append(ward)
        self.wardAddrDict[ward.getGblAddr().getLclAddr()] = (ward, ward.nFree)
        ward.fac = self
        return ward

    def getWards(self, tier=None):
        if tier:
            if tier in self.wardTierDict:
                return self.wardTierDict[tier]
            else:
                return []
        else:
            return [w for w, n in self.wardAddrDict.values()]  # @UnusedVariable

    def getTiers(self):
        return self.wardTierDict.keys()

    def getMsgPayload(self, msgType, patientAgent):
        if issubclass(msgType, ArrivalMsg):
            return patientAgent.ward.getGblAddr()
        elif issubclass(msgType, DepartureMsg):
            return patientAgent.ward.getGblAddr()
        else:
            return super(Facility, self).getMsgPayload(msgType, patientAgent)

    def handleIncomingMsg(self, msgType, payload, timeNow):
        wardGblAddr = payload
        addr = wardGblAddr.getLclAddr()
        try:
            ward, nFree = self.wardAddrDict[addr]  # @UnusedVariable
        except:
            raise RuntimeError("%s: I do not own the ward %s" % (self.name, ward._name))
        if issubclass(msgType, ArrivalMsg):
            pass
        elif issubclass(msgType, DepartureMsg):
            self.wardAddrDict[addr] = (ward, nFree+1)
        else:
            raise RuntimeError('%s: got unknown message type %s' % (self.name,
                                                                    msgType .__name__))
        return timeNow

    def getBedRequestPayload(self, patientAgent, desiredTier):
        return None

    def handleBedRequestResponse(self, ward, payload, timeNow):
        """
        This routine is called in the time slice of the Facility Manager when the manager
        responds to a request for a bed.  If the request was denied, 'ward' will be None.
        The return value of this message becomes the new payload.
        """
        return None

    def handleBedRequestFate(self, ward, payload, timeNow):
        """
        This routine is called in the time slice of the Facility Manager when the manager
        responds to a request for a bed.  If the search for a bed failed, 'ward' will be None.
        """
        pass


class BedRequest(patches.Agent):
    STATE_START = 0
    STATE_ASKWARD = 1
    STATE_DENIEDWARD = 2
    STATE_GOTWARD = 3
    STATE_FAILED = 5
    STATE_MOVING = 6

    def __init__(self, name, patch, tier, homeWardAddr, patientKey,
                 facilityOptions, payload, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.tier = tier                  # the needed CareTier
        self.homeWardAddr = homeWardAddr  # to find our way home at the end of the search
        self.patientKey = patientKey      # to awaken the originating PatientAgent
        self.facilityOptions = facilityOptions[:]  # candidate facilities, in preference order
        self.facilityOptions.reverse()
        self.payload = payload            # useful for derived classes
        self.bedWard = None               # the ward satisfying the request
        self.dest = None                  # the current travel destination
        self.fsmstate = BedRequest.STATE_START

    def run(self, startTime):
        timeNow = startTime
        while True:
            if self.fsmstate == BedRequest.STATE_START:
                if len(self.facilityOptions) == 0:
                    self.fsmstate = BedRequest.STATE_FAILED
                    self.dest = None
                else:
                    self.fsmstate = BedRequest.STATE_MOVING
                    self.dest = self.facilityOptions.pop()
            elif self.fsmstate == BedRequest.STATE_MOVING:
                    addr, final = self.patch.getPathTo(self.dest)
                    if final:
                        self.fsmstate = BedRequest.STATE_ASKWARD
                    timeNow = addr.lock(self)
            elif self.fsmstate == BedRequest.STATE_ASKWARD:
                raise RuntimeError('%s: I SHOULD BE ASLEEP at time %s' % (self.name, timeNow))
            elif self.fsmstate == BedRequest.STATE_GOTWARD:
                addr, final = self.patch.getPathTo(self.homeWardAddr)
                if final:
                    oldWard = addr
                    oldWard.fac.handleBedRequestFate(self.bedWard, self.payload, timeNow)
                    patientAgent = oldWard.fac.holdQueue.awaken(self.patientKey)
                    patientAgent.newWardAddr = self.bedWard
                    break
                timeNow = addr.lock(self)
            elif self.fsmstate == BedRequest.STATE_DENIEDWARD:
                self.fsmstate = BedRequest.STATE_START
            elif self.fsmstate == BedRequest.STATE_FAILED:
                addr, final = self.patch.getPathTo(self.homeWardAddr)
                if final:
                    oldWard = addr
                    oldWard.fac.handleBedRequestFate(None, self.payload, timeNow)
                    patientAgent = oldWard.fac.holdQueue.awaken(self.patientKey)
                    patientAgent.newWardAddr = None
                    break
                timeNow = addr.lock(self)

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['homeWardAddr'] = self.homeWardAddr
        d['facilityOptions'] = self.facilityOptions
        d['tier'] = self.tier
        d['fsmstate'] = self.fsmstate
        d['bedWard'] = self.bedWard
        d['patientKey'] = self.patientKey
        d['dest'] = self.dest
        d['payload'] = self.payload
        return d

    def __setstate__(self, stateDict):
        patches.Agent.__setstate__(self, stateDict)
        self.homeWardAddr = stateDict['homeWardAddr']
        self.facilityOptions = stateDict['facilityOptions']
        self.tier = stateDict['tier']
        self.fsmstate = stateDict['fsmstate']
        self.bedWard = stateDict['bedWard']
        self.patientKey = stateDict['patientKey']
        self.dest = stateDict['dest']
        self.payload = stateDict['payload']


class PatientAgent(Person):

    def __init__(self, name, patch, ward, timeNow=0, debug=False):
        Person.__init__(self, name, patch, ward, debug=debug)
        self.tier = ward.tier
        self.logger = logging.getLogger(__name__ + '.PatientAgent')

    @property
    def ward(self):
        return self.loc

    def handleTierUpdate(self, timeNow):
        return self.ward.tier

    def handleDeath(self, timeNow):
        pass

    def getCandidateFacilityList(self, timeNow, newTier):
        facAddrList = [tpl[1] for tpl in self.patch.serviceLookup('FacRequestQueue')]
        shuffle(facAddrList)
        return facAddrList

    def getNewLocAddr(self, timeNow):
        tier = self.handleTierUpdate(timeNow)
        if tier == self.tier:
            return self.locAddr, timeNow  # things stay the same
        elif tier is None:
            return None, timeNow  # signal death
        else:
            while True:
                facAddrList = self.getCandidateFacilityList(timeNow, tier)
                key = self.ward.fac.holdQueue.getUniqueKey()
                self.patch.launch(BedRequest(self.name + '_bedReq', self.patch,
                                             tier, self.ward.getGblAddr(),
                                             key, facAddrList,
                                             self.ward.fac.getBedRequestPayload(self, tier)),
                                   timeNow)
                if self.debug:
                    self.logger.debug('%s point 3: launched req for tier %s'
                                      % (self.name, tier))
                timeNow = self.ward.fac.holdQueue.lock(self, key=key)
                if self.debug:
                    self.logger.debug('%s point 4: day %s'
                                      % (self.name, timeNow))
                # The BedRequest will have updated the value of self.newLocAddr before exiting.
                # Just be sure not to return it if the new value is None, or the calling loop
                # will interpret it as a signal that the agent has died!
                if self.newLocAddr is None:
                    # Nowhere to go; try again tomorrow
                    # print ('%s is stuck at %s; going back to sleep at %s' %
                    #        (self.name, self.ward, timeNow))
                    if self.debug:
                        self.logger.debug('%s point 8: stuck here day %s'
                                          % (self.name, timeNow))
                    timeNow = self.sleep(1)
                    # state is unchanged
                else:
                    return self.newLocAddr, timeNow

    def handleArrival(self, timeNow):
        """
        An opportunity to do bookkeeping on arrival at self.loc.  The Person agent has already
        locked the interactant self.loc.
        """
        self.tier = self.loc.tier

    def __getstate__(self):
        d = Person.__getstate__(self)
        d['tier'] = self.tier
        return d

    def __setstate__(self, d):
        Person.__setstate__(self, d)
        self.tier = d['tier']
