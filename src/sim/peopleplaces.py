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
import patches

logger = logging.getLogger(__name__)


class Location(patches.MultiInteractant):
    def __init__(self, name, patch, nCapacity, checkInterval=1):
        super(Location, self).__init__(name, nCapacity, patch)
        self.checkInterval = checkInterval

    def getDepartureMsgPayload(self, person):
        return None

    def getArrivalMsgPayload(self, person):
        return None

    def getReqQueueAddr(self):
        raise RuntimeError('Derived classes must implement this method')


class HoldQueue(patches.Interactant):
    def __init__(self, name, patch, debug=False):
        super(HoldQueue, self).__init__(name, patch, debug=debug)
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
        super(Manager, self).__init__(name, patch)
        self.timeless = True
        self.toManage = managementBase
        self.logger = logger.getChild('Manager')

    def handleRequest(self, req, logDebug, timeNow):
        if isinstance(req, SimpleMsg):
            return self.toManage.handleIncomingMsg(req.__class__, req.payload, timeNow)
        else:
            raise RuntimeError("%s unexpectedly got the message %s" % (self.name, req.name))

    def run(self, startTime):
        timeNow = startTime  # @UnusedVariable
        logDebug = self.logger.isEnabledFor(logging.DEBUG)
        while True:
            foundAny = True
            while foundAny:
                foundAny = False
                for rQ in self.toManage.reqQueues:
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
        super(SimpleMsg, self).__init__(name, patch, debug=debug)
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
        super(ManagementBase, self).__init__()
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

    def getAllQueues(self):
        """
        Queues are interactants and must be added to the patch, so this method provides an
        easy-to-access list.
        """
        return [self.holdQueue] + self.reqQueues


class Person(patches.Agent):
    STATE_ATLOC = 0
    STATE_MOVING = 1
    STATE_JUSTARRIVED = 2

    def __init__(self, name, patch, loc, debug=False):
        """
        loc is the current location, an interactant to which the Person instance will be locked.
        """
        super(Person, self).__init__(name, patch, debug=debug)
        self.fsmstate = Person.STATE_ATLOC
        self.loc = loc
        self.loc.lock(self)
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

    def handleDeparture(self, timeNow):
        """
        An opportunity to do bookkeeping on departure from self.loc.  The Person has not yet
        unlocked self.loc .  This call happens even if departure is via death.
        """
        pass

    def handleDeath(self, timeNow):
        """
        The name says it all.  Do any bookkeeping needed to deal with the death of this Person.
        After this method returns, a DepartureMessage will be sent to inform its current location
        of departure and the agent's run method will exit.  handleDeparture will also be called
        for this event.
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
                if self.debug:
                    self.logger.debug('%s point 1: new addr %s vs current %s' %
                                      (self.name, newLocAddr, self.locAddr))
                if newLocAddr is None:
                    if self.debug:
                        self.logger.debug('%s point 5: agent dies' %
                                          self.name)
                    self.handleDeparture(timeNow)
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
                    self.patch.launch(DepartureMsg(self.name + '_depMsg',
                                                   self.patch,
                                                   self.loc.getDepartureMsgPayload(self),
                                                   self.loc.getReqQueueAddr()),
                                      timeNow)
                    self.newLocAddr = newLocAddr
                    self.handleDeparture(timeNow)
                    timeNow = self.loc.unlock(self)
                    self.fsmstate = Person.STATE_MOVING

            elif self.fsmstate == Person.STATE_MOVING:
                self.loc = None
                addr, final = self.patch.getPathTo(self.newLocAddr)
                if final:
                    self.fsmstate = Person.STATE_JUSTARRIVED
                    self.loc = addr
                    self.locAddr = self.newLocAddr
                    if self.debug:
                        self.logger.debug('%s point 6: arrive tier %s day %s'
                                          % (self.name, self.ward.tier, timeNow))
                timeNow = addr.lock(self)
            elif self.fsmstate == Person.STATE_JUSTARRIVED:
                self.handleArrival(timeNow)
                self.patch.launch(ArrivalMsg(self.name + '_arvMsg',
                                             self.patch,
                                             self.loc.getArrivalMsgPayload(self),
                                             self.loc.getReqQueueAddr()),
                                  timeNow)
                self.fsmstate = Person.STATE_ATLOC
                if self.debug:
                    self.logger.debug('%s point 7: day %s'
                                      % (self.name, timeNow))
                timeNow = self.sleep(self.getPostArrivalPauseTime(timeNow))
            if self.debug:
                self.logger.debug('%s point 2: state is now %s day %s' %
                                  (self.name, self.fsmstate, timeNow))

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['newLocAddr'] = self.newLocAddr
        d['fsmstate'] = self.fsmstate
        return d

    def __setstate__(self, d):
        patches.Agent.__setstate__(self, d)
        self.newLocAddr = d['newLocAddr']
        self.fsmstate = d['fsmstate']
