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

from collections import defaultdict
import logging
from greenlet import greenlet

import netinterface
# from pympler import tracker
import agent

logger = logging.getLogger(__name__)


def _substituteClockAgentBreakHook(clockAgent):
    """
    This serves to monkey patch the clock agent's main loop in the 'agent' module so that
    it does not yield its thread, as this module will register a callback for that loop
    that will do so.
    """
    return clockAgent.ownerLoop.sequencer.getTimeNow()

agent._clockAgentBreakHook = _substituteClockAgentBreakHook


"""
To Do:
-Add rollback chain support
-Implement rollback
"""


class MsgTypes():
    GATE = 0


def getCommWorld():
    """Provide easy access to the world to packages that don't want to know about the network"""
    return netinterface.getCommWorld()


class Interactant(agent.Interactant):
    """
    An interactant stays in the same place, and is located by agents via its global
    address.  Agents may only modify interactants if they have locked them.  If more than
    one agent locks the same interactant, all but the first are suspended until the first
    agent unlocks the interactant.  At that point the next agent to lock the interactant
    becomes active.  This process continues until the queue of suspended agents have all
    had their turn.  (However, see MultiInteractant).

    Interactants are sometimes used as queues by agents which lock them at the start of a
    simulation and never lock them.  Other agents which subsequently lock the interactant
    are suspended and enqueued.  The first agent can manipulate the suspended agents and
    then subsequently awaken them, removing them from the lock queue and returning them to
    the list of active agents.
    """
    def __init__(self, name, patch, debug=False):
        agent.Interactant.__init__(self, name, patch.loop, debug)
        self.patch = patch

    def getGblAddr(self):
        return self.patch.group.getGblAddr((self.patch.patchId, self.id))


class MultiInteractant(agent.MultiInteractant):
    """
    A MultiInteractant is a type of Interactant for which multiple agents can hold locks
    and still remain active.  As long as the count of available locks is not exceeded,
    agents can lock the multiinteractant and continue to remain active.  If more than
    count agents lock the multiinteractant, the surplus agents are suspended and enqueued.
    """
    def __init__(self, name, count, patch, debug=False):
        agent.MultiInteractant.__init__(self, name, count, patch.loop, debug)
        self.patch = patch

    def getGblAddr(self):
        return self.patch.group.getGblAddr((self.patch.patchId, self.id))


class Agent(agent.Agent):
    """
    Agents have run methods, and are repeatedly given time slices in which those run
    methods execute.  Agents can be mobile or fixed-location.  A fixed-location agent
    will typically lock one or more Interactants at the beginning of a simulation and
    hold those locks forever, performing some maintenance task on any agent which
    subsequently locks the interactant.  For example, each patch has a GateAgent which
    controls the operation of Gates to and from the patch.

    A mobile agent must support being shipped from processor to processor, and thus
    must be serializable.  The mobile agent thus must implement __getstate__() and
    __setstate__(), as described by the Pickle module.  Typically the run method of a
    mobile agent is a simple finite state machine, with the state explicitly preserved
    by __getstate__() and __setstate__().  Agents which rely on their greenlet nature
    to preserve their state between time slices cannot be mobile, because the
    underlying greenlet is not serializable.
    """
    def __init__(self, name, patch, debug=False):
        agent.Agent.__init__(self, name, patch.loop, debug=debug)
        self.patch = patch

    def reHome(self, newPatch):
        self.ownerLoop = newPatch.loop
        self.parent = newPatch.loop
        self.patch = newPatch
        # print '%s home is now %s' % (self, newPatch)


# class OmniClock(Agent):
#     def __init__(self, ownerPatch):
#         Agent.__init__(self, 'OmniClock', ownerPatch)
#         self.timeless = True
#         self.patch = ownerPatch
#         self.comm = ownerPatch.comm
#         self.minVal = None
#         self.maxVal = None
#         if self.comm.rank == 0:
#             self.mem = np.zeros(2, dtype=np.int64)
#             self.win = MPI.Win.Create(self.mem, disp_unit=1, comm=self.comm)
#         else:
#             self.mem = None
#             self.win = MPI.Win.Create(None, comm=self.comm)
#         self.readMem = np.zeros(2, dtype=np.int64)
#         self.maxInt = np.iinfo(np.int64).max
#         self.writeMem = np.asarray([self.maxInt, self.maxInt], dtype=np.int64)

#     def run(self, startTime):
#         pass
#         timeNow = startTime
#         while True:
#             self.win.Fence()
#             self.win.Get(self.readMem, 0)
#             self.win.Fence()
#             gblTMin, gblTMax = self.getRange()
#             if self.comm.rank == 0 and gblTMax - gblTMin > 1:
#                 print '### Rank 0: global range diverged to (%d, %d) at time %d' % \
#                     (gblTMin, gblTMax, timeNow)
#             if gblTMax > 1000:
#                 print 'rank %d: run complete' % self.comm.rank
#                 MPI.Finalize()
#                 sys.exit(0)
#             if self.comm.rank == 0:
#                 self.writeMem = np.asarray([self.maxInt, self.maxInt])
#                 self.win.Put(self.writeMem, 0)
#             self.win.Fence()
#             timeNow = self.sleep(0)
#             self.win.Accumulate(self.writeMem, 0, op=MPI.MIN)
#
#     def getRange(self):
#         return (np.amin(self.patch.vtime), np.amax(self.patch.vtime))
#
#     def setMyRange(self, myMin, myMax):
#         pass


class GateAgent(Agent):
    def __init__(self, ownerPatch):
        Agent.__init__(self, 'GateAgent', ownerPatch)
        self.cycleCounter = 0
        self.patch = ownerPatch
        self.timeless = True
        self.gateList = []

    def addGate(self, gate):
        gate.lock(self)  # so everyone else is forced into the queue
        self.gateList.append(gate)

    def run(self, startTime):
        timeNow = startTime
        assert timeNow is not None, "Timenow is None"
        while True:
            for gate in self.gateList:
                gate.cycleStart(timeNow)
            timeNow = self.sleep(0)
            assert timeNow is not None, "Timenow is None"
            for gate in self.gateList:
                gate.cycleFinish(timeNow)
            self.cycleCounter += 1


class GateEntrance(Interactant):
    queueBlockSize = 40  # limits network packet size

    def __init__(self, name, ownerPatch, destTag, debug=False):
        Interactant.__init__(self, name, ownerPatch, debug=debug)
        self.destTag = destTag
        self.nInTransit = 0
        self._oldLockQueue = []
        self.logger = logging.getLogger(__name__ + '.GateEntrance')

    def cycleStart(self, timeNow):
        self.logger.debug('%s begins cycleStart; destTag is %s' % (self._name, self.destTag))
        self.nInTransit = len([a for a in self._lockQueue if not a.timeless])
        if self._lockQueue:
            q = self._lockQueue[:]
            while q:
                self.patch.group.enqueue(MsgTypes.GATE, (timeNow, q[:GateEntrance.queueBlockSize]),
                                         self.patch.gblAddr, self.destTag)
                q = q[GateEntrance.queueBlockSize:]
        else:
            self.patch.group.enqueue(MsgTypes.GATE, (timeNow, []),
                                     self.patch.gblAddr, self.destTag)
        self._oldLockQueue = self._lockQueue
        self._lockQueue = []
        self._nEnqueued = 0
        self.logger.debug('%s ends cycleStart' % self._name)

    def cycleFinish(self, timeNow):
        if self._debug:
            self.logger.debug('%s begins cycleFinish' % self._name)
        self.nInTransit = 0
        if not self.patch.group.isLocal(self.destTag):
            for a in self._oldLockQueue:
                a.kill()
        self._oldLockQueue = []
        if self._debug:
            self.logger.debug('%s ends cycleFinish' % self._name)

    def getNWaiting(self):
        return self._nEnqueued + self.nInTransit

    def getWaitingDetails(self):
        """Returns a dict of typeName:nOfThisType entries"""
        result = Interactant.getWaitingDetails(self)
        result['_in_transit'] = self.nInTransit
        return result

    def lock(self, lockingAgent):
        if self._lockingAgent is not None:
            #  This will get enqueued for sending
            if lockingAgent.debug:
                self.logger.debug('%s bound  to gate %s' % (lockingAgent.name, self._name))
        agent.Interactant.lock(self, lockingAgent, debug=False)


class GateExit(Interactant):
    def __init__(self, name, ownerPatch, srcTag, debug=False):
        Interactant.__init__(self, name, ownerPatch, debug=debug)
        self.srcTag = srcTag
        self.partnerMaxDay = 0
        self.logger = logging.getLogger(__name__ + '.GateExit')

    def cycleStart(self, timeNow):
        self.patch.group.expect(self.srcTag, self.patch.gblAddr, self.handleIncoming)

    def cycleFinish(self, timeNow):
        pass

    def handleIncoming(self, msgType, incomingTuple):
        """ This is called by the messaging system to deliver incoming agents """
        if msgType == MsgTypes.GATE:
            senderTime, agentList = incomingTuple
            if self._debug:
                d = {}
                for k in [a.__class__.__name__ for a in agentList]:
                    if k not in d:
                        d[k] = 0
                    d[k] += 1
                self.logger.debug('%s: got time=%s %s agents: %s' % (self._name, senderTime,
                                                                     self.patch.group.nI.vclock,
                                                                     d))
            timeNow = self._ownerLoop.sequencer.getTimeNow()
            for a in agentList:
                a.reHome(self.patch)
                if a.timeless:
                    # Timeless agents always go the the 'current' time
                    if timeNow > senderTime:
                        self.logger.debug('%s jumps forward %d -> %d' %
                                          (a.name, senderTime, timeNow))
                    elif timeNow < senderTime:
                        self.logger.debug('%s jumps backward %d -> %d' %
                                          (a.name, senderTime, timeNow))
                    self._ownerLoop.sequencer.enqueue(a, timeNow)
                else:
                    if timeNow > senderTime:
                        self.logger.critical('%s: MESSAGE FROM THE PAST: %s' %
                                             (self._name, a.name))
                        self.patch.group.nI.comm.Abort()
                    self._ownerLoop.sequencer.enqueue(a, senderTime)
                if a.debug:
                    self.logger.debug('%s materializes at %s' % (a.name, self._name))
        else:
            raise RuntimeError('Unknown message type %s arrived at Gate %s' %
                               (msgType, self._name))


class DateChangeQueue(Interactant):
    pass


class DateChangeMsg(Agent):
    STATE_OUTGOING = 0
    STATE_HOMEWARD = 1
    STATE_TERMINATE = 2

    def __init__(self, name, patch, homeQueueAddr, destQueueAddr, creationVTime,
                 creationDate, debug=True):
        Agent.__init__(self, name, patch, debug=debug)
        self.homeQueueAddr = homeQueueAddr
        self.destQueueAddr = destQueueAddr
        self.creationVTime = creationVTime
        self.creationDate = creationDate
        self.fsmstate = DateChangeMsg.STATE_OUTGOING
        self.timeless = True

    def run(self, startTime):
        timeNow = startTime  # @UnusedVariable
        while True:
            if self.fsmstate == DateChangeMsg.STATE_OUTGOING:
                addr, final = self.patch.getPathTo(self.destQueueAddr)
                if final:
                    self.fsmstate = DateChangeMsg.STATE_HOMEWARD
                timeNow = addr.lock(self)  # @UnusedVariable
            elif self.fsmstate == DateChangeMsg.STATE_HOMEWARD:
                addr, final = self.patch.getPathTo(self.homeQueueAddr)
                if final:
                    self.fsmstate = DateChangeMsg.STATE_TERMINATE
                timeNow = addr.lock(self)  # @UnusedVariable
            elif self.fsmstate == DateChangeMsg.STATE_TERMINATE:
                break

    def __getstate__(self):
        d = Agent.__getstate__(self)
        d['homeQueueAddr'] = self.homeQueueAddr
        d['destQueueAddr'] = self.destQueueAddr
        d['creationVTime'] = self.creationVTime.vec
        d['creationVTimeRank'] = self.creationVTime.rank
        d['creationDate'] = self.creationDate
        d['fsmstate'] = self.fsmstate
        return d

    def __setstate__(self, stateDict):
        Agent.__setstate__(self, stateDict)
        self.homeQueueAddr = stateDict['homeQueueAddr']
        self.destQueueAddr = stateDict['destQueueAddr']
        vec = stateDict['creationVTime']
        self.creationVTime = netinterface.VectorClock(vec.shape[0], stateDict['creationVTimeRank'],
                                                      vec=vec)
        self.creationDate = stateDict['creationDate']
        self.fsmstate = stateDict['fsmstate']


class DateChangeAgent(Agent):
    def __init__(self, name, patch):
        Agent.__init__(self, name, patch)
        self.timeless = True
        self.inputQueue = DateChangeQueue(name+'_rQ', patch)
        self.inputQueue.lock(self)
        self.counter = 0
        self.mostRecentBusyVTime = None
        self.msgDict = {}
        self.logger = logging.getLogger(__name__ + '.DateChangeAgent')

    def run(self, startTime):
        timeNow = startTime
        logDebug = self.logger.isEnabledFor(logging.DEBUG)
        while True:
            dWT = self.patch.doneWithToday()
            bumpTime = False
            if logDebug:
                self.logger.debug('%s: dWT is %s; %d in queue' %
                                  (self.name, dWT, len(self.inputQueue._lockQueue)))
            if dWT:
                nInGroup = 0
                for nm, destAddr in self.patch.serviceLookup('DateChangeQueue'):
                    if not self.patch.isLocal(destAddr):
                        newMsg = DateChangeMsg('%s_msg_%s_%d' % (self.name, nm, self.counter),
                                               self.patch, self.inputQueue.getGblAddr(),
                                               destAddr, self.patch.group.nI.vclock.copy(),
                                               timeNow)
                        self.counter += 1
                        nInGroup += 1
                        self.patch.launch(newMsg, timeNow)
                if nInGroup > 0:
                    tickNow = self.patch.group.nI.vclock.vec[self.patch.group.nI.comm.rank]
                    self.msgDict[tickNow] = (nInGroup, 0)
                    if logDebug:
                        self.logger.debug('%s spawned and launched %d' % (self.name, nInGroup))
                else:
                    bumpTime = True  # we are alone, so just change date
            else:
                self.mostRecentBusyVTime = self.patch.group.nI.vclock.copy()
            while self.inputQueue._lockQueue:
                msg = self.inputQueue._lockQueue[0]
                if logDebug:
                    self.logger.debug('%s found %s in input queue' % (self.name, msg.name))
                if isinstance(msg, DateChangeMsg):
                    if msg.homeQueueAddr == self.inputQueue.getGblAddr():
                        if logDebug:
                            self.logger.debug('DateChangeMsg %s made it home!'
                                              ' Started %s %s, now %s %s'
                                              % (msg.name, msg.creationVTime, msg.creationDate,
                                                 self.patch.group.nI.vclock, timeNow))
                        msgTick = msg.creationVTime.vec[self.patch.group.nI.comm.rank]
                        if logDebug:
                            self.logger.debug('%s: dWT = %s, msgDict %s' %
                                              (self.name, dWT, self.msgDict))
                        if dWT and msgTick in self.msgDict:
                            nSent, nSeen = self.msgDict[msgTick]
                            if logDebug:
                                self.logger.debug('%s considering date change' % self.name)
                            if (msg.creationVTime.before(self.patch.group.nI.vclock)
                                    and msg.creationDate == timeNow
                                    and self.mostRecentBusyVTime.before(msg.creationVTime)):
                                nSeen += 1
                                if logDebug:
                                    self.logger.debug('%s counted date change '
                                                      'msg %d of %d for tick %s' %
                                                      (self.name, nSeen, nSent, msgTick))
                                if nSeen >= nSent:
                                    if logDebug:
                                        self.logger.debug('%s will bump time' % self.name)
                                    bumpTime = True
                                    self.msgDict = {}
                                else:
                                    self.msgDict[msgTick] = (nSent, nSeen)
                            else:
                                if logDebug:
                                    self.logger.debug('%s: change rejected: %s %s %s' %
                                                      (self.name,
                                                       msg.creationVTime.before(self.patch.group
                                                                                .nI.vclock),
                                                       msg.creationDate == timeNow,
                                                       self.mostRecentBusyVTime
                                                       .before(msg.creationVTime)))
                                del self.msgDict[msgTick]
                            msg.fsmstate = DateChangeMsg.STATE_TERMINATE
                            self.inputQueue.awaken(msg)
                        else:
                            if logDebug:
                                self.logger.debug('%s missed chance for date change' % self.name)
                            msg.fsmstate = DateChangeMsg.STATE_TERMINATE
                            self.inputQueue.awaken(msg)
                    else:
                        if ((dWT and msg.creationDate == timeNow)
                                or msg.creationDate < timeNow):
                            if logDebug:
                                self.logger.debug('%s Passing %s forward' % (self.name, msg.name))
                            msg.fsmstate = DateChangeMsg.STATE_HOMEWARD
                            self.inputQueue.awaken(msg)
                        else:
                            if logDebug:
                                self.logger.debug('%s canceling %s: %s %d %d' %
                                                  (self.name, msg.name, dWT, msg.creationDate,
                                                   timeNow))
                            msg.fsmstate = DateChangeMsg.STATE_TERMINATE
                            self.inputQueue.awaken(msg)
                else:
                    raise RuntimeError("%s unexpectedly got the message %s" %
                                       (self.name, msg.name))
            if bumpTime:
                self.logger.info('%s BUMPING TIME' % self.name)
                self.patch.loop.sequencer.bumpTime()
            timeNow = self.sleep(0)


class Patch(object):
    counter = 0

    def _createPerTickCB(self):
        def tickFun(thisAgent, timeLastTick, timeNow):
            # Force the current patch to exit so that the next patch gets
            # a time slice.
            thisAgent.ownerLoop.sequencer.enqueue(thisAgent, timeNow)
            self.group.switch(timeNow)
        return tickFun

    def __init__(self, group, name=None, patchId=None):
        if patchId is None:
            self.patchId = Patch.counter
            Patch.counter += 1
        else:
            self.patchId = patchId
        self.group = group
        self.gblAddr = group.getGblAddr(self.patchId)
        if name is None:
            self.name = "Patch_%s" % str(self.gblAddr)
        else:
            self.name = name
        self.logger = logging.getLogger(__name__ + '.Patch')
        self.loop = agent.MainLoop(self.name + '.loop')
        self.gateAgent = GateAgent(self)
        self.dateChangeAgent = DateChangeAgent(self.name + '_DateChangeAgent', self)
        self.outgoingGateDict = {}
        self.incomingGateDict = {}
        self.interactantDict = {}  # Does not include gates
        self.loop.addPerTickCallback(self._createPerTickCB())
        self.addAgents([self.gateAgent, self.dateChangeAgent])
        self.addInteractants([self.dateChangeAgent.inputQueue])

    def doneWithToday(self):
        """
        If the only agents still in today's sequencer loop are marked 'timeless', this
        will return True; otherwise False.  Note that this condition can be reversed if
        a new agent is inserted into today's loop.
        """
        for iact in (self.outgoingGateDict.values() + self.incomingGateDict.values()
                     + self.interactantDict.values()):
            if (iact._lockingAgent is not None and iact._lockingAgent.timeless
                    and iact.getNWaiting()):
                self.logger.debug('doneWithToday is false because %s has %d waiting: %s' %
                                  (iact, iact.getNWaiting(), iact.getWaitingDetails()))
                return False
        # Gates are not in interactantDict so we have not checked them, but they are timeless
        # by definition.
        return all([a.timeless for a in
                    self.loop.sequencer._timeQueues[self.loop.sequencer._timeNow]])

    def addGateFrom(self, otherPatchTag):
        gateExit = GateExit(("%s.GateExit_%s" % (self.name, otherPatchTag)),
                            self, otherPatchTag, debug=False)
        self.gateAgent.addGate(gateExit)
        self.incomingGateDict[otherPatchTag] = gateExit
        return gateExit

    def addGateTo(self, otherPatchTag):
        gateEntrance = GateEntrance(("%s.GateEntrance_%s" % (self.name, otherPatchTag)),
                                    self, otherPatchTag, debug=False)
        self.gateAgent.addGate(gateEntrance)
        self.outgoingGateDict[otherPatchTag] = gateEntrance
        return gateEntrance

    def addAgents(self, agentList):
        self.loop.addAgents(agentList)

    def addInteractants(self, interactantList):
        for iact in interactantList:
            if isinstance(iact, GateEntrance):
                if iact.destTag not in self.outgoingGateDict:
                    self.outgoingGateDict[iact.destTag] = iact
            elif isinstance(iact, GateExit):
                if iact.srcTag not in self.incomingGateDict:
                    self.incomingGateDict[iact.srcTag] = iact
            else:
                self.interactantDict[iact.getGblAddr()] = iact

    def getAllLocalInteractants(self):
        """
        Returns a generator providing links to everything added with self.addInteractants,
        except for Gates.
        """
        return self.interactantDict.itervalues()

    def __str__(self):
        return '<%s>' % self.name

    def launch(self, agent, startTime):
        agent.parent = self.loop
        self.loop.sequencer.enqueue(agent, startTime)

    def serviceLookup(self, typeNameStr, patchAddr=None):
        if patchAddr is None:
            return self.group.worldInteractants[typeNameStr][:]
        else:
            return [(nm, addr)
                    for nm, addr in self.group.worldInteractants[typeNameStr]
                    if addr.getPatchAddr() == patchAddr]

    def isLocal(self, gblAddr):
        """Is the address local to this patch?"""
        return (netinterface.GblAddr.tupleGetPatchAddr(gblAddr) == self.gblAddr)

    def getPathTo(self, gblAddr):
        if gblAddr in self.interactantDict:
            return (self.interactantDict[gblAddr], True)
        else:
            assert not self.isLocal(gblAddr), \
                "%s: Unknown supposedly-local address %s" % (self.name, gblAddr)
            patchAddr = netinterface.GblAddr.tupleGetPatchAddr(gblAddr)
            if patchAddr in self.outgoingGateDict:
                return (self.outgoingGateDict[patchAddr], False)
            else:
                raise RuntimeError("%s: No path to correct patch for address %s" % (self.name,
                                                                                    gblAddr))


def greenletTrace(event, args):
    if event == 'switch':
        origin, target = args
        # Handle a switch from origin to target.
        # Note that callback is running in the context of target
        # greenlet and any exceptions will be passed as if
        # target.throw() was used instead of a switch.
        print 'TRACE switch %s -> %s (parent %s)' % (origin, target, target.parent)
        return
    if event == 'throw':
        origin, target = args
        # Handle a throw from origin to target.
        # Note that callback is running in the context of target
        # greenlet and any exceptions will replace the original, as
        # if target.throw() was used with the replacing exception.
        print 'TRACE throw %s -> %s' % (origin, target)
        return


def output_func(o):
    import types
    if hasattr(o, 'name'):
        return '%s:%s' % (str(type(o)), o.name)
    elif isinstance(o, types.MethodType):
        return 'instance method %s.%s' % (output_func(o.im_self), o.im_func.__name__)
    else:
        return str(type(o))


class PatchGroup(greenlet):

    def createPerEventCallback(self):
        def evtFun(mainLoop, scalarTimeNow):
            # self.nI.vclock.incr()
            pass
        return evtFun

    def __init__(self, comm, name=None, sync=True, trace=False, deterministic=False,
                 printCensus=False):
        if trace:
            greenlet.settrace(greenletTrace)

        self.patches = []
        self.nI = netinterface.NetworkInterface(comm, sync=sync, deterministic=deterministic)
        if name is None:
            self.name = 'PatchGroup_%d' % comm.rank
        else:
            self.name = name
        self.outgoingDict = {}
        self.outstandingSendReqs = []
        self.outstandingRecvReqs = []
        self.expectFrom = set()
        self.clientGateExits = {}
        self.sync = sync
        self.deterministic = deterministic
        self.printCensus = printCensus
        self.prevTraceCB = None
        self.stopNow = False
        self.logger = logging.getLogger(__name__ + '.PatchGroup')

    def setTrace(self):
        self.prevTraceCB = greenlet.gettrace()
        greenlet.settrace(greenletTrace)

    def unsetTrace(self):
        greenlet.settrace(self.prevTraceCB)

    @property
    def vclock(self):
        print '!!!!!!!!!!!!!!!! is vclock needed globally?'
        return self.nI.vclock

    @property
    def vtime(self):
        print '!!!!!!!!!!!!!!!! is vtime needed globally?'
        return self.nI.vclock.vec

    def barrier(self):
        """
        It is probably never necessary for a user-level program to call barrier() explicitly.
        """
        self.nI.barrier()

    def getGblAddr(self, lclId):
        return self.nI.getGblAddr(lclId)

    def addPatch(self, patch):
        patch.loop.parent = self
        self.patches.append(patch)
        patch.loop.freezeDate()  # No new days until I say so
        patch.loop.addPerEventCallback(self.createPerEventCallback())
        return patch

    def run(self):
        # tr = tracker.SummaryTracker()
        logDebug = self.logger.isEnabledFor(logging.DEBUG)
        while True:
            if logDebug:
                self.logger.debug('%s: new pass of run' % (self.name))
            for p in self.patches[:]:
                if logDebug:
                    self.logger.debug('%s: running patch %s: %s agents at time %s' %
                                      (self.name, p.name, p.loop.sequencer.getNWaitingNow(),
                                       p.loop.sequencer.getTimeNow()))
                reply = p.loop.switch()  # @UnusedVariable
                if self.printCensus:
                    p.loop.printCensus(tickNum=self.nI.vclock.vec[self.nI.comm.rank])

            if logDebug:
                self.logger.debug('%s: finish last recv' % self.name)
            self.nI.finishRecv()
            if logDebug:
                self.logger.debug('%s: finish last send' % self.name)
            self.nI.finishSend()
            if self.stopNow:
                self.logger.debug('%s Sending done signal' % self.name)
                if self.nI.sendDoneSignal():
                    self.logger.debug('%s: everyone is done' % self.name)
                    return '%s claims all done' % self.name
            if logDebug:
                self.logger.debug('%s: start recv' % self.name)
            self.nI.startRecv()
            if logDebug:
                self.logger.debug('%s: start send' % self.name)
            self.nI.startSend()
            if logDebug:
                self.logger.debug('%s: finished networking' % self.name)

    def __str__(self):
        return '<%s>' % self.name

    def enqueue(self, msgType, thing, srcTag, destTag):
        self.nI.enqueue(msgType, thing, srcTag, destTag)

    def expect(self, srcAddr, destAddr, handleIncoming):
        """
        The handleIncoming is a callback with the signature

            handleIncoming(msgType, incomingTuple)

        There is no way to drop a rank from the expected source set because one can
        never be sure there is no straggler message from that rank
        """
        self.nI.expect(srcAddr, destAddr, handleIncoming)

    def shareInteractantDirectories(self, patchList):
        myInteractants = defaultdict(list)
        for p in patchList:
            pId = p.patchId
            myInteractants['_'].append(self.getGblAddr(pId))
            for iact in p.interactantDict.values():
                info = iact.getInfo()
                classNm = iact.__class__.__name__
                myInteractants[classNm].append((info, self.getGblAddr((pId, iact.id))))
        gblAllInteractants = defaultdict(list)
        gblAllPatches = []
        if self.deterministic:
            # Be meticulous, so that things get inserted in consistent order
            dList = self.nI.comm.allgather(myInteractants)
            l = [(d['_'], d) for d in dList]
            l.sort()
            l = [b for a, b in l]  # @UnusedVariable
            for d in l:
                gblAllPatches.extend(d['_'])
                subL = [(k, v) for k, v in d.items() if k != '_']
                subL.sort()
                for k, v in subL:
                    gblAllInteractants[k].extend(v)
        else:
            for d in self.nI.comm.allgather(myInteractants):
                for k, v in d.items():
                    if k == '_':
                        gblAllPatches.extend(v)
                    else:
                        gblAllInteractants[k].extend(v)
        return gblAllInteractants, gblAllPatches

    def isLocal(self, gblAddr):
        return self.nI.isLocal(gblAddr)

    def start(self):
        # Collect remote geometry information.  This includes an implicit barrier
        self.worldInteractants, self.allPatches = self.shareInteractantDirectories(self.patches)

        # Build the global gate network
        for localP in self.patches:
            for friend in self.allPatches:
                if (localP.gblAddr != friend):
                    localP.addGateTo(friend)
                    localP.addGateFrom(friend)

        self.stopNow = False
        return self.switch()

    def stop(self):
        self.stopNow = True

    def doneWithToday(self):
        return all([p.doneWithToday() for p in self.patches])
