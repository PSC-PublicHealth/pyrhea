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
from greenlet import greenlet
import agent
import netinterface
#from pympler import tracker

"""
To Do:
-Add rollback chain support
-Implement rollback
-Add day management logic
"""


class MsgTypes():
    GATE = 0
    ENDOFDAY = 1
    ANTI_ENDOFDAY = 2
    SINGLE_ENDOFDAY = 3
    SINGLE_ANTI_ENDOFDAY = 4


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
        agent.Agent.__init__(self, name, patch.loop, debug=False)
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
    queueBlockSize = 4  # limits network packet size

    def __init__(self, name, ownerPatch, destTag, debug=False):
        Interactant.__init__(self, name, ownerPatch, debug=debug)
        self.destTag = destTag
        self.nInTransit = 0
        self._oldLockQueue = []

    def cycleStart(self, timeNow):
        if self._debug:
            print '%s begins cycleStart; destTag is %s' % (self._name, self.destTag)
        self.nInTransit = len(self._lockQueue)
        if self._lockQueue:
            q = self._lockQueue[:]
            while q:
                self.patch.group.enqueue(MsgTypes.GATE, (timeNow, q[:GateEntrance.queueBlockSize]),
                                         self.patch.tag, self.destTag)
                q = q[GateEntrance.queueBlockSize:]
        else:
            self.patch.group.enqueue(MsgTypes.GATE, (timeNow, []),
                                     self.patch.tag, self.destTag)
        self._oldLockQueue = self._lockQueue
        self._lockQueue = []
        self._nEnqueued = 0
        if self._debug:
            print '%s ends cycleStart' % self._name

    def cycleFinish(self, timeNow):
        if self._debug:
            print '%s begins cycleFinish' % self._name
        self.nInTransit = 0
        if not self.patch.group.isLocal(self.destTag):
            for a in self._oldLockQueue:
                a.kill()
        self._oldLockQueue = []
        if self._debug:
            print '%s ends cycleFinish' % self._name

    def getNWaiting(self):
        return self._nEnqueued + self.nInTransit

    def lock(self, lockingAgent):
        if self._lockingAgent is not None:
            #  This will get enqueued for sending
            if lockingAgent.debug:
                print '%s bound  to gate %s' % (lockingAgent.name, self._name)
        agent.Interactant.lock(self, lockingAgent, debug=False)


class GateExit(Interactant):
    def __init__(self, name, ownerPatch, srcTag, debug=False):
        Interactant.__init__(self, name, ownerPatch, debug=debug)
        self.srcTag = srcTag
        self.partnerEndOfDay = False
        self.partnerMaxDay = 0

    def cycleStart(self, timeNow):
        if self._debug:
            print '%s begins cycleStart; tag is %s' % (self._name, self.srcTag)
        self.patch.group.expect(self.srcTag, self.patch.tag, self.handleIncoming)
        if self._debug:
            print '%s ends cycleStart' % self._name

    def cycleFinish(self, timeNow):
        if self._debug:
            print '%s begins cycleFinish' % self._name

        if self._debug:
            print '%s ends cycleFinish' % self._name

    def handleIncoming(self, msgType, incomingTuple):
        """ This is called by the messaging system to deliver incoming agents """
        if msgType == MsgTypes.GATE:
            senderTime, agentList = incomingTuple
            # print '%s: got time=%s' % (self._name, senderTime)
            # print '%s: agentList: %s' % (self._name, [str(a) for a in agentList])
            timeNow = self._ownerLoop.sequencer.getTimeNow()
            if timeNow > senderTime and agentList:
                print '%s: MESSAGE FROM THE PAST' % self._name
            senderTime = timeNow
            for a in agentList:
                a.reHome(self.patch)
                self._ownerLoop.sequencer.enqueue(a, senderTime)
                if a.debug:
                    print '%s materializes at %s' % (a.name, self._name)
        elif msgType == MsgTypes.SINGLE_ENDOFDAY:
            srcAddr, partnerDay = incomingTuple  # @UnusedVariable
            # print ('%s got EOD %s from %s %s when partnerEOD = %s' %
            #        (self._name, partnerDay, self.srcTag, srcAddr, self.partnerEndOfDay))
            self.partnerMaxDay = partnerDay
            self.partnerEndOfDay = True
        elif msgType == MsgTypes.SINGLE_ANTI_ENDOFDAY:
            srcAddr, partnerDay = incomingTuple  # @UnusedVariable
            # print ('%s got Anti-EOD %d from %s when partnerEOD = %s' %
            #         (self._name, partnerDay, self.srcTag, self.partnerEndOfDay))
            self.partnerMaxDay = partnerDay
            self.partnerEndOfDay = False
        else:
            raise RuntimeError('Unknown message type %s arrived at Gate %s' %
                               (msgType, self._name))


class Patch(object):
    counter = 0

    def _createPerTickCB(self):
        def tickFun(thisAgent, timeLastTick, timeNow):
            if self.endOfDay:
                if self.loop.sequencer.doneWithToday():
                    # Still waiting for partners
                    pass
                else:
                    # Some incoming event has knocked us out of end-of-day
                    self.endOfDay = False
                    for g in self.outgoingGates:
                        self.group.sendGateAntiEOD(self.tag, g.destTag, timeNow)
            else:
                if self.loop.sequencer.doneWithToday():
                    # Newly end-of-day
                    self.endOfDay = True
                    # print ('%s sending out endOfDay with timeNow %s at %s' %
                    #        (self.name, timeNow, self.group.nI.vclock.vec))
                    for g in self.outgoingGates:
                        self.group.sendGateEOD(self.tag, g.destTag, timeNow)
                else:
                    # still not done with the day
                    # print ('%s sequencer.doneWithToday is false at timeNow %s at %s' %
                    #        (self.name, timeNow, self.group.nI.vclock.vec))
                    pass
            allPartnersEOD = all([g.partnerEndOfDay for g in self.incomingGates])
            if len(self.incomingGates) > 0:
                minPartnerDay = min([g.partnerMaxDay for g in self.incomingGates])
                if (self.loop.sequencer.doneWithToday() and self.endOfDay
                        and allPartnersEOD and minPartnerDay >= timeNow):
                    if self.stillEndOfDay:
                        # print '%s: bumping day!' % self.name
                        self.loop.sequencer.bumpIfAllTimeless()
                        timeNow += 1
                        self.endOfDay = False
                        # print '%s: stillEndOfDay -> False' % self.name
                        self.stillEndOfDay = False
                        for g in self.outgoingGates:
                            self.group.sendGateAntiEOD(self.tag, g.destTag, timeNow)
                    else:
                        # print '%s: stillEndOfDay -> True' % self.name
                        self.stillEndOfDay = True
                else:
                    # print '%s: stillEndOfDay -> False' % self.name
                    self.stillEndOfDay = False
            else:
                if (self.loop.sequencer.doneWithToday() and self.endOfDay
                        and allPartnersEOD):
                    if self.stillEndOfDay:
                        # print '%s: bumping day!' % self.name
                        self.loop.sequencer.bumpIfAllTimeless()
                        timeNow += 1
                        self.endOfDay = False
                        # print '%s: stillEndOfDay -> False' % self.name
                        self.stillEndOfDay = False
                        for g in self.outgoingGates:
                            self.group.sendGateAntiEOD(self.tag, g.destTag, timeNow)
                    else:
                        # print '%s: stillEndOfDay -> True' % self.name
                        self.stillEndOfDay = True
                else:
                    # print '%s: stillEndOfDay -> False' % self.name
                    self.stillEndOfDay = False
            # And now we force the current patch to exit so that the next patch gets
            # a time slice.
            thisAgent.ownerLoop.sequencer.enqueue(thisAgent, timeNow)
            self.group.switch(timeNow)
        return tickFun

    def __init__(self, group, name=None):
        self.patchId = Patch.counter
        Patch.counter += 1
        self.group = group
        self.tag = group.getGblAddr(self.patchId)
        if name is None:
            self.name = "Patch_%s" % str(self.tag)
        self.loop = agent.MainLoop(self.name + '.loop')
        self.gateAgent = GateAgent(self)
        self.outgoingGates = []
        self.incomingGates = []
        self.interactants = []  # Does not include gates
        self.loop.addPerTickCallback(self._createPerTickCB())
        self.loop.addAgents([self.gateAgent])
        self.endOfDay = False
        self.stillEndOfDay = False

    def addGateFrom(self, otherPatchTag):
        gateExit = GateExit(("%s.GateExit_%s" % (self.name, otherPatchTag)),
                            self, otherPatchTag, debug=False)
        self.gateAgent.addGate(gateExit)
        self.incomingGates.append(gateExit)
        return gateExit

    def addGateTo(self, otherPatchTag):
        gateEntrance = GateEntrance(("%s.GateEntrance_%s" % (self.name, otherPatchTag)),
                                    self, otherPatchTag, debug=False)
        self.gateAgent.addGate(gateEntrance)
        self.outgoingGates.append(gateEntrance)
        return gateEntrance

    def addAgents(self, agentList):
        self.loop.addAgents(agentList)

    def addInteractants(self, interactantList):
        for iact in interactantList:
            if isinstance(iact, GateEntrance):
                if iact not in self.outgoingGates:
                    self.outgoingGates.append(iact)
            elif isinstance(iact, GateExit):
                if iact not in self.incomingGates:
                    self.incomingGates.append(iact)
            else:
                self.interactants.append(iact)

    def __str__(self):
        return '<%s>' % self.name

    def launch(self, agent, startTime):
        agent.parent = self.loop
        self.loop.sequencer.enqueue(agent, startTime)

    def serviceLookup(self, typeNameStr):
        return self.group.worldInteractants[typeNameStr][:]

    def isLocal(self, gblAddr):
        """Is the address local to this patch?"""
        return (gblAddr.getPatchAddr() == self.tag)

    def getPathTo(self, gblAddr):
        if self.isLocal(gblAddr):
            for itr in self.interactants:
                if itr.getGblAddr() == gblAddr:
                    return (itr, True)
            raise RuntimeError("%s: Unknown supposedly-local address %s" % (self.name, gblAddr))
        else:
            patchAddr = gblAddr.getPatchAddr()
            for g in self.outgoingGates:
                if g.destTag == patchAddr:
                    return (g, False)
            raise RuntimeError("%s: No path to right patch for address %s" % (self.name, gblAddr))


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
            self.nI.vclock.incr()
        return evtFun

    def __init__(self, comm, name=None, sync=True, trace=False, deterministic=False):
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
        self.endOfDay = False
        self.prevTraceCB = None
        self.stopNow = False

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
        while True:
            # print '######## %s: new pass of run' % (self.name)
            for p in self.patches[:]:
                # print '######## %s: running patch %s: %s agents at time %s' % \
                #     (self.name, p.name, p.loop.sequencer.getNWaitingNow(),
                #      p.loop.sequencer.getTimeNow())
                reply = p.loop.switch()  # @UnusedVariable
                # p.loop.printCensus()
#                if p.loop.stopNow:
#                     if reply:
#                         print '%s stopped running: %s' % (p.name, reply)
#                     else:
#                         print '%s stopped running' % p.name
#                    self.patches.remove(p)
#                    self.setTrace()

            # print '######### %s: finish last recv' % self.name
            self.nI.finishRecv()
            # print '######### %s: finish last send' % self.name
            self.nI.finishSend()
            # print ('######### %s: checking done; seen %d, unfinished patches = %d' %
            #        (self.name, self.nI.doneSignalsSeen,
            #         len([1 for p in self.patches if p.loop.stopNow])))
            if self.stopNow:
                print '%s Sending done signal' % self.name
                if self.nI.sendDoneSignal():
                    return '%s claims all done' % self.name
            # print '######### %s: start recv' % self.name
            self.nI.startRecv()
            # print '######### %s: start send' % self.name
            self.nI.startSend()
            # print '######### %s: finished networking' % self.name

    def __str__(self):
        return '<%s>' % self.name

    def enqueue(self, msgType, thing, srcTag, destTag):
        self.nI.enqueue(msgType, thing, srcTag, destTag)

    def expect(self, srcAddr, destAddr, handleIncoming):
        """
        The handleIncoming is a callback with the signature

            handleIncoming(incomingTuple)

        There is no way to drop a rank from the expected source set because one can
        never be sure there is no straggler message from that rank
        """
        self.nI.expect(srcAddr, destAddr, handleIncoming)

    def shareInteractantDirectories(self):
        myInteractants = defaultdict(list)
        for p in self.patches:
            pId = p.patchId
            myInteractants['_'].append(self.getGblAddr(pId))
            for iact in p.interactants:
                nm = iact._name
                classNm = iact.__class__.__name__
                myInteractants[classNm].append((nm, self.getGblAddr((pId, iact.id))))
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
        self.worldInteractants, self.allPatches = self.shareInteractantDirectories()

        # Build the global gate network
        for localP in self.patches:
            for friend in self.allPatches:
                if localP.tag != friend:
                    localP.addGateTo(friend)
                    localP.addGateFrom(friend)

        self.stopNow = False
        return self.switch()

    def stop(self):
        self.stopNow = True

    def doneWithToday(self):
        return all([p.loop.sequencer.doneWithToday() for p in self.patches])

    def sendGateEOD(self, srcAddr, destAddr, currentDay):
        self.enqueue(MsgTypes.SINGLE_ENDOFDAY, (srcAddr, currentDay), srcAddr, destAddr)

    def sendGateAntiEOD(self, srcAddr, destAddr, currentDay):
        self.enqueue(MsgTypes.SINGLE_ANTI_ENDOFDAY, (srcAddr, currentDay), srcAddr, destAddr)
