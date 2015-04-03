#! /usr/bin/env python

_rhea_svn_id_ = "$Id$"

from greenlet import greenlet
import math
import agent
import netinterface

"""
To Do:
-Add day management logic
-Add rollback chain support
-Implement rollback
"""


def getCommWorld():
    """Provide easy access to the world to packages that don't want to know about the network"""
    return netinterface.getCommWorld()


class GblAddr(netinterface.GblAddr):
    pass


class Interactant(agent.Interactant):
    def __init__(self, name, patch, debug=False):
        agent.Interactant.__init__(self, name, patch.loop, debug)
        self.patch = patch


class Agent(agent.Agent):
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
    def __init__(self, name, ownerPatch, destTag, debug=False):
        Interactant.__init__(self, name, ownerPatch, debug=debug)
        self.destTag = destTag
        self.nInTransit = 0

    def cycleStart(self, timeNow):
        if self._debug:
            print '%s begins cycleStart; destTag is %s' % (self._name, self.destTag)
        self.patch.group.enqueue((timeNow, self._lockQueue), self.destTag)
        self.nInTransit = len(self._lockQueue)
        self._lockQueue = []
        self._nEnqueued = 0
        if self._debug:
            print '%s ends cycleStart' % self._name

    def cycleFinish(self, timeNow):
        if self._debug:
            print '%s begins cycleFinish' % self._name
        self.nInTransit = 0
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

    def handleIncoming(self, senderTime, agentList):
        """ This is called by the messaging system to deliver incoming agents """
        print '%s: got time=%s' % (self._name, senderTime)
        print '%s: agentList: %s' % (self._name, [str(a) for a in agentList])
        timeNow = self._ownerLoop.sequencer.getTimeNow()
        if timeNow > senderTime:
            print '%s: MESSAGE FROM THE PAST' % self._name
            senderTime = timeNow
        for a in agentList:
            a.reHome(self.patch)
            self._ownerLoop.sequencer.enqueue(a, senderTime)
            if a.debug:
                print '%s materializes at %s' % (a.name, self._name)


class Patch(object):
    counter = 0

    def _createPerTickCB(self):
        def tickFun(thisAgent, timeLastTick, timeNow):
            tMin, tMax = thisAgent.ownerLoop.sequencer.getTimeRange()
            worldTMax = self.group.nI.vclock.max()
            worldTMin = self.group.nI.vclock.min()
            print '%s: tick! time change %s -> %s in range (%s, %s), world range (%s, %s)' % \
                (self.name, timeLastTick, timeNow, tMin, tMax, worldTMin, worldTMax)
            # thisAgent.ownerLoop.printCensus()
            if tMin > worldTMin:
                self.loop.freezeTime()
                print '%s: time frozen since %s > %s' % (self.name, tMin, worldTMin)
            else:
                if thisAgent.ownerLoop.timeFrozen:
                    print '%s: time unfrozen' % self.name
                thisAgent.ownerLoop.unfreezeTime()
            # And now we force the current patch to exit so that the next patch gets
            # a time slice.
            thisAgent.ownerLoop.sequencer.enqueue(thisAgent, timeNow)
            self.group.switch(timeNow)
        return tickFun

    def __init__(self, group, name=None):
        self.patchId = Patch.counter
        Patch.counter += 1
        self.tag = group.getGblAddr(self.patchId)
        if name is None:
            self.name = "Patch_%s" % self.tag
        self.loop = agent.MainLoop(self.name + '.loop')
        self.gateAgent = GateAgent(self)
        self.loop.addPerTickCallback(self._createPerTickCB())
        self.loop.addAgents([self.gateAgent])

    def addGateFrom(self, otherPatchTag):
        gateExit = GateExit(("%s.GateExit_%s" % (self.name, otherPatchTag)),
                            self, otherPatchTag, debug=False)
        self.gateAgent.addGate(gateExit)
        return gateExit

    def addGateTo(self, otherPatchTag):
        gateEntrance = GateEntrance(("%s.GateEntrance_%s" % (self.name, otherPatchTag)),
                                    self, otherPatchTag, debug=False)
        self.gateAgent.addGate(gateEntrance)
        return gateEntrance

    def addAgents(self, agentList):
        self.loop.addAgents(agentList)

    def __str__(self):
        return '<%s>' % self.name


class PatchGroup(greenlet):
    def __init__(self, name, comm, sync=True):
        self.patches = []
        self.comm = comm
        self.nI = netinterface.NetworkInterface(self.comm, sync=sync)
        self.name = name
        self.outgoingDict = {}
        self.outstandingSendReqs = []
        self.outstandingRecvReqs = []
        self.expectFrom = set()
        self.clientGateExits = {}
        self.sync = sync

    @property
    def vclock(self):
        print '!!!!!!!!!!!!!!!! is vclock needed globally?'
        return self.nI.vclock

    @property
    def vtime(self):
        print '!!!!!!!!!!!!!!!! is vtime needed globally?'
        return self.nI.vclock.vec

    def barrier(self):
        self.nI.barrier()

    def getGblAddr(self, lclId):
        return self.nI.getGblAddr(lclId)

    def addPatch(self, patch):
        patch.group = self
        patch.loop.parent = self
        self.patches.append(patch)

    def run(self):
        while True:
            for p in self.patches:
                # print '######## %s: running patch %s' % (self.name, p.name)
                reply = p.loop.switch()  # @UnusedVariable
            # print '######### %s: finish last send' % self.name
            self.nI.finishSend()
            # print '######### %s: finish last recv' % self.name
            self.nI.finishRecv()
            # print '######### %s: start recv' % self.name
            if self.doneWithToday():
                print ('!!!!!!!!!!!!!!!!!!!! %s done with today at %s!' %
                       (self.name, self.nI.vclock.vec))
            self.nI.startRecv()
            # print '######### %s: start send' % self.name
            self.nI.startSend()
            # print '######### %s: finished networking' % self.name

    def __str__(self):
        return '<%s>' % self.name

    def enqueue(self, thing, destTag):
        self.nI.enqueue(thing, destTag)

    def expect(self, srcAddr, destAddr, handleIncoming):
        """
        the handleIncoming is a callback with the signature

            handleIncoming(time, agentList)

        There is no way to drop a rank from the expected source set because one can
        never be sure there is no straggler message from that rank
        """
        self.nI.expect(srcAddr, destAddr, handleIncoming)

    def doneWithToday(self):
        return all([p.loop.sequencer.doneWithToday() for p in self.patches])
