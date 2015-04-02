#! /usr/bin/env python

_rhea_svn_id_ = "$Id$"

import sys
from mpi4py import MPI
from greenlet import greenlet
import numpy as np
import math
import agent
from random import randint


class Interactant(agent.Interactant):
    def __init__(self, name, patch, debug=False):
        agent.Interactant.__init__(self, name, patch.loop, debug)
        self.patch = patch


class Agent(agent.Agent):
    def __init__(self, name, ownerPatch, debug=False):
        agent.Agent.__init__(self, name, ownerPatch.loop, debug=False)
        self.patch = ownerPatch

    def reHome(self, newPatch):
        self.ownerLoop = newPatch.loop
        self.parent = newPatch.loop
        self.patch = newPatch
        print '%s home is now %s' % (self, newPatch)


class VectorClock(object):
    def __init__(self, commSize, rank):
        self.rank = rank
        self.vec = np.zeros(commSize, dtype=np.int32)

    def incr(self):
        self.vec[self.rank] += 1

    def merge(self, foreignVec):
        """ This operation does not include incrememting the local time """
        self.vec = np.maximum(self.vec, foreignVec)


class OmniClock(Agent):
    def __init__(self, ownerPatch):
        Agent.__init__(self, 'OmniClock', ownerPatch)
        self.timeless = True
        self.patch = ownerPatch
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

    def run(self, startTime):
        pass
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

    def getRange(self):
        return (np.amin(self.patch.vtime), np.amax(self.patch.vtime))

    def setMyRange(self, myMin, myMax):
        pass


class GateAgent(Agent):
    def __init__(self, ownerPatch):
        Agent.__init__(self, 'GateAgent', ownerPatch)
        self.cycleCounter = 0
        self.patch = ownerPatch
        self.comm = self.patch.comm
        self.timeless = True
        self.gateList = []

    def addGate(self, gate):
        gate.lock(self)  # so everyone else is forced into the queue
        self.gateList.append(gate)

    def run(self, startTime):
        timeNow = startTime
        assert timeNow is not None, "Timenow is None"
        while True:
            vTimeNow = self.patch.vtime
            for gate in self.gateList:
                gate.cycleStart(timeNow, vTimeNow)
            timeNow = self.sleep(0)
            assert timeNow is not None, "Timenow is None"
            for gate in self.gateList:
                gate.cycleFinish(timeNow)
            self.cycleCounter += 1


class GateEntrance(Interactant):
    def __init__(self, name, ownerPatch, destTag, debug=False):
        Interactant.__init__(self, name, ownerPatch, debug=debug)
        self.destTag = destTag
        self.comm = ownerPatch.comm
        self.nInTransit = 0

    def cycleStart(self, timeNow, vTimeNow):
        if self._debug:
            print '%s begins cycleStart; destTag is %s' % (self._name, self.destTag)
        print '%s: sending vtime %s' % (self._name, vTimeNow)
        self.patch.group.enqueue((timeNow, vTimeNow, self._lockQueue),
                                 self.destTag)
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

    def cycleStart(self, timeNow, vTimeNow):
        if self._debug:
            print '%s begins cycleStart; tag is %s' % (self._name, self.srcTag)
        self.patch.group.expect(self.srcTag, self.patch.tag, self)
        if self._debug:
            print '%s ends cycleStart' % self._name

    def cycleFinish(self, timeNow):
        if self._debug:
            print '%s begins cycleFinish' % self._name

        if self._debug:
            print '%s ends cycleFinish' % self._name

    def handleIncoming(self, senderTime, senderVTime, agentList):
        """ This is called by the messaging system to deliver incoming agents """
        print '%s: got time=%s vtime=%s' % (self._name, senderTime, senderVTime)
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
            self.clock.setMyRange(tMin, tMax)
            self.group.vclock.incr()
            worldTMin, worldTMax = self.clock.getRange()  # @UnusedVariable
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
            print '%s: vtime is %s' % (self.name, self.vtime)
            # And now we force the current patch to exit so that the next patch gets
            # a time slice.
            thisAgent.ownerLoop.sequencer.enqueue(thisAgent, timeNow)
            self.group.switch(timeNow)
        return tickFun

    @staticmethod
    def patchTag(patchRank, patchId):
        return (100000*patchId) + patchRank

    @staticmethod
    def parseTag(patchTag):
        patchRank = patchTag % 100000
        patchId = (patchTag - patchRank) / 100000
        return patchRank, patchId

    def __init__(self, comm, name=None):
        self.comm = comm
        assert math.log10(comm.size) < 5.0, "Gate tags were built assuming rank < 100000"
        self.patchId = Patch.counter
        Patch.counter += 1
        self.creationRank = comm.rank
        if name is None:
            self.name = "Patch_%d_%d" % (self.creationRank, self.patchId)
        self.tag = self.patchTag(self.comm.rank, self.patchId)
        self.loop = agent.MainLoop(self.name + '.loop')
        self.clock = OmniClock(self)
        self.gateAgent = GateAgent(self)
        self.loop.addPerTickCallback(self._createPerTickCB())
        # self.loop.addAgents([self.clock, self.gateAgent])
        self.loop.addAgents([self.gateAgent])

    @property
    def vtime(self):
        return self.group.vclock.vec

    def addGateFrom(self, otherPatchRank, otherPatchId):
        gateExit = GateExit(("%s.GateExit_%d_%d" % (self.name, otherPatchRank, otherPatchId)),
                            self,
                            Patch.patchTag(otherPatchRank, otherPatchId),
                            debug=False)
        self.gateAgent.addGate(gateExit)
        return gateExit

    def addGateTo(self, otherPatchRank, otherPatchId):
        gateEntrance = GateEntrance(("%s.GateEntrance_%d_%d" % (self.name, otherPatchRank,
                                                                otherPatchId)),
                                    self,
                                    Patch.patchTag(otherPatchRank, otherPatchId),
                                    debug=False)
        self.gateAgent.addGate(gateEntrance)
        return gateEntrance

    def addAgents(self, agentList):
        self.loop.addAgents(agentList)

    def __str__(self):
        return '<%s>' % self.name


class PatchGroup(greenlet):
    def __init__(self, name, comm):
        self.patches = []
        self.comm = comm
        self.vclock = VectorClock(self.comm.size, self.comm.rank)
        self.name = name
        self.outgoingDict = {}
        self.outstandingSendReqs = []
        self.outstandingRecvReqs = []
        self.expectFrom = set()
        self.clientGateExits = {}

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
            self.finishSend()
            # print '######### %s: finish last recv' % self.name
            self.finishRecv()
            # print '######### %s: start recv' % self.name
            self.startRecv()
            # print '######### %s: start send' % self.name
            self.startSend()
            # print '######### %s: finished networking' % self.name

    def __str__(self):
        return '<%s>' % self.name

    def enqueue(self, thing, destTag):
        destRank, destId = Patch.parseTag(destTag)  # @UnusedVariable
        if destRank not in self.outgoingDict:
            self.outgoingDict[destRank] = []
        self.outgoingDict[destRank].append((destTag, thing))

    def expect(self, srcTag, destTag, gateExit):
        """
        There is no way to drop a rank from the expected source set because one can
        never be sure there is no straggler message from that rank
        """
        srcRank, srcTag = Patch.parseTag(srcTag)  # @UnusedVariable
        self.expectFrom.add(srcRank)
        self.clientGateExits[destTag] = gateExit

    def startRecv(self):
        for srcRank in self.expectFrom:
            self.outstandingRecvReqs.append(self.comm.irecv(None, srcRank, MPI.ANY_TAG))

    def finishRecv(self):
        while True:
            if not self.outstandingRecvReqs:
                # print '######## %s nothing to recv' % self.name
                break
            # print ('######## %s entering recv testany on %d outstanding reqs' %
            #        (self.name, len(self.outstandingRecvReqs)))
            idx, flag, msg = MPI.Request.testany(self.outstandingRecvReqs)
            if idx >= 0:
                self.outstandingRecvReqs.pop(idx)
            if flag:
                for destTag, partTpl in msg:
                    tm, vtm, agentList = partTpl
                    # print ('######## %s: got %d agents for %s, time=%s, vtime=%s' %
                    #        (self.name, len(agentList), destTag, tm, vtm))
                    self.vclock.merge(vtm)
                    if agentList:
                        self.clientGateExits[destTag].handleIncoming(tm, vtm, agentList)
            else:
                # print '######## %s empty recv queue' % self.name
                break
        self.outstandingRecvReqs = []

    def startSend(self):
        for destRank, msgList in self.outgoingDict.items():
            bigCargo = []
            for destTag, cargo in msgList:
                bigCargo.append((destTag, cargo))
            self.outstandingSendReqs.append(self.comm.isend(bigCargo, destRank))
            # print '######### %s sent %s to %s' % (self.name, len(bigCargo), destRank)
            del self.outgoingDict[destRank]

    def finishSend(self):
        # print ('######## %s entering send waitall on %d requests' %
        #        (self.name, len(self.outstandingSendReqs)))
        result = MPI.Request.waitall(self.outstandingSendReqs)  # @UnusedVariable
        # print '######## %s finished send waitall; result was %s' % (self.name, result)
        self.outstandingSendReqs = []


class TestAgent(Agent):
    def __init__(self, name, ownerPatch, debug=False):
        Agent.__init__(self, name, ownerPatch, debug)

    def run(self, startTime):
        timeNow = startTime
        while True:
            fate = randint(0, len(self.patch.interactants))
            if fate == len(self.patch.interactants):
                whichGate = randint(0, len(self.patch.gates)-1)
                gate = self.patch.gates[whichGate]
                if self.debug:
                    print '%s is jumping to %s at %s' % (self.name, whichGate, timeNow)
                timeNow = gate.lock(self)
                timeNow = gate.unlock(self)  # but it's no longer going to be there
            else:
                timeNow = self.patch.interactants[fate].lock(self)
                if self.debug:
                    print 'progress for %s' % self.name
                timeNow = self.sleep(1)
                timeNow = self.patch.interactants[fate].unlock(self)
        if self.debug:
            return '%s is exiting at %s' % (self, timeNow)

    def __getstate__(self):
        d = agent.Agent.__getstate__(self)
        return d

    def __setstate__(self, stateDict):
        agent.Agent.__setstate__(self, stateDict)


class TestPatch(Patch):
    def __init__(self, comm, name=None):
        Patch.__init__(self, comm, name)
        self.interactants = []
        self.gates = []

    def setInteractants(self, interactantList):
        self.interactants = interactantList

    def addGate(self, gate):
        self.gates.append(gate)


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


def describeSelf():
    print """
    This main provides diagnostics. -t, -v and -d for trace, verbose and debug respectively.
    """


def main():
    trace = False
    verbose = False
    debug = False

    for a in sys.argv[1:]:
        if a == '-v':
            verbose = True
        elif a == '-d':
            debug = True
        elif a == '-t':
            trace = True
        else:
            describeSelf()
            sys.exit('unrecognized argument %s' % a)

    comm = MPI.COMM_WORLD
    rank = comm.rank
    print 'Hello from %s' % rank
    if trace:
        greenlet.settrace(greenletTrace)

    patchGroup = PatchGroup('PatchGroup_%d' % rank, comm)
    nPatches = 2
#     if rank == 0:
#         nPatches = 2
    for j in xrange(nPatches):
        print 'rank %d point 0' % rank

        patch = TestPatch(comm)
        print 'rank %d point 0b' % rank

        patch.setInteractants([Interactant('%s_%d_%d' % (nm, rank, j), patch)
                               for nm in ['SubA', 'SubB', 'SubC']])

        leftGateEntrance = patch.addGateTo((rank-1) % comm.size,
                                           (j+1) % nPatches)
        patch.addGate(leftGateEntrance)

        leftGateExit = patch.addGateFrom((rank-1) % comm.size,  # @UnusedVariable
                                         (j-1) % nPatches)

        rightGateEntrance = patch.addGateTo((rank+1) % comm.size,
                                            (j-1) % nPatches)
        patch.addGate(rightGateEntrance)

        rightGateExit = patch.addGateFrom((rank+1) % comm.size,  # @UnusedVariable
                                          (j-1) % nPatches)

        allAgents = []
        for i in xrange(100):
            allAgents.append(TestAgent('Agent_%d_%d_%d' % (rank, j, i),
                                       patch, debug=debug))

        patch.addAgents(allAgents)
        patchGroup.addPatch(patch)
    print 'Rank %s reaches barrier' % rank
    comm.Barrier()
    print 'Rank %s leaves barrier' % rank
    patchGroup.switch()
    print '%d all done (from main)' % rank


############
# Main hook
############

if __name__ == "__main__":
    main()
