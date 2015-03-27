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


class Agent(agent.Agent):
    def __init__(self, name, ownerPatch, debug=False):
        agent.Agent.__init__(self, name, ownerPatch.loop, debug=False)
        self.patch = ownerPatch

    def reHome(self, newPatch):
        self.ownerLoop = newPatch.loop
        self.parent = newPatch.loop
        self.patch = newPatch
        print '%s home is now %s' % (self, newPatch)


class OmniClock(Agent):
    def __init__(self, ownerPatch):
        agent.Agent.__init__(self, 'OmniClock', ownerPatch.loop)
        self.timeless = True
        self.patch = ownerPatch
        self.comm = ownerPatch.comm
        self.minVal = None
        self.maxVal = None
        if self.comm.rank == 0:
            self.mem = np.zeros(2, dtype=np.int64)
            self.win = MPI.Win.Create(self.mem, disp_unit=1, comm=self.comm)
        else:
            self.mem = None
            self.win = MPI.Win.Create(None, comm=self.comm)
        self.readMem = np.zeros(2, dtype=np.int64)
        self.maxInt = np.iinfo(np.int64).max
        self.writeMem = np.asarray([self.maxInt, self.maxInt], dtype=np.int64)

    def run(self, startTime):
        timeNow = startTime
        while True:
            self.win.Fence()
            self.win.Get(self.readMem, 0)
            self.win.Fence()
            gblTMin, gblTMax = self.getRange()
            if self.comm.rank == 0 and gblTMax - gblTMin > 1:
                print '### Rank 0: global range diverged to (%d, %d) at time %d' % \
                    (gblTMin, gblTMax, timeNow)
            if gblTMax > 1000:
                print 'rank %d: run complete' % self.comm.rank
                MPI.Finalize()
                sys.exit(0)
            if self.comm.rank == 0:
                self.writeMem = np.asarray([self.maxInt, self.maxInt])
                self.win.Put(self.writeMem, 0)
            self.win.Fence()
            timeNow = self.sleep(0)
            self.win.Accumulate(self.writeMem, 0, op=MPI.MIN)

    def getRange(self):
        return (self.readMem[0], -self.readMem[1])

    def setMyRange(self, myMin, myMax):
        self.writeMem[0] = myMin
        self.writeMem[1] = -myMax


class GateAgent(Agent):
    def __init__(self, ownerPatch):
        Agent.__init__(self, 'GateAgent', ownerPatch)
        self.cycleCounter = 0
        self.patch = ownerPatch
        self.comm = self.patch.comm
        self.mem = np.zeros(self.comm.size, dtype=np.int64)
        self.win = MPI.Win.Create(self.mem, disp_unit=4, comm=self.comm)
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
    def __init__(self, name, toRank, ownerPatch, debug=False, srcTag=None):
        Interactant.__init__(self, name, ownerPatch, debug=debug)
        if srcTag is None:
            self.srcTag = self.comm.rank
        else:
            self.srcTag = srcTag
        self.toRank = toRank
        self.comm = ownerPatch.comm
        self.sendRequest = None
        self.nInTransit = 0

    def cycleStart(self, timeNow):
        if self._debug:
            print '%s begins cycleStart; tag is %s' % (self._name, self.srcTag)
        print '%s outgoing agentList: %s' % (self._name, [str(a) for a in self._lockQueue])
        self.sendRequest = self.comm.isend([timeNow, self._lockQueue],
                                           self.toRank, self.srcTag)
        self.nInTransit = len(self._lockQueue)
        self._lockQueue = []
        self._nEnqueued = 0
        if self._debug:
            print '%s ends cycleStart' % self._name

    def cycleFinish(self, timeNow):
        if self._debug:
            print '%s begins cycleFinish' % self._name
        self.sendRequest.wait()
        self.sendRequest = None
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
    def __init__(self, name, fromRank, ownerPatch, debug=False, srcTag=None):
        Interactant.__init__(self, name, ownerPatch, debug=debug)
        self.fromRank = fromRank
        self.comm = ownerPatch.comm
        if srcTag is None:
            self.srcTag = self.fromRank
        else:
            self.srcTag = srcTag
        self.rcvRequest = None

    def cycleStart(self, timeNow):
        if self._debug:
            print '%s begins cycleStart; tag is %s' % (self._name, self.srcTag)
        self.rcvRequest = self.comm.irecv(None, self.fromRank, tag=self.srcTag)
        if self._debug:
            print '%s ends cycleStart' % self._name

    def cycleFinish(self, timeNow):
        if self._debug:
            print '%s begins cycleFinish' % self._name
        incoming = self.rcvRequest.wait()
        # Format of the received object matches that encoded by GateEntrance.cycleStart()
        senderTime = incoming[0]
        agentList = incoming[1]
        print '%s: agentList: %s' % (self._name, [str(a) for a in agentList])
        if senderTime != timeNow and len(agentList) > 0:
            print "%s has time mismatch, %s vs. %s" % (self._name, senderTime, timeNow)
        if senderTime > timeNow:
            raise RuntimeError('Message from the past!')
        for a in agentList:
            a.reHome(self.patch)
            self._ownerLoop.sequencer.enqueue(a, timeNow)
            if a.debug:
                print '%s materializes at %s' % (a.name, self._name)
        if self._debug:
            print '%s ends cycleFinish' % self._name


def describeSelf():
    print """This main provides diagnostics. -v and -d for verbose and debug respectively."""


class Patch(object):
    counter = 0

    def _createPerTickCB(self):
        def tickFun(thisLoop, timeLastTick, timeNow):
            tMin, tMax = thisLoop.sequencer.getTimeRange()
            self.clock.setMyRange(tMin, tMax)
            worldTMin, worldTMax = self.clock.getRange()  # @UnusedVariable
            print '%s: tick! time change %s -> %s in range (%s, %s), world range (%s, %s)' % \
                (self.name, timeLastTick, timeNow, tMin, tMax, worldTMin, worldTMax)
            # thisLoop.printCensus()
            if tMin > worldTMin:
                self.loop.freezeTime()
                print '%s: time frozen since %s > %s' % (self.name, tMin, worldTMin)
            else:
                if thisLoop.timeFrozen:
                    print '%s: time unfrozen' % self.name
                thisLoop.unfreezeTime()
        return tickFun

    def _patchTag(self, patchRank, patchId):
        return (100000*patchId) + patchRank

    def __init__(self, comm, name=None):
        self.comm = comm
        assert math.log10(comm.size) < 5.0, "Gate tags were built assuming rank < 100000"
        self.patchId = Patch.counter
        Patch.counter += 1
        self.creationRank = comm.rank
        if name is None:
            self.name = "Patch_%d_%d" % (self.creationRank, self.patchId)
        self.loop = agent.MainLoop(self.name + '.loop')
        self.clock = OmniClock(self)
        self.gateAgent = GateAgent(self)
        self.loop.addPerTickCallback(self._createPerTickCB())
        self.loop.addAgents([self.clock, self.gateAgent])

    def addGateFrom(self, otherPatchRank, otherPatchId):
        gateExit = GateExit(("%s.GateExit_%d_%d" % (self.name, otherPatchRank, otherPatchId)),
                            otherPatchRank, self, debug=True,
                            srcTag=self._patchTag(self.comm.rank, self.patchId))
        self.gateAgent.addGate(gateExit)
        return gateExit

    def addGateTo(self, otherPatchRank, otherPatchId):
        gateEntrance = GateEntrance(("%s.GateEntrance_%d_%d" % (self.name, otherPatchRank,
                                                                otherPatchId)),
                                    otherPatchRank, self, debug=True,
                                    srcTag=self._patchTag(otherPatchRank, otherPatchId))
        self.gateAgent.addGate(gateEntrance)
        return gateEntrance

    def addAgents(self, agentList):
        self.loop.addAgents(agentList)

    def __str__(self):
        return '<%s>' % self.name


class PatchGroup(greenlet):
    def __init__(self, name, stepsPerPatch=10):
        self.patches = []
        self.stepsPerPatch = stepsPerPatch
        self.name = name

    def addPatch(self, patch):
        patch.loop.parent = self
        self.patches.append(patch)

    def run(self):
        while True:
            for p in self.patches:
                reply = p.loop.switch(limit=self.stepsPerPatch)

    def __str__(self):
        return '<%s>' % self.name


class TestAgent(Agent):
    def __init__(self, name, ownerPatch, debug=False):
        Agent.__init__(self, name, ownerPatch, debug)
        self.interactants = ownerPatch.interactants
        self.gateEntrance = ownerPatch.gateEntrance

    def run(self, startTime):
        timeNow = startTime
        while True:
            fate = randint(0, len(self.interactants))
            if fate == len(self.interactants):
                if self.debug:
                    print '%s is jumping at %s' % (self.name, timeNow)
                self.gateEntrance.lock(self)
                timeNow = self.sleep(0)  # yield thread
                self.gateEntrance.unlock(self)  # but it's no longer going to be gateEntrance
            else:
                timeNow = self.interactants[fate].lock(self)
                if self.debug:
                    print 'progress for %s' % self.name
                timeNow = self.sleep(1)
                timeNow = self.interactants[fate].unlock(self)
        if self.debug:
            return '%s is exiting at %s' % (self, timeNow)

    def reHome(self, newOwnerPatch):
        agent.Agent.reHome(self, newOwnerPatch)
        self.interactants = newOwnerPatch.interactants
        self.gateEntrance = newOwnerPatch.gateEntrance

    def __getstate__(self):
        d = agent.Agent.__getstate__(self)
        return d

    def __setstate__(self, stateDict):
        agent.Agent.__setstate__(self, stateDict)


class TestPatch(Patch):
    def __init__(self, comm, name=None):
        Patch.__init__(self, comm, name)
        self.interactants = []
        self.gateEntrance = None

    def setInteractants(self, interactantList):
        self.interactants = interactantList

    def setGateEntrance(self, gateEntrance):
        self.gateEntrance = gateEntrance


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


def main():
    global verbose, debug

    for a in sys.argv[1:]:
        if a == '-v':
            verbose = True
        elif a == '-d':
            debug = True
        else:
            describeSelf()
            sys.exit('unrecognized argument %s' % a)

    comm = MPI.COMM_WORLD
    rank = comm.rank
    print 'Hello from %s' % rank
#     if rank == 0:
#         greenlet.settrace(greenletTrace)

    patchGroup = PatchGroup('PatchGroup_%d' % rank)
    for j in xrange(3):

        patch = TestPatch(comm)

        patch.setInteractants([Interactant('%s_%d_%d' % (nm, rank, j), patch)
                               for nm in ['SubA', 'SubB', 'SubC']])

        toRank = (comm.rank + 1) % comm.size
        patch.setGateEntrance(patch.addGateTo(toRank, (j+1) % 3))

        fromRank = (comm.rank + comm.size - 1) % comm.size
        gateExit = patch.addGateFrom(fromRank, (j-1) % 3)  # @UnusedVariable

        allAgents = []
        for i in xrange(100):
            allAgents.append(TestAgent('Agent_%d_%d_%d' % (rank, j, i),
                                       patch, debug=False))

        patch.addAgents(allAgents)
        patchGroup.addPatch(patch)
    patchGroup.switch()
    print '%d all done (from main)' % rank


############
# Main hook
############

if __name__ == "__main__":
    main()
