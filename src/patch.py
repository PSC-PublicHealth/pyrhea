#! /usr/bin/env python

_rhea_svn_id_ = "$Id$"

import sys
from mpi4py import MPI
import numpy as np
import agent
from random import randint


class OmniClock(agent.Agent):
    def __init__(self, comm, ownerLoop):
        agent.Agent.__init__(self, 'OmniClock', ownerLoop)
        self.timeless = True
        self.comm = comm
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


class GateAgent(agent.Agent):
    def __init__(self, comm, ownerLoop):
        agent.Agent.__init__(self, 'GateAgent', ownerLoop)
        self.cycleCounter = 0
        self.comm = comm
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


class GateIn(agent.Interactant):
    def __init__(self, name, toRank, comm, ownerLoop, debug=False):
        agent.Interactant.__init__(self, name, ownerLoop, debug=debug)
        self.toRank = toRank
        self.comm = comm
        self.sendRequest = None
        self.nInTransit = 0

    def cycleStart(self, timeNow):
        if self._debug:
            print '%s begins cycleStart' % self._name
        self.sendRequest = self.comm.isend([timeNow, self._lockQueue],
                                           self.toRank, self.comm.rank)
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


class GateOut(agent.Interactant):
    def __init__(self, name, fromRank, comm, ownerLoop, debug=False):
        agent.Interactant.__init__(self, name, ownerLoop, debug=debug)
        self.fromRank = fromRank
        self.comm = comm
        self.rcvRequest = None

    def cycleStart(self, timeNow):
        if self._debug:
            print '%s begins cycleStart' % self._name
        self.rcvRequest = self.comm.irecv(None, self.fromRank, tag=self.fromRank)
        if self._debug:
            print '%s ends cycleStart' % self._name

    def cycleFinish(self, timeNow):
        if self._debug:
            print '%s begins cycleFinish' % self._name
        incoming = self.rcvRequest.wait()
        # Format of the received object matches that encoded by GateIn.cycleStart()
        senderTime = incoming[0]
        agentList = incoming[1]
        if senderTime != timeNow and len(agentList) > 0:
            print "%s has time mismatch, %s vs. %s" % (self._name, senderTime, timeNow)
        if senderTime > timeNow:
            raise RuntimeError('Message from the past!')
        for a in agentList:
            a.reHome(self._ownerLoop)
            self._ownerLoop.sequencer.enqueue(a, timeNow)
            if a.debug:
                print '%s materializes at %s' % (a.name, self._name)
        if self._debug:
            print '%s ends cycleFinish' % self._name


def describeSelf():
    print """This main provides diagnostics. -v and -d for verbose and debug respectively."""


def createPerTickCB(comm, omniClock):
    def tickFun(myLoop, timeLastTick, timeNow):
        tMin, tMax = myLoop.sequencer.getTimeRange()
        omniClock.setMyRange(tMin, tMax)
        worldTMin, worldTMax = omniClock.getRange()  # @UnusedVariable
        print 'rank %d: tick! time change %s -> %s in range (%s, %s), world range (%s, %s)' % \
            (comm.rank, timeLastTick, timeNow, tMin, tMax, worldTMin, worldTMax)
        # myLoop.printCensus()
        if tMin > worldTMin:
            myLoop.freezeTime()
            print 'rank %d: time frozen since %s > %s' % (comm.rank, tMin, worldTMin)
        else:
            if myLoop.timeFrozen:
                print 'rank %d: time unfrozen' % comm.rank
            myLoop.unfreezeTime()
    return tickFun

gateIn = None

gateOut = None

interactants = None


class TestAgent(agent.Agent):
    def run(self, startTime):
        global gateIn, gateOut, interactants
        timeNow = startTime
        while True:
            fate = randint(0, len(interactants))
            if fate == len(interactants):
                if self.debug:
                    print '%s is jumping at %s' % (self.name, timeNow)
                gateIn.lock(self)
                timeNow = self.sleep(0)  # yield thread
                gateIn.unlock(self)  # but it's no longer going to be gateIn
            else:
                timeNow = interactants[fate].lock(self)
                if self.debug:
                    print 'progress for %s' % self.name
                timeNow = self.sleep(1)
                timeNow = interactants[fate].unlock(self)
        if self.debug:
            return '%s is exiting at %s' % (self, timeNow)

    def __getstate__(self):
        d = agent.Agent.__getstate__(self)
        return d

    def __setstate__(self, stateDict):
        agent.Agent.__setstate__(self, stateDict)


def main():
    global verbose, debug, gateIn, gateOut, interactants

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

    mainLoop = agent.MainLoop('Mainloop_%s' % rank)

    interactants = [agent.Interactant(nm, mainLoop) for nm in ['SubA', 'SubB', 'SubC']]

    allAgents = []
    for i in xrange(100):
        allAgents.append(TestAgent('Agent_%d_%d' % (rank, i), mainLoop))

    omniClock = OmniClock(comm, mainLoop)
    gateAgent = GateAgent(comm, mainLoop)
    toRank = (comm.rank + 1) % comm.size
    gateIn = GateIn(("GateIn_%d_%d" % (comm.rank, toRank)), toRank, comm, mainLoop, debug=False)
    gateAgent.addGate(gateIn)
    fromRank = (comm.rank + comm.size - 1) % comm.size
    gateOut = GateOut(("GateOut_%d_%d" % (fromRank, comm.rank)), fromRank, comm, mainLoop, debug=False)
    gateAgent.addGate(gateOut)
    mainLoop.addPerTickCallback(createPerTickCB(comm, omniClock))
    mainLoop.addAgents([omniClock, gateAgent])
    mainLoop.addAgents(allAgents)
    mainLoop.switch()
    print '%d all done' % rank

#
#     if rank == 0:
#         data = {'a': 7, 'b': 3.14}
#         comm.send(data, dest=1, tag=11)
#     else:
#         data = comm.recv(source=0, tag=11)
#         print '%s says %s' % (rank, data)

############
# Main hook
############

if __name__ == "__main__":
    main()
