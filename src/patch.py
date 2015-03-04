#! /usr/bin/env python

_rhea_svn_id_="$Id$"

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
#         self.mem = np.zeros(2, dtype=np.int64)
#         self.win = MPI.Win.Create(self.mem, disp_unit=8, comm=self.comm)
        self.readMem = np.zeros(2, dtype=np.int64)
        self.maxInt = np.iinfo(np.int64).max
        self.writeMem = np.asarray([self.maxInt, self.maxInt], dtype=np.int64)

    def run(self, startTime):
        timeNow = startTime
        while True:
            self.win.Fence()
            self.win.Get(self.readMem, 0)
            self.win.Fence()
            print 'rank %d: %s has time range %s at local time %s' % \
                (self.comm.rank, self.name, self.getRange(), timeNow)
            if self.getRange()[0] > 100:
                print 'rank %d: run complete' % self.comm.rank
                MPI.Finalize()
                sys.exit(0)
            if self.comm.rank == 0:
                self.writeMem = np.asarray([self.maxInt, self.maxInt])
                self.win.Put(self.writeMem, 0)
            self.win.Fence()
            timeNow = self.ownerLoop.sleep(self, 0)
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

    def run(self, startTime):
        timeNow = startTime
        while True:
            if self.cycleCounter:  # don't do this on pass 0
                self.win.Fence()   # win is now closed
                print '%s at time %s cycle %d says %s' % (self.comm.rank, timeNow,
                                                          self.cycleCounter, self.mem)
            self.win.Fence()  # re-open the window
            v = np.ones(1, dtype=np.int64)
            for i in xrange(self.comm.size):
                self.win.Accumulate(v, i, target=self.comm.rank, op=MPI.SUM)
            timeNow = self.ownerLoop.sleep(self, 0)  # yield thread
            self.cycleCounter += 1


def describeSelf():
    print """This main provides diagnostics. -v and -d for verbose and debug respectively."""


def createPerTickCB(comm, omniClock):
    def tick(myLoop):
        tMin, tMax = myLoop.sequencer.getTimeRange()
        print 'tick! from %d: %s %s' % (comm.rank, tMin, tMax)
        omniClock.setMyRange(tMin, tMax)
    return tick


def main():
    global verbose, debug

    mainLoop = agent.MainLoop()

    interactants = [agent.Interactant(nm, mainLoop) for nm in ['SubA', 'SubB', 'SubC']]

    class TestAgent(agent.Agent):
        def run(self, startTime):
            timeNow = startTime
            while True:
                # print '%s new iter' % self
                fate = randint(0, len(interactants))
                if fate == len(interactants):
                    # print 'no lock for %s at %s' % (self.name, timeNow)
                    timeNow = self.ownerLoop.sleep(self, 0)  # yield thread
                else:
                    timeNow = interactants[fate].lock(self)
                    timeNow = self.ownerLoop.sleep(self, 1)
                    timeNow = interactants[fate].unlock(self)
            return '%s is exiting at %s' % (self, timeNow)

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
    allAgents = []
    for i in xrange(100):
        allAgents.append(TestAgent('Agent_%d' % i, mainLoop))

    omniClock = OmniClock(comm, mainLoop)
    mainLoop.addPerTickCallback(createPerTickCB(comm, omniClock))
    mainLoop.addAgents([omniClock])
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
