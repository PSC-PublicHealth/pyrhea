#! /usr/bin/env python

_rhea_svn_id_ = "$Id$"

import sys
from mpi4py import MPI
from agentworld import *
from random import randint


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


class FacilityManager(Agent):
    def __init__(self, name, mainLoop):
        Agent.__init__(self, name, mainLoop)
        self.wardList = []

    def addWard(self, ward):
        self.wardList.append(ward)


class Ward(Interactant):
    def __init__(self, name, mainLoop, nBeds):
        Interactant.__init__(self, name, mainLoop)
        self._nBeds = nBeds


class PatientAgent(Agent):
    pass


class BedRequest(Agent):
    pass


class TestAgent(Agent):
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
        d = Agent.__getstate__(self)
        return d

    def __setstate__(self, stateDict):
        Agent.__setstate__(self, stateDict)

def describeSelf():
    print """This should write some documentation"""

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
    assert comm.size == 3, "This program runs only with np=3"

    mainLoop = MainLoop('Mainloop_%s' % rank)
    omniClock = OmniClock(comm, mainLoop)

    gateAgent = GateAgent(comm, mainLoop)
    for i in xrange(comm.size):
        if i != rank:
            gateIn = GateIn(("GateIn_%d_%d" % (rank, i)), i, comm, mainLoop, debug=False)
            gateAgent.addGate(gateIn)
            gateOut = GateOut(("GateOut_%d_%d" % (i, rank)), i, comm, mainLoop, debug=True)
            gateAgent.addGate(gateOut)

    allAgents = [omniClock, gateAgent]

    for tier, nWards, nBedsPerWard in [(1, 3, 20),
                                       (2, 2, 100),
                                       (3, 1, 1000)]:
        tierManager = FacilityManager('Rank%d_Tier%d_Manager' % (rank, tier), mainLoop)
        allAgents.append(tierManager)
        for i in xrange(nWards):
            tierManager.addWard(Ward('Rank%d_Tier%d_Ward%d' % (rank, tier, i),
                                     mainLoop,
                                     nBedsPerWard))

#     interactants = [Interactant(nm, mainLoop) for nm in ['SubA', 'SubB', 'SubC']]
#
#     for i in xrange(100):
#         allAgents.append(TestAgent('Agent_%d_%d' % (rank, i), mainLoop, debug=True))

    mainLoop.addPerTickCallback(createPerTickCB(comm, omniClock))
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
