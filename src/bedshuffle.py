#! /usr/bin/env python

_rhea_svn_id_ = "$Id$"

import sys
from random import randint

import patches


class FacilityManager(patches.Agent):
    def __init__(self, name, patch):
        patches.Agent.__init__(self, name, patch)
        self.wardList = []

    def addWard(self, ward):
        self.wardList.append(ward)


class Ward(patches.Interactant):
    def __init__(self, name, patch, nBeds):
        patches.Interactant.__init__(self, name, patch)
        self._nBeds = nBeds


class PatientAgent(patches.Agent):
    pass


class BedRequest(patches.Agent):
    pass


class TestAgent(patches.Agent):
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
        d = patches.Agent.__getstate__(self)
        return d

    def __setstate__(self, stateDict):
        patches.Agent.__setstate__(self, stateDict)


def describeSelf():
    print """This should write some documentation"""


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

    comm = patches.getCommWorld()
    rank = comm.rank

    patchGroup = patches.PatchGroup(comm, trace=trace)
    nPatches = 2
    for j in xrange(nPatches):

        patch = patches.Patch(patchGroup)

        patch.setInteractants([patches.Interactant('%s_%d_%d' % (nm, rank, j), patch)
                               for nm in ['SubA', 'SubB', 'SubC']])

        # Fully connect everything
        for r in xrange(comm.size):
            for jj in xrange(nPatches):
                if not (r == comm.rank and j == jj):
                    friend = patches.GblAddr(r, jj)
                    patch.addGateTo(friend)
                    patch.addGateFrom(friend)

        allAgents = []
        for i in xrange(100):
            allAgents.append(TestAgent('Agent_%d_%d_%d' % (rank, j, i),
                                       patch, debug=debug))

        patch.addAgents(allAgents)
        patchGroup.addPatch(patch)
    print 'Rank %s reaches barrier' % rank
    patchGroup.barrier()
    print 'Rank %s leaves barrier' % rank
    patchGroup.start()
    print '%d all done (from main)' % rank


############
# Main hook
############

if __name__ == "__main__":
    main()

