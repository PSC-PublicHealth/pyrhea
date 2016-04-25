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

"""
This provides a simple test routine for patches.
"""

import sys
from random import choice, seed

import patches
import logging

logger = logging.getLogger(__name__)


class TestInteractant(patches.Interactant):
    pass


class TestAgent(patches.Agent):
    STATE_ATTARGET = 0
    STATE_MOVING = 1

    def __init__(self, name, patch, timeNow=0, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.targetAddr = None
        self.target = None
        self.fsmstate = TestAgent.STATE_ATTARGET

    def run(self, startTime):
        timeNow = startTime
        while True:
            if self.fsmstate == TestAgent.STATE_ATTARGET:
                timeNow = self.sleep(1)
                candidates = self.patch.serviceLookup('TestInteractant')
                if not candidates:
                    raise RuntimeError('Found no TestInteractants to visit!')
                nm, newTargetAddr = choice(candidates)
                self.targetAddr = newTargetAddr
                if self.target:
                    timeNow = self.target.unlock(self)
                    self.target = None
                self.fsmstate = TestAgent.STATE_MOVING
                if self.debug:
                    logger.debug('%s leaving for %s at %s' % (self.name, nm, timeNow))
            elif self.fsmstate == TestAgent.STATE_MOVING:
                nxtI, final = self.patch.getPathTo(self.targetAddr)
                if final:
                    self.fsmstate = TestAgent.STATE_ATTARGET
                    self.target = nxtI
                    if self.debug:
                        logger.debug('%s arrives at %s day %s' % (self.name, nxtI._name, timeNow))
                timeNow = nxtI.lock(self)
        if self.debug:
            return '%s is exiting at %s' % (self, timeNow)

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['fsmstate'] = self.fsmstate
        d['targetAddr'] = self.targetAddr
        d['target'] = self.target  # which should be None while in motion
        return d

    def __setstate__(self, d):
        patches.Agent.__setstate__(self, d)
        self.fsmstate = d['fsmstate']
        self.targetAddr = d['targetAddr']
        self.target = d['target']  # which is currently None


def createPerDayCB(patch, runDurationDays):
    def perDayCB(loop, timeNow):
        if timeNow > runDurationDays:
            patch.group.stop()
    return perDayCB


def describeSelf():
    print """
    This main provides diagnostics. -t and -d for trace and debug respectively.
    """


def main():
    trace = False
    debug = False
    deterministic = False

    for a in sys.argv[1:]:
        if a == '-d':
            debug = True
        elif a == '-t':
            trace = True
        elif a == '--deterministic':
            deterministic = True
        else:
            describeSelf()
            sys.exit('unrecognized argument %s' % a)

    if debug:
        logLevel = 'DEBUG'
    else:
        logLevel = 'INFO'

    comm = patches.getCommWorld()
    rank = comm.rank
    logging.basicConfig(format="%%(levelname)s:%%(name)s:rank%s:%%(message)s" % rank,
                        level=logLevel)

    if deterministic:
        seed(1234)

    patchGroup = patches.PatchGroup(comm, trace=trace, deterministic=deterministic)
    nPatches = 2
    for j in xrange(nPatches):

        patch = patches.Patch(patchGroup)

        patch.addInteractants([TestInteractant('%s_%d_%d' % (nm, rank, j), patch)
                               for nm in ['SubA', 'SubB', 'SubC']])

        allAgents = []
        for i in xrange(1000):
            debugThis = (i % 50 == 0)
            allAgents.append(TestAgent('Agent_%d_%d_%d' % (rank, j, i),
                                       patch, debug=debugThis))

        patch.addAgents(allAgents)
        patch.loop.addPerDayCallback(createPerDayCB(patch, 1000))
        patchGroup.addPatch(patch)
    logger.info('starting main loop')
    msg = patchGroup.start()
    logger.info('%d all done (from main) with msg %s' % (rank, msg))
    logging.shutdown()


############
# Main hook
############

if __name__ == "__main__":
    main()
