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
from random import choice, seed, random

import patches
import peopleplaces
import logging

logger = logging.getLogger(__name__)


class LocTypeBase(peopleplaces.Location):
    def __init__(self, name, patch, capacity):
        super(LocTypeBase, self).__init__(name, patch, capacity)
        self.grp = None

    def getReqQueueAddr(self):
        return self.grp.reqQueues[0].getGblAddr()

    def getNClients(self):
        return len(self._lockingAgentList)


class LocType0(LocTypeBase):
    pass


class LocType1(LocTypeBase):
    pass


class LocType2(LocTypeBase):
    pass


class LocType3(LocTypeBase):
    pass


class LocType4(LocTypeBase):
    pass


class LocType5(LocTypeBase):
    pass


class LocType6(LocTypeBase):
    pass

locTypeCycle = [LocType0, LocType1, LocType2, LocType3, LocType4, LocType5, LocType6]


class Walker(peopleplaces.Person):
    def getNewLocAddr(self, timeNow):
        """
        This method is called once each time the Person agent is active and returns a tuple
        of the form (newLocGlobalAddr, updatedTimeNow).  Returning a newLocGlobalAddr of
        self.locAddr (that is, the current value) indicates that the Person stays attached
        to the same location for this time slice.  Returning a newLocGlobalAddr
        of None signals 'death' and will result in self.handleDeath being called and the agent's
        thread exiting.  This method may not return for a long time as the agent waits for
        a new location to become available.
        """
        wantLocType = locTypeCycle[timeNow % len(locTypeCycle)]
        if random() < 0.01:
            return None, timeNow
        if isinstance(self.loc, wantLocType):
            return self.locAddr, timeNow
        else:
            facAddrList = [tpl[1] for tpl in self.patch.serviceLookup(wantLocType.__name__)
                           if tpl[1] != self.locAddr]
            newAddr = choice(facAddrList)
            return newAddr, timeNow

    def handleArrival(self, timeNow):
        """
        An opportunity to do bookkeeping on arrival at self.loc.  The Person agent has already
        locked the interactant self.loc.
        """
        self.loc.grp.altPop += 1

    def handleDeparture(self, timeNow):
        """
        An opportunity to do bookkeeping on arrival at self.loc.  The Person agent has already
        locked the interactant self.loc.
        """
        self.loc.grp.altPop -= 1

    def handleDeath(self, timeNow):
        """
        The name says it all.  Do any bookkeeping needed to deal with the death of this Person.
        After this method returns, a DepartureMessage will be sent to inform its current location
        of departure and the agent's run method will exit.
        """
        self.loc.grp.nDead += 1

    def getPostArrivalPauseTime(self, timeNow):
        """
        This allows the insertion of an extra pause on arrival at a new location.  It is useful
        for desynchronizing the activity cycles of agents at a location which does not sample
        all agents daily.
        """
        return 0

    def __getstate__(self):
        d = super(Walker, self).__getstate__()
        return d

    def __setstate__(self, d):
        super(Walker, self).__setstate__(d)


class LocManager(peopleplaces.Manager):
    pass


class LocGroup(peopleplaces.ManagementBase):
    def __init__(self, name, patch):
        super(LocGroup, self).__init__(name, patch, LocManager)
        self.locs = []
        self.altPop = 0
        self.nDead = 0

    def addLocs(self, locList):
        self.locs.extend(locList)
        for loc in locList:
            loc.grp = self


class MyPatch(patches.Patch):
    """
    Specialized so that we can collect some display data
    """
    def __init__(self, group, name=None, patchId=None):
        super(MyPatch, self).__init__(group, name, patchId)
        self.locGroups = []
        self.termsList = []  # used in printing output

    def addAgents(self, agentList):
        super(MyPatch, self).addAgents(agentList)
        for a in agentList:
            if isinstance(a, LocManager):
                self.locGroups.append(a.toManage)
            elif isinstance(a, Walker):
                a.loc.grp.altPop += 1


def createPerTickCB(patch, runDurationDays):
    def perTickCB(loop, timeNow, newTimeNow):
        assert isinstance(patch, MyPatch), 'You forgot to use MyPatch instances for patches'
        # Print the output terms from the last tick *before* the date change
        if timeNow != newTimeNow:
            print 'time is %s -> %s' % (timeNow, newTimeNow)
            for nm, ct, altPop, nDead in patch.termsList[-2]:
                print ('%s: %2d %2d %2d' % (nm, ct, altPop, nDead)),
            print ""   # newline
        terms = []
        for gp in patch.locGroups:
            terms.append((gp.name,
                          sum([loc.getNClients() for loc in gp.locs]),
                          gp.altPop, gp.nDead))
        terms.sort()
        # Keep the last two terms
        patch.termsList = patch.termsList[-1:].append(terms)
        if newTimeNow > runDurationDays:
            patch.group.stop()
    return perTickCB


def describeSelf():
    print """
    This main provides diagnostics. -t and -d for trace and debug respectively.
    """


def main():
    trace = False
    debug = False
    deterministic = False
    locCapacity = 100
    agentsPerPatch = 35
    locsPerPatch = 5
    patchesPerRank = 2
    runDuration = 30

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
    for j in xrange(patchesPerRank):

        patch = MyPatch(patchGroup)

        locList = []
        itrList = []
        groupList = []
        for i, LocTp in enumerate(locTypeCycle):
            locGroup = LocGroup('Grp_%d_%d' % (j, i), patch)
            theseLocs = [LocTp('loc_%s_%d_%d_%d' % (LocTp.__name__, j, i, k), patch, locCapacity)
                         for k in xrange(locsPerPatch)]
            locGroup.addLocs(theseLocs)
            locList.extend(theseLocs)
            itrList.extend(theseLocs + locGroup.getAllQueues())
            groupList.append(locGroup)
        patch.addInteractants(itrList)
        patch.addAgents([gp.manager for gp in groupList])
        agentList = []
        for i in xrange(agentsPerPatch):
            agentList.append(Walker('walker_%d_%d' % (j, i), patch, choice(locList)))
        patch.addAgents(agentList)

        # Use a PerTick callback rather than PerDay to make sure we catch the exact edge of the day
        patch.loop.addPerTickCallback(createPerTickCB(patch, runDuration))
        patchGroup.addPatch(patch)
    logger.info('starting main loop')
    msg = patchGroup.start()
    logger.info('%d all done (from main) with msg "%s"' % (rank, msg))
    logging.shutdown()


############
# Main hook
############

if __name__ == "__main__":
    main()
