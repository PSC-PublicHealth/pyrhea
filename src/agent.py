#! /usr/bin/env python

_rhea_svn_id_="$Id$"

import sys
import types
from greenlet import greenlet
from random import randint
import weaklist

try:
    import faulthandler
    faulthandler.enable()
    if sys.excepthook != sys.__excepthook__:
        print("Warning: 3rd party exception hook is active")
        if sys.excepthook.__name__ == 'apport_excepthook':
            print("         killing Ubuntu's Apport hook")
            sys.excepthook = sys.__excepthook__
except:
    pass

class Sequencer(object):

    def __init__(self, name):
        self._timeQueues = {}
        self._timeNow = 0
        self._name = name

    def __iter__(self):
        while self._timeQueues:
            todayQueue = self._timeQueues[self._timeNow]
            if todayQueue:
                yield (todayQueue.pop(0), self._timeNow)
            else:
                if self._timeNow in self._timeQueues:
                    del self._timeQueues[self._timeNow]
                self._timeNow += 1

    def enqueue(self, agent, whenInfo=0):
        assert isinstance(whenInfo, types.IntType)
        if whenInfo < self._timeNow:
            raise RuntimeError('%s: Cannot go backwards in time from %d to %d' %
                               (self._name, self._timeNow, whenInfo))
        else:
            if whenInfo not in self._timeQueues:
                self._timeQueues[whenInfo] = []
            self._timeQueues[whenInfo].append(agent)

    def getTimeNow(self):
        return self._timeNow

    def getNWaitingNow(self):
        return len([a for a in self._timeQueues[self._timeNow] if not a.timeless])

    def getTimeRange(self):
        t = self._timeNow
        tMax = max(self._timeQueues.keys())
        tMin = None
        for iact in Interactant.getLiveList():
            if (iact._lockingAgent is not None and iact._lockingAgent.timeless
                    and iact.getNWaiting()):
                tMin = t
                # print '%s has %d waiting' % (iact, iact.getNWaiting())
                break
        while tMin is None:
            if self._timeQueues[t]:
                for a in self._timeQueues[t]:
                    if not a.timeless:
                        tMin = t
                        break
            assert t <= tMax, "Lost in time"
            t += 1
        return (tMin, tMax)

    def bumpIfAllTimeless(self):
        """
        If all the agents for 'today' are marked timeless, shift them to 'tomorrow' and move
        time forward by a day.
        """
        for iact in Interactant.getLiveList():
            if (iact._lockingAgent is not None and iact._lockingAgent.timeless
                    and iact.getNWaiting()):
                return
        if all([a.timeless for a in self._timeQueues[self._timeNow]]):
            print '%s: bump time %s -> %s' % (self._name, self._timeNow, self._timeNow+1)
            oldDay = self._timeQueues[self._timeNow]
            del self._timeQueues[self._timeNow]
            self._timeNow += 1
            self._timeQueues[self._timeNow].extend(oldDay)


class Agent(greenlet):
    def __init__(self, name, ownerLoop, debug=False):
        self.name = name
        self.ownerLoop = ownerLoop
        self.timeless = False
        self.debug = debug

    def run(self, startTime):
        raise RuntimeError('Derived class must subclass this method!')

    def __getstate__(self):
        return {'name': self.name, 'timeless': self.timeless,
                'debug': self.debug}

    def __setstate__(self, stateDict):
        for k, v in stateDict.items():
            setattr(self, k, v)

    def sleep(self, deltaTime):
        return self.ownerLoop.sleep(self, deltaTime)

    def __str__(self):
        return '<%s>' % self.name


class Interactant():
    _liveInstances = weaklist.WeakList()

    @classmethod
    def getLiveList(cls):
        return cls._liveInstances

    def __init__(self, name, ownerLoop, debug=False):
        self._name = name
        self._ownerLoop = ownerLoop
        self._lockingAgent = None
        self._lockQueue = []
        self._debug = debug
        self._nEnqueued = 0  # counts only things which are not 'timeless'
        self._liveInstances.append(self)

    def getNWaiting(self):
        """This returns the count of waiting agents for which timeless is false"""
        return self._nEnqueued

    def __str__(self):
        return '<%s>' % self._name

    def lock(self, lockingAgent, debug=True):
        timeNow = self._ownerLoop.sequencer.getTimeNow()
        if ((self._lockingAgent is None and not self._lockQueue)
                or self._lockingAgent == lockingAgent):
            self._lockingAgent = lockingAgent
            if debug and self._debug and lockingAgent.debug:
                print '%s fast lock of %s' % (lockingAgent, self._name)
            return timeNow
        else:
            self._lockQueue.append(lockingAgent)
            if not lockingAgent.timeless:
                self._nEnqueued += 1
            if debug and self._debug and lockingAgent.debug:
                print '%s slow lock of %s (%d in queue)' % \
                    (lockingAgent, self._name, self._nEnqueued)
            timeNow = self._ownerLoop.switch('%s is %d in %s queue' %
                                             (lockingAgent, len(self._lockQueue), self._name))
            return timeNow

    def unlock(self, oldLockingAgent):
        if self._lockingAgent != oldLockingAgent:
            raise RuntimeError('%s is not the lock of %s' % (oldLockingAgent, self._name))
        timeNow = self._ownerLoop.sequencer.getTimeNow()
        if self._lockQueue:
            newAgent = self._lockQueue.pop(0)
            if not newAgent.timeless:
                self._nEnqueued -= 1
            if self._debug:
                print '%s unlock of %s awakens %s (%d still in queue)' % \
                    (self._name, oldLockingAgent, newAgent, self._nEnqueued)
            self._lockingAgent = newAgent
            self._ownerLoop.sequencer.enqueue(newAgent, timeNow)
            self._ownerLoop.sequencer.enqueue(oldLockingAgent, timeNow)
            timeNow = self._ownerLoop.switch("%s and %s enqueued" % (newAgent, oldLockingAgent))
        else:
            if self._debug:
                print '%s fast unlock of %s' % (self._name, oldLockingAgent)
            self._lockingAgent = None
        return timeNow


class MainLoop(greenlet):
    class ClockAgent(Agent):
        def __init__(self, ownerLoop):
            Agent.__init__(self, 'ClockAgent', ownerLoop)
            self.timeless = True

        def run(self, timeNow):
            while True:
                if not self.ownerLoop.timeFrozen:
                    self.ownerLoop.sequencer.bumpIfAllTimeless()
                newTimeNow = self.sleep(0)  # yield thread
                for cb in self.ownerLoop.perTickCallbacks:
                    cb(self, timeNow, newTimeNow)
                if newTimeNow != timeNow:
                    print '%s ClockAgent: time is now %s' % (self.ownerLoop.name, newTimeNow)
                    timeNow = newTimeNow

    def __init__(self, name=None, safety=None):
        self.newAgents = [MainLoop.ClockAgent(self)]
        self.perTickCallbacks = []
        self.safety = safety  # After how many ticks to bail, if any
        assert safety is None or isinstance(safety, types.IntType)
        if name is None:
            self.name = 'MainLoop'
        else:
            self.name = name
        self.sequencer = Sequencer(self.name + ".Sequencer")
        self.timeFrozen = False

    def addAgents(self, agentList):
        assert all([a.ownerLoop == self for a in agentList]), \
            "%s: Tried to add a foreign agent!" % self.name
        self.newAgents.extend(agentList)

    def addPerTickCallback(self, cb):
        self.perTickCallbacks.append(cb)

    def freezeTime(self):
        self.timeFrozen = True

    def unfreezeTime(self):
        self.timeFrozen = False

    def run(self):
        counter = 0
        for a in self.newAgents:
            a.parent = self  # so dead agents return here
            self.sequencer.enqueue(a)
        self.newAgents = []
        for agent, timeNow in self.sequencer:
            print '%s Stepping %s at %d' % (self.name, agent, timeNow)
            reply = agent.switch(timeNow)  # @UnusedVariable
            # print 'Stepped %s at %d; reply was %s' % (agent, timeNow, reply)
            counter += 1
            if self.safety is not None and counter > self.safety:
                print '%s: finished my block' % self.name
                self.parent.switch(counter)
                counter = 0

    def sleep(self, agent, nDays):
        assert isinstance(nDays, types.IntType), 'nDays should be an integer'
        assert nDays >= 0, 'No sleeping for negative time'
        self.sequencer.enqueue(agent, self.sequencer.getTimeNow() + nDays)
        return self.switch('%s: %s sleep %d days' % (self.name, agent, nDays))

    def printCensus(self):
        print '%s: Census at time %s:' % (self.name, self.sequencer.getTimeNow())
        for iact in Interactant.getLiveList():
            print '    %s : %d' % (iact._name, iact.getNWaiting())
        print '    main loop now : %d' % self.sequencer.getNWaitingNow()

    def __str__(self):
        return '<%s>' % self.name


def describeSelf():
    print """This main provides diagnostics. -v and -d for verbose and debug respectively."""


def main():
    "This is a simple test routine which takes csv files as arguments"
    global verbose, debug

    mainLoop = MainLoop(safety=10000)

    interactants = [Interactant(nm, mainLoop) for nm in ['SubA', 'SubB', 'SubC']]

    class TestAgent(Agent):
        def run(self, startTime):
            timeNow = startTime
            while True:
                print '%s new iter' % self
                fate = randint(0, len(interactants))
                if fate == len(interactants):
                    print 'no lock for %s at %s' % (self.name, timeNow)
                    timeNow = self.sleep(0)  # yield thread
                else:
                    timeNow = interactants[fate].lock(self)
                    timeNow = self.sleep(1)
                    timeNow = interactants[fate].unlock(self)
            return '%s is exiting' % self

    for a in sys.argv[1:]:
        if a == '-v':
            verbose = True
        elif a == '-d':
            debug = True
        else:
            describeSelf()
            sys.exit('unrecognized argument %s' % a)

    allAgents = []
    for i in xrange(20000):
        allAgents.append(TestAgent('Agent_%d' % i, mainLoop))

    mainLoop.addAgents(allAgents)
    mainLoop.switch()
    print 'all done'

############
# Main hook
############

if __name__ == "__main__":
    main()
