#! /usr/bin/env python

_rhea_svn_id_="$Id$"

import sys
import types
from greenlet import greenlet
from random import randint


class Sequencer(object):

    def __init__(self):
        self._timeQueues = {}
        self._timeNow = 0

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
            raise RuntimeError('Cannot go backwards in time from %d to %d' %
                               (self.timeNow, whenInfo))
        else:
            if whenInfo not in self._timeQueues:
                self._timeQueues[whenInfo] = []
            self._timeQueues[whenInfo].append(agent)

    def getTimeNow(self):
        return self._timeNow

    def getTimeRange(self):
        return (self._timeNow, max(self._timeQueues.keys()))

    def bumpIfAllTimeless(self):
        """
        If all the agents for 'today' are marked timeless, shift them to 'tomorrow' and move
        time forward by a day.
        """
        if all([a.timeless for a in self._timeQueues[self._timeNow]]):
            print 'Bump!'
            oldDay = self._timeQueues[self._timeNow]
            del self._timeQueues[self._timeNow]
            self._timeNow += 1
            self._timeQueues[self._timeNow].extend(oldDay)


class Agent(greenlet):
    def __init__(self, name, ownerLoop):
        self.name = name
        self.ownerLoop = ownerLoop
        self.timeless = False

    def run(self, startTime):
        raise RuntimeError('Derived class must subclass this method!')

    def serialize(self):
        raise RuntimeError('agent cannot serialize')

    @staticmethod
    def deserialize(agentRep):
        raise RuntimeError('agent cannot deserialize')

    def sleep(self, deltaTime):
        self.ownerLoop.sleep(self, deltaTime)

    def __str__(self):
        return '<%s>' % self.name


class Interactant():
    def __init__(self, name, ownerLoop, debug=False):
        self._name = name
        self._ownerLoop = ownerLoop
        self._lockingAgent = None
        self._lockQueue = []
        self._debug = debug

    def lock(self, lockingAgent):
        timeNow = self._ownerLoop.sequencer.getTimeNow()
        if ((self._lockingAgent is None and not self._lockQueue)
                or self._lockingAgent == lockingAgent):
            self._lockingAgent = lockingAgent
            if self._debug:
                print '%s fast lock of %s' % (lockingAgent, self._name)
            return timeNow
        else:
            self._lockQueue.append(lockingAgent)
            if self._debug:
                print '%s slow lock of %s' % (lockingAgent, self._name)
            timeNow = self._ownerLoop.switch('%s is %d in %s queue' %
                                             (lockingAgent, len(self._lockQueue), self._name))
            return timeNow

    def unlock(self, oldLockingAgent):
        if self._lockingAgent != oldLockingAgent:
            raise RuntimeError('%s is not the lock of %s' % (oldLockingAgent, self._name))
        timeNow = self._ownerLoop.sequencer.getTimeNow()
        if self._lockQueue:
            newAgent = self._lockQueue.pop(0)
            if self._debug:
                print '%s unlock of %s awakens %s' % (self._name, oldLockingAgent, newAgent)
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
                self.ownerLoop.sequencer.bumpIfAllTimeless()
                newTimeNow = self.sleep(0)  # yield thread
                if newTimeNow != timeNow:
                    # print 'ClockAgent: time is now %s' % newTimeNow
                    timeNow = newTimeNow
                for cb in self.ownerLoop.perTickCallbacks:
                    cb(self.ownerLoop)

    def __init__(self, safety=None):
        self.sequencer = Sequencer()
        self.newAgents = [MainLoop.ClockAgent(self)]
        self.perTickCallbacks = []
        self.safety = safety  # After how many ticks to bail, if any
        assert safety is None or isinstance(safety, types.IntType)

    def addAgents(self, agentList):
        assert all([a.ownerLoop == self for a in agentList]), "Tried to add a foreign agent!"
        self.newAgents.extend(agentList)

    def addPerTickCallback(self, cb):
        self.perTickCallbacks.append(cb)

    def run(self):
        counter = 0
        for a in self.newAgents:
            a.parent = self  # so dead agents return here
            self.sequencer.enqueue(a)
        for agent, timeNow in self.sequencer:
            # print 'Stepping %s at %d' % (agent, timeNow)
            reply = agent.switch(timeNow)
            # print 'Stepped %s at %d; reply was %s' % (agent, timeNow, reply)
            counter += 1
            if self.safety is not None and counter > self.safety:
                print 'Safety exit!'
                break

    def sleep(self, agent, nDays):
        assert isinstance(nDays, types.IntType), 'nDays should be an integer'
        assert nDays >= 0, 'No sleeping for negative time'
        self.sequencer.enqueue(agent, self.sequencer.getTimeNow() + nDays)
        return self.switch('%s sleep %d days' % (agent, nDays))


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
