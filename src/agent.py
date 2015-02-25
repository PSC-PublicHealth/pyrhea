#! /usr/bin/env python

_rhea_svn_id_="$Id$"

import abc
import types
from greenlet import greenlet
from random import randint


class AgentException(Exception):
    pass


class AgentSleepNowException(AgentException):
    pass


class Sequencer(object):

    def __init__(self):
        self._todayQueue = []
        self._tomorrowQueue = []
        self._timeNow = 0

    def __iter__(self):
        while True:
            if self._todayQueue:
                yield (self._todayQueue.pop(0), self._timeNow)
            elif self._tomorrowQueue:
                self._todayQueue = self._tomorrowQueue
                self._tomorrowQueue = []
                self._timeNow += 1
            else:
                break

    def enqueue(self, agent, whenInfo=0):
        assert isinstance(whenInfo, types.IntType)
        if whenInfo < self._timeNow:
            raise RuntimeError('Cannot go backwards in time to %d' % self.timeNow)
        elif whenInfo == self._timeNow:
            self._todayQueue.append(agent)
        else:
            self._tomorrowQueue.append(agent)

    def getTimeNow(self):
        return self._timeNow


class Agent(greenlet):
    def __init__(self, name):
        self.name = name

    def run(self, startTime):
        timeNow = startTime
        while True:
            print '%s new iter' % self
            fate = randint(0, len(interactants))
            if fate == len(interactants):
                print 'no lock for %s at %s' % (self.name, timeNow)
                timeNow = mainLoop.sleep(self, 0)  # yield thread
            else:
                timeNow = interactants[fate].lock(self)
                timeNow = mainLoop.sleep(self, 1)
                timeNow = interactants[fate].unlock(self)
        return '%s is exiting' % self

    def serialize(self):
        raise RuntimeError('agent cannot serialize')

    @staticmethod
    def deserialize(agentRep):
        raise RuntimeError('agent cannot deserialize')

    def __str__(self):
        return '<%s>' % self.name


class Interactant():
    def __init__(self, name):
        self._name = name
        self._lockingAgent = None
        self._lockQueue = []

    def lock(self, lockingAgent):
        timeNow = mainLoop.sequencer.getTimeNow()
        if ((self._lockingAgent is None and not self._lockQueue)
                or self._lockingAgent == lockingAgent):
            self._lockingAgent = lockingAgent
            print '%s fast lock of %s' % (lockingAgent, self._name)
            return timeNow
        else:
            self._lockQueue.append(lockingAgent)
            print '%s slow lock of %s' % (lockingAgent, self._name)
            timeNow = mainLoop.switch('%s is %d in %s queue' %
                                      (lockingAgent, len(self._lockQueue), self._name))
            return timeNow

    def unlock(self, oldLockingAgent):
        if self._lockingAgent != oldLockingAgent:
            raise RuntimeError('%s is not the lock of %s' % (oldLockingAgent, self._name))
        timeNow = mainLoop.sequencer.getTimeNow()
        if self._lockQueue:
            newAgent = self._lockQueue.pop(0)
            print '%s unlock of %s awakens %s' % (self._name, oldLockingAgent, newAgent)
            self._lockingAgent = newAgent
            mainLoop.sequencer.enqueue(newAgent, timeNow)
            mainLoop.sequencer.enqueue(oldLockingAgent, timeNow)
            timeNow = mainLoop.switch("%s and %s enqueued" % (newAgent, oldLockingAgent))
        else:
            print '%s fast unlock of %s' % (self._name, oldLockingAgent)
            self._lockingAgent = None
        return timeNow


interactants = [Interactant(nm) for nm in ['SubA', 'SubB', 'SubC']]


class MainLoop(greenlet):
    def __init__(self):
        self.sequencer = Sequencer()
        self.newAgents = []

    def addAgents(self, agentList):
        self.newAgents.extend(agentList)

    def run(self):
        counter = 0
        for a in self.newAgents:
            a.parent = self  # so dead agents return here
            self.sequencer.enqueue(a)
        for agent, timeNow in self.sequencer:
            print 'Stepping %s at %d' % (agent, timeNow)
            reply = agent.switch(timeNow)
            print 'Stepped %s at %d; reply was %s' % (agent, timeNow, reply)
            counter += 1
            if counter > 100000:
                print 'Safety exit!'
                break

    def sleep(self, agent, nDays):
        assert isinstance(nDays, types.IntType), 'nDays should be an integer'
        assert nDays >= 0, 'No sleeping for negative time'
        self.sequencer.enqueue(agent, self.sequencer.getTimeNow() + nDays)
        return self.switch('%s sleep %d days' % (agent, nDays))

allAgents = []
for i in xrange(20000):
    allAgents.append(Agent('Agent_%d' % i))

mainLoop = MainLoop()
mainLoop.addAgents(allAgents)
mainLoop.switch()
print 'all done'
