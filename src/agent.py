#! /usr/bin/env python

_rhea_svn_id_="$Id$"

import sys
import types
from greenlet import greenlet
#from gevent import Greenlet as greenlet
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
        if self._timeNow in self._timeQueues:
            return len([a for a in self._timeQueues[self._timeNow] if not a.timeless])
        else:
            return 0

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
        if self.doneWithToday():
            # print '%s: bump time %s -> %s' % (self._name, self._timeNow, self._timeNow+1)
            oldDay = self._timeQueues[self._timeNow]
            del self._timeQueues[self._timeNow]
            self._timeNow += 1
            if self._timeNow not in self._timeQueues:
                self._timeQueues[self._timeNow] = []
            self._timeQueues[self._timeNow].extend(oldDay)

    def doneWithToday(self):
        """
        If the only agents still in today's sequencer loop are marked 'timeless', this
        will return True; otherwise False.  Note that this condition can be reversed if
        a new agent is inserted into today's loop.
        """
        for iact in Interactant.getLiveList():
            if (iact._lockingAgent is not None and iact._lockingAgent.timeless
                    and iact.getNWaiting()):
                # print ('doneWithToday is false because %s has %d waiting: %s' %
                #        (iact, iact.getNWaiting(), iact.getWaitingDetails()))
                return False
        return all([a.timeless for a in self._timeQueues[self._timeNow]])


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

    def kill(self):
        """
        Cause this greenlet to throw GreenletExit.  If it is not the current greenlet,
        this greenlet's parent is set to the current greenlet before throwing the exception,
        causing execution to return to the current (calling) greenlet after the exception
        is thrown.  If the greenlet is the current greenlet, execution passes to its parent.
        """
        if self != greenlet.getcurrent():
            self.parent = greenlet.getcurrent()
        self.throw()


class Interactant():
    counter = 0
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
        self.id = Interactant.counter
        Interactant.counter += 1

    def getNWaiting(self):
        """This returns the count of waiting agents for which timeless is false"""
        return self._nEnqueued

    def getWaitingDetails(self):
        """Returns a dict of typeName:nOfThisType entries"""
        result = {}
        for a in self._lockQueue:
            nm = type(a).__name__
            if nm in result:
                result[nm] += 1
            else:
                result[nm] = 1
        return result

    def __str__(self):
        return '<%s>' % self._name

    def lock(self, lockingAgent, debug=True):
        """
        Agents always lock interactants before modifying their state.  This can be thought of as
        'docking with' the interactant.  Only one agent can hold a lock for a normal interactant
        and be active; any agent that subsequently tries to lock the same interactant will be
        suspended until the firs agent unlocks the interactant.  (However, see MultiInteractant).
        Thus interactants can serve as queues; agents can enqueue themselves by locking the
        interactant and some other agent can modify them while they are in the locked state.
        """
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
            # print '%s locked %s; nEnqueued = %d' % (self, lockingAgent, self._nEnqueued)
            if debug and self._debug and lockingAgent.debug:
                print '%s slow lock of %s (%d in queue)' % \
                    (lockingAgent, self._name, self._nEnqueued)
            timeNow = self._ownerLoop.switch('%s is %d in %s queue' %
                                             (lockingAgent, len(self._lockQueue), self._name))
            return timeNow

    def unlock(self, oldLockingAgent):
        """
        This method will typically be called by an active agent which holds a lock on the
        interactant.  The lock is broken, causing the first agent which is suspended waiting
        for a lock to become active.
        """
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

    def awaken(self, agent):
        """
        The agent is expected to be sleeping in this interactant's _lockQueue.  Calling
        awaken(agent) removes the agent from _lockQueue and re-inserts it into the
        main loop, so that it will resume running in its turn.  Essentially, the agent
        has been unlocked from the interactant in such a way that it will never become
        active in a locked state.  The agent which calls awaken on a locked agent does
        not yield its thread; the awakened agent simply joins the queue to become active
        in its turn.  It is an error to awaken an agent which is not suspended in the
        interactant's wait queue.
        """
        timeNow = self._ownerLoop.sequencer.getTimeNow()
        if agent not in self._lockQueue:
            raise RuntimeError("%s does not hold %s in its lock queue; cannot awaken" %
                               (self._name, agent.name))
        # print 'Interactant %s removes %s' % (self, agent)
        self._lockQueue.remove(agent)
        if not agent.timeless:
            self._nEnqueued -= 1
        # print 'Interactant %s nEnqueued %d' % (self, self._nEnqueued)
        if self._debug:
            print ('%s removes %s from lock queue and awakens it (%d still in queue)' %
                   (self._name, agent.name, self._nEnqueued))
        self._ownerLoop.sequencer.enqueue(agent, timeNow)
        return agent

    def isLocked(self, agent):
        """
        Returns True if this interactant is currently locked by the given agent, whether the
        agent is active or has been suspended in the interactant's lock wait queue.  It is
        unlikely that a scenario will arise in which an agent will ever have to test whether
        it holds a lock, but the method exists for completeness of the API.
        """
        return (self._lockingAgent == agent or agent in self._lockQueue)


class MultiInteractant(Interactant):
    """
    A MultiInteractant functions like a generic Interactant, except that more than one
    agent can lock it simultaneously and yet remain active.
    """

    def __init__(self, name, count, ownerLoop, debug=False):
        """
        Create an interactant that can simulaneously hold locks for 'count' active
        agents.
        """
        Interactant.__init__(self, name, ownerLoop, debug)
        self._nLocks = count
        self._lockingAgentList = []
        self._debug = debug

    def lock(self, lockingAgent, debug=False):
        """
        Works like the lock() method of a standard Interactant, except that the first
        'count' agents to lock the interactant remain active.
        """
        timeNow = self._ownerLoop.sequencer.getTimeNow()
        if lockingAgent in self._lockingAgentList:
            if debug and self._debug and lockingAgent.debug:
                print '%s already locked by %s' % (self._name, lockingAgent)
            return timeNow
        elif len(self._lockingAgentList) < self._nLocks:
            self._lockingAgentList.append(lockingAgent)
            if debug and self._debug and lockingAgent.debug:
                print '%s fast locked by %s' % (self._name, lockingAgent)
            return timeNow
        else:
            self._lockQueue.append(lockingAgent)
            if not lockingAgent.timeless:
                self._nEnqueued += 1
            if debug and self._debug and lockingAgent.debug:
                print '%s slow lock by %s (%d in queue)' % \
                    (self._name, lockingAgent, self._nEnqueued)
            timeNow = self._ownerLoop.switch('%s is %d in %s queue' %
                                             (lockingAgent, len(self._lockQueue), self._name))
            return timeNow

    def unlock(self, oldLockingAgent):
        if oldLockingAgent not in self._lockingAgentList:
            raise RuntimeError('%s is not a lock of %s' % (oldLockingAgent, self._name))
        timeNow = self._ownerLoop.sequencer.getTimeNow()
        self._lockingAgentList.remove(oldLockingAgent)
        if self._lockQueue:
            newAgent = self._lockQueue.pop(0)
            if not newAgent.timeless:
                self._nEnqueued -= 1
            if self._debug:
                print '%s unlock of %s awakens %s (%d still in queue)' % \
                    (self._name, oldLockingAgent, newAgent, self._nEnqueued)
            self._lockingAgentList.append(newAgent)
            self._ownerLoop.sequencer.enqueue(newAgent, timeNow)
            self._ownerLoop.sequencer.enqueue(oldLockingAgent, timeNow)
            timeNow = self._ownerLoop.switch("%s and %s enqueued" % (newAgent, oldLockingAgent))
        else:
            if self._debug:
                print '%s fast unlock of %s' % (self._name, oldLockingAgent)
        return timeNow

    def isLocked(self, agent):
        return (agent in self._lockingAgentList or agent in self._lockQueue)

    def __str__(self):
        return '<%s (%d of %d)>' % (self._name, len(self._lockingAgentList), self._nLocks)

    @property
    def nFree(self):
        return self._nLocks - len(self._lockingAgentList)


class MainLoop(greenlet):
    class ClockAgent(Agent):
        def __init__(self, ownerLoop):
            Agent.__init__(self, 'ClockAgent', ownerLoop)
            self.timeless = True

        def run(self, timeNow):
            while True:
                if not self.ownerLoop.dateFrozen:
                    self.ownerLoop.sequencer.bumpIfAllTimeless()
                newTimeNow = self.sleep(0)  # yield thread
                for cb in self.ownerLoop.perTickCallbacks:
                    cb(self, timeNow, newTimeNow)
                if newTimeNow != timeNow:
                    print '###### %s ClockAgent: time is now %s' % (self.ownerLoop.name, newTimeNow)
                    timeNow = newTimeNow

    @staticmethod
    def everyEventCB(loop, timeNow):
        loop.counter += 1
        if loop.counter > loop.safety:
            print '%s: safety exit' % loop.name
            loop.parent.switch(loop.counter)
            loop.counter = 0

    def __init__(self, name=None, safety=None):
        self.newAgents = [MainLoop.ClockAgent(self)]
        self.perTickCallbacks = []
        self.perEventCallbacks = []
        self.safety = safety  # After how many ticks to bail, if any
        assert safety is None or isinstance(safety, types.IntType)
        if name is None:
            self.name = 'MainLoop'
        else:
            self.name = name
        self.sequencer = Sequencer(self.name + ".Sequencer")
        self.dateFrozen = False
        self.counter = 0
        if self.safety is not None:
            self.addPerEventCallback(MainLoop.everyEventCB)

    def addAgents(self, agentList):
        assert all([a.ownerLoop == self for a in agentList]), \
            "%s: Tried to add a foreign agent!" % self.name
        self.newAgents.extend(agentList)

    def addPerTickCallback(self, cb):
        self.perTickCallbacks.append(cb)

    def addPerEventCallback(self, cb):
        self.perEventCallbacks.append(cb)

    def freezeDate(self):
        self.dateFrozen = True

    def unfreezeDate(self):
        self.dateFrozen = False

    def run(self):
        for a in self.newAgents:
            a.parent = self  # so dead agents return here
            self.sequencer.enqueue(a)
        self.newAgents = []
        for agent, timeNow in self.sequencer:
            # print '%s Stepping %s at %d' % (self.name, agent, timeNow)
            for cb in self.perEventCallbacks:
                cb(self, timeNow)
            reply = agent.switch(timeNow)  # @UnusedVariable
            # print 'Stepped %s at %d; reply was %s' % (agent, timeNow, reply)
        print 'sequencer ran dry'

    def sleep(self, agent, nDays):
        assert isinstance(nDays, types.IntType), 'nDays should be an integer'
        assert nDays >= 0, 'No sleeping for negative time'
        self.sequencer.enqueue(agent, self.sequencer.getTimeNow() + nDays)
        return self.switch('%s: %s sleep %d days' % (self.name, agent, nDays))

    def printCensus(self):
        print '%s: Census at time %s:' % (self.name, self.sequencer.getTimeNow())
        for iact in Interactant.getLiveList():
            print '    %s : %s' % (iact._name, iact.getWaitingDetails())
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
