#! /usr/bin/env python

_rhea_svn_id_ = "$Id$"

import sys
from random import randint, shuffle
import types
from collections import defaultdict

import patches
from netinterface import GblAddr


class Ward(patches.MultiInteractant):
    def __init__(self, name, patch, tier, nBeds):
        patches.MultiInteractant.__init__(self, name, nBeds, patch)
        self.tier = tier
        self.fac = None


class BedRequestQueue(patches.Interactant):
    pass


class HoldQueue(patches.Interactant):
    def __init__(self, name, patch, debug=False):
        patches.Interactant.__init__(self, name, patch, debug=debug)
        self.keyCounter = 0
        self.heldDict = {}

    def getUniqueKey(self):
        key = self.keyCounter
        self.keyCounter += 1
        return key

    def lock(self, lockingAgent, key=None, debug=True):
        if key is not None:
            self.heldDict[key] = lockingAgent
        return patches.Interactant.lock(self, lockingAgent)

    def awaken(self, agentOrKey):
        if isinstance(agentOrKey, patches.Agent):
            agent = agentOrKey
        elif agentOrKey in self.heldDict:
            agent = self.heldDict[agentOrKey]
            del self.heldDict[agentOrKey]
        else:
            raise RuntimeError("%s does not know the agent or key %s" % (self._name, agentOrKey))
        return patches.Interactant.awaken(self, agent)


class FacilityManager(patches.Agent):
    def __init__(self, name, patch, facility):
        patches.Agent.__init__(self, name, patch)
        self.fac = facility
        self.timeless = True

    def run(self, startTime):
        timeNow = startTime
        while True:
            while self.fac.reqQueue._lockQueue:
                req = self.fac.reqQueue._lockQueue[0]
                tier = req.tier
                if tier not in self.fac.wardDict:
                    self.fac.wardDict[tier] = []
                if isinstance(req, BedRequest):
                    ward = None
                    for idx, (w, n) in enumerate(self.fac.wardDict[tier]):
                        if n > 0:
                            self.fac.wardDict[tier].pop(idx)
                            self.fac.wardDict[tier].append((w, n-1))
                            ward = w
                    if ward is None:
                        # print '%s: no beds available for %s' % (self.name, req.name)
                        req.fsmstate = BedRequest.STATE_DENIEDWARD
                    else:
                        # print '%s: found a bed for %s' % (self.name, req.name)
                        req.bedWard = ward.getGblAddr()
                        req.fsmstate = BedRequest.STATE_GOTWARD
                    self.fac.reqQueue.awaken(req)
                elif isinstance(req, DepartureMsg):
                    ward = None
                    for idx, (w, n) in enumerate(self.fac.wardDict[tier]):
                        if w.getGblAddr() == req.wardAddr:
                            ward = w
                            self.fac.wardDict[tier].pop(idx)
                            self.fac.wardDict[tier].append((w, n+1))
                    if ward is None:
                        raise RuntimeError("%s: I do not own the ward at %s" %
                                           self.name, self.wardAddr)
                    else:
                        # print '%s: incremented ward %s' % (self.name, ward._name)
                        pass
                    self.fac.reqQueue.awaken(req)
                else:
                    raise RuntimeError("%s unexpectedly got the message %s" %
                                       (self.name, req.name))
            timeNow = self.sleep(0)


class Facility(object):
    def __init__(self, name, patch):
        self.name = name
        self.manager = FacilityManager(name + '_Mgr', patch, self)
        self.reqQueue = BedRequestQueue(name+'_rQ', patch)
        self.reqQueue.lock(self.manager)
        self.holdQueue = HoldQueue(name+'_hQ', patch)
        self.holdQueue.lock(self.manager)
        self.wardDict = {}  # wards by tier

    def addWard(self, ward):
        if ward.tier not in self.wardDict:
            self.wardDict[ward.tier] = []
        self.wardDict[ward.tier].append((ward, ward.nFree))
        ward.fac = self
        return ward


class BedRequest(patches.Agent):
    STATE_START = 0
    STATE_ASKWARD = 1
    STATE_DENIEDWARD = 2
    STATE_GOTWARD = 3
    STATE_FAILED = 5
    STATE_MOVING = 6

    def __init__(self, name, patch, tier, homeWardAddr, patientKey,
                 facilityOptions, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.homeWardAddr = homeWardAddr
        self.patientKey = patientKey
        self.facilityOptions = facilityOptions
        self.tier = tier
        self.bedWard = None
        self.fsmstate = BedRequest.STATE_START
        self.dest = None

    def run(self, startTime):
        timeNow = startTime
        while True:
            if self.fsmstate == BedRequest.STATE_START:
                if len(self.facilityOptions) == 0:
                    self.fsmstate = BedRequest.STATE_FAILED
                    self.dest = None
                else:
                    self.fsmstate = BedRequest.STATE_MOVING
                    self.dest = self.facilityOptions.pop()
            elif self.fsmstate == BedRequest.STATE_MOVING:
                    addr, final = self.patch.getPathTo(self.dest)
                    if final:
                        self.fsmstate = BedRequest.STATE_ASKWARD
                    timeNow = addr.lock(self)
            elif self.fsmstate == BedRequest.STATE_ASKWARD:
                print '%s: I SHOULD BE ASLEEP at time %s' % (self.name, timeNow)
            elif self.fsmstate == BedRequest.STATE_GOTWARD:
                addr, final = self.patch.getPathTo(self.homeWardAddr)
                if final:
                    oldWard = addr
                    patientAgent = oldWard.fac.holdQueue.awaken(self.patientKey)
                    patientAgent.newWardAddr = self.bedWard
                    break
                timeNow = addr.lock(self)
            elif self.fsmstate == BedRequest.STATE_DENIEDWARD:
                self.fsmstate = BedRequest.STATE_START
            elif self.fsmstate == BedRequest.STATE_FAILED:
                addr, final = self.patch.getPathTo(self.homeWardAddr)
                if final:
                    oldWard = addr
                    patientAgent = oldWard.fac.holdQueue.awaken(self.patientKey)
                    patientAgent.newWardAddr = None
                    break
                timeNow = addr.lock(self)

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['homeWardAddr'] = self.homeWardAddr
        d['facilityOptions'] = self.facilityOptions
        d['tier'] = self.tier
        d['fsmstate'] = self.fsmstate
        d['bedWard'] = self.bedWard
        d['patientKey'] = self.patientKey
        d['dest'] = self.dest
        return d

    def __setstate__(self, stateDict):
        patches.Agent.__setstate__(self, stateDict)
        self.homeWardAddr = stateDict['homeWardAddr']
        self.facilityOptions = stateDict['facilityOptions']
        self.tier = stateDict['tier']
        self.fsmstate = stateDict['fsmstate']
        self.bedWard = stateDict['bedWard']
        self.patientKey = stateDict['patientKey']
        self.dest = stateDict['dest']


class DepartureMsg(patches.Agent):
    STATE_MOVING = 0
    STATE_ARRIVED = 1

    def __init__(self, name, patch, tier, wardAddr, destAddr, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.tier = tier
        self.wardAddr = wardAddr
        self.destAddr = destAddr
        self.fsmstate = DepartureMsg.STATE_MOVING

    def run(self, startTime):
        timeNow = startTime
        while True:
            if self.fsmstate == DepartureMsg.STATE_MOVING:
                    addr, final = self.patch.getPathTo(self.destAddr)
                    if final:
                        self.fsmstate = DepartureMsg.STATE_ARRIVED
                    timeNow = addr.lock(self)
            elif self.fsmstate == DepartureMsg.STATE_ARRIVED:
                break  # we are done

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['tier'] = self.tier
        d['fsmstate'] = self.fsmstate
        d['wardAddr'] = self.wardAddr
        d['destAddr'] = self.destAddr
        return d

    def __setstate__(self, stateDict):
        patches.Agent.__setstate__(self, stateDict)
        self.tier = stateDict['tier']
        self.fsmstate = stateDict['fsmstate']
        self.wardAddr = stateDict['wardAddr']
        self.destAddr = stateDict['destAddr']


class PatientAgent(patches.Agent):
    STATE_ATWARD = 0
    STATE_HEALED = 1
    STATE_MOVING = 2

    def __init__(self, name, patch, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.careTier = 0  # start healthy
        self.ward = None
        self.newWardAddr = None
        self.fsmstate = PatientAgent.STATE_ATWARD

    def run(self, startTime):
        timeNow = startTime
        while True:
            if self.fsmstate == PatientAgent.STATE_ATWARD:
                assert self.ward is not None, \
                    "%s: I should have been assigned to a ward" % self.name
                if self.careTier == 0:
                    timeNow = self.sleep(randint(30, 60))
                    self.careTier = 2
                elif self.careTier == 1:
                    timeNow = self.sleep(randint(14, 21))
                    self.careTier = 0
                elif self.careTier == 2:
                    timeNow = self.sleep(randint(2, 6))
                    self.careTier = 1
                else:
                    raise RuntimeError('%s: unknown care tier %s' % (self.name, self.careTier))
                self.fsmstate = PatientAgent.STATE_HEALED
            elif self.fsmstate == PatientAgent.STATE_HEALED:
                assert self.ward is not None, \
                    ("%s: I should have been assigned to a ward" % self.name)
                if self.ward.tier != self.careTier:
                    print ('%s wants a tier %d ward at %s' %
                           (self.name, self.careTier, timeNow))
                    facAddrList = [tpl[1] for tpl in self.patch.serviceLookup('BedRequestQueue')]
                    shuffle(facAddrList)
                    key = self.ward.fac.holdQueue.getUniqueKey()
                    self.patch.launch(BedRequest(self.name + '_bedReq', self.patch,
                                                 self.careTier, self.ward.getGblAddr(), key,
                                                 facAddrList),
                                      timeNow)
                    timeNow = self.ward.fac.holdQueue.lock(self, key=key)
                    # print '%s is awake with new addr %s!' % (self.name, self.newWardAddr)
                    if self.newWardAddr is None:
                        # Nowhere to go; try again tomorrow
                        print '%s is stuck here; going back to sleep' % self.name
                        timeNow = self.sleep(1)
                        self.fsmstate = PatientAgent.STATE_HEALED
                    else:
                        self.patch.launch(DepartureMsg(self.name + '_depMsg', self.patch,
                                                       self.ward.tier,
                                                       self.ward.getGblAddr(),
                                                       self.ward.fac.reqQueue.getGblAddr()),
                                          timeNow)
                        timeNow = self.ward.unlock(self)
                        self.fsmstate = PatientAgent.STATE_MOVING
                else:
                    self.fsmstate = PatientAgent.STATE_ATWARD
            elif self.fsmstate == PatientAgent.STATE_MOVING:
                self.ward = None
                addr, final = self.patch.getPathTo(self.newWardAddr)
                if final:
                    self.fsmstate = PatientAgent.STATE_ATWARD
                    self.ward = addr
                    print '%s arrived at new ward %s' % (self.name, addr._name)
                timeNow = addr.lock(self)

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['careTier'] = self.careTier
        d['ward'] = self.ward
        d['newWardAddr'] = self.newWardAddr
        d['fsmstate'] = self.fsmstate
        return d

    def __setstate__(self, stateDict):
        patches.Agent.__setstate__(self, stateDict)
        self.careTier = stateDict['careTier']
        self.ward = stateDict['ward']
        self.newWardAddr = stateDict['newWardAddr']
        self.fsmstate = stateDict['fsmstate']


class TestPatch(patches.Patch):
    def serviceLookup(self, typeNameStr):
        return self.group.worldInteractants[typeNameStr][:]

    def isLocal(self, gblAddr):
        """Is the address local to this patch?"""
        return (gblAddr.getPatchAddr() == self.tag)

    def getPathTo(self, gblAddr):
        if self.isLocal(gblAddr):
            for itr in self.interactants:
                if itr.getGblAddr() == gblAddr:
                    return (itr, True)
            raise RuntimeError("%s: Unknown supposedly-local address %s" % (self.name, gblAddr))
        else:
            patchAddr = gblAddr.getPatchAddr()
            for g in self.outgoingGates:
                if g.destTag == patchAddr:
                    return (g, False)
            raise RuntimeError("%s: No path to right patch for address %s" % (self.name, gblAddr))


class TestPatchGroup(patches.PatchGroup):

    def shareInteractantDirectories(self):
        myInteractants = defaultdict(list)
        for p in self.patches:
            pId = p.patchId
            for iact in p.interactants:
                nm = iact._name
                classNm = iact.__class__.__name__
                myInteractants[classNm].append((nm, self.getGblAddr((pId, iact.id))))
        gblAllInteractants = defaultdict(list)
        for d in self.nI.comm.allgather(myInteractants):
            for k, v in d.items():
                gblAllInteractants[k].extend(v)
        return gblAllInteractants

    def isLocal(self, gblAddr):
        return self.nI.isLocal(gblAddr)

    def start(self):
        self.worldInteractants = self.shareInteractantDirectories()
        patches.PatchGroup.start(self)


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
    patchGroup = TestPatchGroup(comm, trace=trace)
    nPatches = 2
    for j in xrange(nPatches):  # @UnusedVariable
        patch = TestPatch(patchGroup)
        patchGroup.addPatch(patch)
        facility = Facility('Facility_%s' % patch.tag, patch)
        ward0 = facility.addWard(Ward('Ward_%s_Tier0' % patch.tag, patch, 0, 1000))
        ward1 = facility.addWard(Ward('Ward_%s_Tier1' % patch.tag, patch, 1, 100))
        ward2 = facility.addWard(Ward('Ward_%s_Tier2' % patch.tag, patch, 2, 20))
        iList = [facility.reqQueue, facility.holdQueue, ward0, ward1, ward2]
        allAgents = [facility.manager]

        # Fully connect everything
        for r in xrange(comm.size):
            for jj in xrange(nPatches):
                friend = GblAddr(r, jj)
                if patch.tag != friend:  # whoops, that's me!
                    patch.addGateTo(friend)
                    patch.addGateFrom(friend)

        for i in xrange(10):
            a = PatientAgent('PatientAgent_%s_%d' % (patch.tag, i),
                             patch, debug=debug)
            ward0.lock(a)
            a.ward = ward0
            allAgents.append(a)

        patch.addInteractants(iList)
        patch.addAgents(allAgents)
    patchGroup.barrier()
    patchGroup.start()
    print '%s all done (from main)' % patchGroup.name


############
# Main hook
############

if __name__ == "__main__":
    main()

