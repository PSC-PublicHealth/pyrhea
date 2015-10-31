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

import logging
from random import shuffle
import patches

logger = logging.getLogger(__name__)


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

    def lock(self, lockingAgent, key=None, debug=False):
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
        self.logger = logger.getChild('FacilityManager')

    def findAvailableBed(self, tier):
        if tier in self.fac.wardTierDict:
            for w in self.fac.wardTierDict[tier]:
                addr = w.getGblAddr().getLclAddr()
                ww, n = self.fac.wardAddrDict[addr]  # @UnusedVariable
                if n > 0:
                    self.fac.wardAddrDict[addr] = (w, n-1)
                    # print 'alloc: ward %s now has %d free beds of %d' % (w._name, n-1, w._nLocks)
                    return w
            return None
        else:
            return None

    def handleBedRequest(self, req, logDebug, timeNow):
        payload = req.payload
        tier = req.tier
        ward = self.findAvailableBed(tier)
        req.payload = self.fac.handleBedRequestResponse(ward, payload, timeNow)
        if ward is None:
            if logDebug:
                if tier in self.fac.wardTierDict:
                    self.logger.debug('%s: no beds available for %s' % (self.name, req.name))
                else:
                    self.logger.debug('%s: no such ward for %s' % (self.name, req.name))
            req.fsmstate = BedRequest.STATE_DENIEDWARD
        else:
            req.fsmstate = BedRequest.STATE_GOTWARD
            req.bedWard = ward.getGblAddr()
            if logDebug:
                self.logger.debug('%s: found a bed for %s' % (self.name, req.name))

    def handleArrivalMsg(self, req, logDebug, timeNow):
        payload = req.payload
        addr = req.wardAddr.getLclAddr()
        ward, nFree = self.fac.wardAddrDict[addr]  # @UnusedVariable
        self.fac.handleWardArrival(ward, payload, timeNow)
        if logDebug:
            self.logger.debug('%s: arrival at ward %s' % (self.name, ward._name))

    def handleDepartureMsg(self, req, logDebug, timeNow):
        payload = req.payload
        addr = req.wardAddr.getLclAddr()
        try:
            ward, n = self.fac.wardAddrDict[addr]
        except:
            raise RuntimeError("%s: I do not own the ward %s" % (self.name, ward._name))
        self.fac.wardAddrDict[addr] = (ward, n+1)
        # print 'depart: ward %s now has %d free beds of %d' % (ward._name, n+1, ward._nLocks)
        self.fac.handleWardDeparture(ward, payload, timeNow)
        if logDebug:
            self.logger.debug('%s: incremented ward %s' % (self.name, ward._name))

    def run(self, startTime):
        timeNow = startTime  # @UnusedVariable
        logDebug = self.logger.isEnabledFor(logging.DEBUG)
        while True:
            while self.fac.reqQueue._lockQueue:
                req = self.fac.reqQueue._lockQueue[0]
                if isinstance(req, BedRequest):
                    self.handleBedRequest(req, logDebug, timeNow)
                elif isinstance(req, ArrivalMsg):
                    self.handleArrivalMsg(req, logDebug, timeNow)
                elif isinstance(req, DepartureMsg):
                    self.handleDepartureMsg(req, logDebug, timeNow)
                else:
                    raise RuntimeError("%s unexpectedly got the message %s" %
                                       (self.name, req.name))
                self.fac.reqQueue.awaken(req)
            timeNow = self.sleep(0)  # @UnusedVariable


class Facility(object):
    def __init__(self, name, patch, managerClass=None):
        if managerClass is None:
            managerClass = FacilityManager
        self.name = name
        self.manager = managerClass(name + '_Mgr', patch, self)
        self.reqQueue = BedRequestQueue(name+'_rQ', patch)
        self.reqQueue.lock(self.manager)
        self.holdQueue = HoldQueue(name+'_hQ', patch)
        self.holdQueue.lock(self.manager)
        self.wardTierDict = {}  # lists of wards by tier
        self.wardAddrDict = {}  # dict of (ward, nFree) by lclAddr

    def addWard(self, ward):
        if ward.tier not in self.wardTierDict:
            self.wardTierDict[ward.tier] = []
        self.wardTierDict[ward.tier].append(ward)
        self.wardAddrDict[ward.getGblAddr().getLclAddr()] = (ward, ward.nFree)
        ward.fac = self
        return ward

    def getWards(self, tier=None):
        if tier:
            if tier in self.wardTierDict:
                return self.wardTierDict[tier]
            else:
                return []
        else:
            return [w for w, n in self.wardAddrDict.values()]  # @UnusedVariable

    def getTiers(self):
        return self.wardTierDict.keys()

    def getArrivalMsgPayload(self, patientAgent):
        return None

    def handleWardArrival(self, ward, payload, timeNow):
        """
        This method is called in the time slice of the facility manager when the manager
        learns of the patient's arrival.
        """
        pass

    def getDepartureMsgPayload(self, patientAgent):
        return None

    def handleWardDeparture(self, ward, payload, timeNow):
        """
        This routine is called in the time slice of the Facility Manager when the manager
        learns of the patient's departure.
        """
        pass

    def getBedRequestPayload(self, patientAgent, desiredTier):
        return None

    def handleBedRequestResponse(self, ward, payload, timeNow):
        """
        This routine is called in the time slice of the Facility Manager when the manager
        responds to a request for a bed.  If the request was denied, 'ward' will be None.
        The return value of this message becomes the new payload.
        """
        return None

    def handleBedRequestFate(self, ward, payload, timeNow):
        """
        This routine is called in the time slice of the Facility Manager when the manager
        responds to a request for a bed.  If the search for a bed failed, 'ward' will be None.
        """
        pass


class BedRequest(patches.Agent):
    STATE_START = 0
    STATE_ASKWARD = 1
    STATE_DENIEDWARD = 2
    STATE_GOTWARD = 3
    STATE_FAILED = 5
    STATE_MOVING = 6

    def __init__(self, name, patch, tier, homeWardAddr, patientKey,
                 facilityOptions, payload, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.tier = tier                  # the needed CareTier
        self.homeWardAddr = homeWardAddr  # to find our way home at the end of the search
        self.patientKey = patientKey      # to awaken the originating PatientAgent
        self.facilityOptions = facilityOptions  # candidate facilities, in preference order
        self.payload = payload            # useful for derived classes
        self.bedWard = None               # the ward satisfying the request
        self.dest = None                  # the current travel destination
        self.fsmstate = BedRequest.STATE_START

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
                raise RuntimeError('%s: I SHOULD BE ASLEEP at time %s' % (self.name, timeNow))
            elif self.fsmstate == BedRequest.STATE_GOTWARD:
                addr, final = self.patch.getPathTo(self.homeWardAddr)
                if final:
                    oldWard = addr
                    oldWard.fac.handleBedRequestFate(self.bedWard, self.payload, timeNow)
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
                    oldWard.fac.handleBedRequestFate(None, self.payload, timeNow)
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
        d['payload'] = self.payload
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
        self.payload = stateDict['payload']


class ArrivalMsg(patches.Agent):
    STATE_MOVING = 0
    STATE_ARRIVED = 1

    def __init__(self, name, patch, payload, wardAddr, destAddr, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.payload = payload
        self.wardAddr = wardAddr
        self.destAddr = destAddr
        self.fsmstate = ArrivalMsg.STATE_MOVING

    def run(self, startTime):
        timeNow = startTime  # @UnusedVariable
        while True:
            if self.fsmstate == ArrivalMsg.STATE_MOVING:
                addr, final = self.patch.getPathTo(self.destAddr)
                if final:
                    self.fsmstate = ArrivalMsg.STATE_ARRIVED
                timeNow = addr.lock(self)  # @UnusedVariable
            elif self.fsmstate == ArrivalMsg.STATE_ARRIVED:
                break  # we are done

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['payload'] = self.payload
        d['fsmstate'] = self.fsmstate
        d['wardAddr'] = self.wardAddr
        d['destAddr'] = self.destAddr
        return d

    def __setstate__(self, stateDict):
        patches.Agent.__setstate__(self, stateDict)
        self.payload = stateDict['payload']
        self.fsmstate = stateDict['fsmstate']
        self.wardAddr = stateDict['wardAddr']
        self.destAddr = stateDict['destAddr']


class DepartureMsg(patches.Agent):
    STATE_MOVING = 0
    STATE_ARRIVED = 1

    def __init__(self, name, patch, payload, wardAddr, destAddr, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.payload = payload
        self.wardAddr = wardAddr
        self.destAddr = destAddr
        self.fsmstate = DepartureMsg.STATE_MOVING

    def run(self, startTime):
        timeNow = startTime  # @UnusedVariable
        while True:
            if self.fsmstate == DepartureMsg.STATE_MOVING:
                addr, final = self.patch.getPathTo(self.destAddr)
                if final:
                    self.fsmstate = DepartureMsg.STATE_ARRIVED
                timeNow = addr.lock(self)  # @UnusedVariable
            elif self.fsmstate == DepartureMsg.STATE_ARRIVED:
                break  # we are done

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['payload'] = self.payload
        d['fsmstate'] = self.fsmstate
        d['wardAddr'] = self.wardAddr
        d['destAddr'] = self.destAddr
        return d

    def __setstate__(self, stateDict):
        patches.Agent.__setstate__(self, stateDict)
        self.payload = stateDict['payload']
        self.fsmstate = stateDict['fsmstate']
        self.wardAddr = stateDict['wardAddr']
        self.destAddr = stateDict['destAddr']


class PatientAgent(patches.Agent):
    STATE_ATWARD = 0
    STATE_MOVING = 1
    STATE_JUSTARRIVED = 2

    def __init__(self, name, patch, ward, timeNow=0, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.ward = ward
        self.tier = ward.tier
        self.newWardAddr = None
        self.fsmstate = PatientAgent.STATE_ATWARD
        self.logger = logging.getLogger(__name__ + '.PatientAgent')

    def getPostArrivalPauseTime(self, timeNow):
        return 0

    def handleTierUpdate(self, timeNow):
        return self.ward.tier

    def handlePatientDeath(self, timeNow):
        pass

    def getCandidateFacilityList(self, timeNow, newTier):
        facAddrList = [tpl[1] for tpl in self.patch.serviceLookup('BedRequestQueue')]
        shuffle(facAddrList)
        return facAddrList

    def run(self, startTime):
        timeNow = startTime
        timeNow = self.sleep(self.getPostArrivalPauseTime(timeNow))
        while True:
            if self.debug:
                self.logger.debug('%s point 0: status %s state %s day %s' %
                                  (self.name, str(self._status), self.fsmstate, timeNow))
            if self.fsmstate == PatientAgent.STATE_ATWARD:
                tier = self.handleTierUpdate(timeNow)
                if self.debug:
                    self.logger.debug('%s point 1: status %s tiers %s vs %s' %
                                      (self.name, str(self._status), tier, self.ward.tier))
                if tier is None:
                    self.handlePatientDeath(timeNow)
                    self.patch.launch(DepartureMsg(self.name + '_depMsg',
                                                   self.patch,
                                                   self.ward.fac.getDepartureMsgPayload(self),
                                                   self.ward.getGblAddr(),
                                                   self.ward.fac.reqQueue.getGblAddr()),
                                      timeNow)
                    timeNow = self.ward.unlock(self)
                    break
                elif self.ward.tier == tier:
                    if self.debug:
                        self.logger.debug('%s point 9: status %s' % (self.name, str(self._status)))
                    timeNow = self.sleep(self.ward.checkInterval)
                else:
                    # print ('%s wants a tier %s ward at %s' %
                    #        (self.name, CareTier.names[tier], timeNow))
                    facAddrList = self.getCandidateFacilityList(timeNow, tier)
                    key = self.ward.fac.holdQueue.getUniqueKey()
                    self.patch.launch(BedRequest(self.name + '_bedReq', self.patch,
                                                 tier, self.ward.getGblAddr(),
                                                 key, facAddrList,
                                                 self.ward.fac.getBedRequestPayload(self, tier)),
                                      timeNow)
                    if self.debug:
                        self.logger.debug('%s point 3: launched req for tier %s'
                                          % (self.name, tier))
                    timeNow = self.ward.fac.holdQueue.lock(self, key=key)
                    if self.debug:
                        self.logger.debug('%s point 4: status %s day %s'
                                          % (self.name, str(self._status), timeNow))
                    if self.newWardAddr is None:
                        # Nowhere to go; try again tomorrow
                        # print ('%s is stuck at %s; going back to sleep at %s' %
                        #        (self.name, self.ward, timeNow))
                        if self.debug:
                            self.logger.debug('%s point 8: status %s day %s'
                                              % (self.name, str(self._status), timeNow))
                        timeNow = self.sleep(1)
                        # state is unchanged
                    else:
                        if self.debug:
                            self.logger.debug('%s point 5: status %s day %s'
                                              % (self.name, str(self._status), timeNow))
                        rQAddr = self.ward.fac.reqQueue.getGblAddr()
                        self.patch.launch(DepartureMsg(self.name + '_depMsg',
                                                       self.patch,
                                                       self.ward.fac.getDepartureMsgPayload(self),
                                                       self.ward.getGblAddr(),
                                                       rQAddr),
                                          timeNow)
                        timeNow = self.ward.unlock(self)
                        self.fsmstate = PatientAgent.STATE_MOVING
            elif self.fsmstate == PatientAgent.STATE_MOVING:
                self.ward = None
                addr, final = self.patch.getPathTo(self.newWardAddr)
                if final:
                    self.fsmstate = PatientAgent.STATE_JUSTARRIVED
                    self.ward = addr
                    self.tier = self.ward.tier
                    if self.debug:
                        self.logger.debug('%s point 6: arrive tier %s day %s'
                                          % (self.name, self.ward.tier, timeNow))
                timeNow = addr.lock(self)
            elif self.fsmstate == PatientAgent.STATE_JUSTARRIVED:
                rQAddr = self.ward.fac.reqQueue.getGblAddr()
                self.patch.launch(ArrivalMsg(self.name + '_arvMsg',
                                             self.patch,
                                             self.ward.fac.getArrivalMsgPayload(self),
                                             self.ward.getGblAddr(),
                                             rQAddr),
                                  timeNow)
                self.fsmstate = PatientAgent.STATE_ATWARD
                if self.debug:
                    self.logger.debug('%s point 7: status %s day %s'
                                      % (self.name, str(self._status), timeNow))
                timeNow = self.sleep(self.getPostArrivalPauseTime(timeNow))
            if self.debug:
                self.logger.debug('%s point 2: status %s'
                                  % (self.name, str(self._status)))

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['tier'] = self.tier
        d['ward'] = self.ward
        d['newWardAddr'] = self.newWardAddr
        d['fsmstate'] = self.fsmstate
        return d

    def __setstate__(self, d):
        patches.Agent.__setstate__(self, d)
        self.tier = d['tier']
        self.ward = d['ward']
        self.newWardAddr = d['newWardAddr']
        self.fsmstate = d['fsmstate']
