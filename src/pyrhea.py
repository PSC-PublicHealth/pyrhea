#! /usr/bin/env python

_rhea_svn_id_ = "$Id$"

import sys
import math
from random import randint, random, shuffle, seed, choice
import patches
import pyrheabase
from pyrheautils import enum, namedtuple


CareTier = enum('HOME', 'REHAB', 'ASSISTED', 'NURSING', 'HOSP', 'ICU')

DiagClassA = enum('HEALTHY', 'SICK')

DiagClassB = enum('CLEAR', 'COLONIZED', 'INFECTED')

TreatmentProtocol = enum('NORMAL')  # for things like patient isolation

PatientStatus = namedtuple('PatientStatus',
                           ['diagClassA',           # one of DiagClassA
                            'diagACountdown',       # days until better
                            'diagClassB',           # one of DiagClassB
                            ],
                           field_types=[DiagClassA, None, DiagClassB])

PatientDiagnosis = namedtuple('PatientDiagnosis',
                              ['diagClassA',           # one of DiagClassA
                               'diagClassB',           # one of DiagClassB
                               ],
                              field_types=[DiagClassA, DiagClassB])

healthyGetSickProbPerDay = 0.01


class Ward(pyrheabase.Ward):
    def __init__(self, name, patch, tier, nBeds):
        pyrheabase.Ward.__init__(self, name, patch, tier, nBeds)
        self.checkInterval = 1  # check health daily


class CommunityWard(Ward):
    """This 'ward' type represents being out in the community"""

    def __init__(self, name, patch, nBeds):
        Ward.__init__(self, name, patch, CareTier.HOME, nBeds)
        self.checkInterval = 30  # check health daily


class Facility(pyrheabase.Facility):
    def diagnose(self, patientStatus):
        """
        This provides a way to introduce false positive or false negative diagnoses.  The
        only way in which patient status affects treatment policy or ward is via diagnosis.
        """
        return PatientDiagnosis(patientStatus.diagClassA, patientStatus.diagClassB)

    def prescribe(self, patientDiagnosis, patientTreatment):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            if patientDiagnosis.diagClassB in [DiagClassB.CLEAR, DiagClassB.COLONIZED]:
                return (CareTier.HOME, TreatmentProtocol.NORMAL)
            elif patientDiagnosis.diagClassB == DiagClassB.INFECTED:
                return (CareTier.HOSP, TreatmentProtocol.NORMAL)
            else:
                raise RuntimeError('Unknown DiagClassB %s' % str(patientDiagnosis.diagClassB))
        elif patientDiagnosis.diagClassA == DiagClassA.SICK:
            if patientDiagnosis.diagClassB in [DiagClassB.CLEAR, DiagClassB.COLONIZED]:
                return (CareTier.HOSP, TreatmentProtocol.NORMAL)
            elif patientDiagnosis.diagClassB == DiagClassB.INFECTED:
                return (CareTier.ICU, TreatmentProtocol.NORMAL)
            else:
                raise RuntimeError('Unknown DiagClassB %s' % str(patientDiagnosis.diagClassB))
        else:
            raise RuntimeError('Unknown DiagClassA %s' % str(patientDiagnosis.diagClassA))


PatientState = enum('ATWARD', 'MOVING')


class PatientAgent(patches.Agent):

    def __init__(self, name, patch, ward, timeNow=0, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.ward = ward
        self.tier = ward.tier
        self.newWardAddr = None
        self.fsmstate = PatientState.ATWARD
        self._status = PatientStatus(DiagClassA.HEALTHY, 0, DiagClassB.CLEAR)
        self._diagnosis = self.ward.fac.diagnose(self._status)
        newTier, self._treatment = self.ward.fac.prescribe(self._diagnosis,  # @UnusedVariable
                                                           TreatmentProtocol.NORMAL)
        self.lastUpdateTime = timeNow

    def printSummary(self):
        print '%s as of %s' % (self.name, self.lastUpdateTime)
        print '  status: %s ' % str(self._status)
        print '  diagnosis: %s' % str(self._diagnosis)
        print '  treatment: %s' % TreatmentProtocol.names[self._treatment]

    def updateDiseaseState(self, treatment, timeNow):
        """This should embody healing, community-acquired infection, etc."""
        dT = timeNow - self.lastUpdateTime
        self.lastUpdateTime = timeNow
        diagA = self._status.diagClassA
        diagACountdown = self._status.diagACountdown
        diagB = self._status.diagClassB
        if self._status.diagClassA == DiagClassA.HEALTHY:
            pBar = 1.0 - healthyGetSickProbPerDay
            gotSick = (random() > math.pow(pBar, dT))
            if gotSick:
                diagA = DiagClassA.SICK
                diagACountdown = randint(3, 10)
            else:
                pass  # no update; still healthy
        elif self._status.diagClassA == DiagClassA.SICK:
            diagACountdown -= dT
            if diagACountdown <= 0:
                diagA = DiagClassA.HEALTHY
                diagACountdown = 0
        else:
            raise RuntimeError('Unknown DiagClassA %s' % str(diagA))

        if diagB == DiagClassB.CLEAR:
            if random() < 0.1:
                diagB = DiagClassB.COLONIZED
        elif diagB == DiagClassB.COLONIZED:
            r = random()
            if r < 0.3:
                diagB = DiagClassB.CLEAR
            elif r < 0.7:
                diagB = DiagClassB.COLONIZED
            else:
                diagB = DiagClassB.INFECTED
        elif diagB == DiagClassB.INFECTED:
            if random() < 0.2:
                diagB = DiagClassB.COLONIZED
            else:
                diagB = DiagClassB.INFECTED
        else:
            raise RuntimeError('Unknown DiagClassB %s' % str(diagB))

        self._status = PatientStatus(diagA, diagACountdown, diagB)

    def updateEverything(self, timeNow):
        self.updateDiseaseState(self._treatment, timeNow)
        self._diagnosis = self.ward.fac.diagnose(self._status)
        newTier, self._treatment = self.ward.fac.prescribe(self._diagnosis, self._treatment)
        self.lastUpdateTime = timeNow
        return newTier

    def run(self, startTime):
        timeNow = startTime
        while True:
            if self.fsmstate == PatientState.ATWARD:
                tier = self.updateEverything(timeNow)
                if self.ward.tier == tier:
                    timeNow = self.sleep(self.ward.checkInterval)
                else:
                    print ('%s wants a tier %s ward at %s' %
                           (self.name, CareTier.names[tier], timeNow))
                    facAddrList = [tpl[1] for tpl in self.patch.serviceLookup('BedRequestQueue')]
                    shuffle(facAddrList)
                    key = self.ward.fac.holdQueue.getUniqueKey()
                    self.patch.launch(pyrheabase.BedRequest(self.name + '_bedReq', self.patch,
                                                            tier, self.ward.getGblAddr(),
                                                            key, facAddrList),
                                      timeNow)
                    timeNow = self.ward.fac.holdQueue.lock(self, key=key)
                    # print '%s is awake with new addr %s!' % (self.name, self.newWardAddr)
                    if self.newWardAddr is None:
                        # Nowhere to go; try again tomorrow
                        print '%s is stuck here; going back to sleep' % self.name
                        timeNow = self.sleep(1)
                        # state is unchanged
                    else:
                        rQAddr = self.ward.fac.reqQueue.getGblAddr()
                        self.patch.launch(pyrheabase.DepartureMsg(self.name + '_depMsg',
                                                                  self.patch, self.ward.tier,
                                                                  self.ward.getGblAddr(),
                                                                  rQAddr),
                                          timeNow)
                        timeNow = self.ward.unlock(self)
                        self.fsmstate = PatientState.MOVING
            elif self.fsmstate == PatientState.MOVING:
                self.ward = None
                addr, final = self.patch.getPathTo(self.newWardAddr)
                if final:
                    self.fsmstate = PatientState.ATWARD
                    self.ward = addr
                    print '%s arrived at new ward %s' % (self.name, addr._name)
                timeNow = addr.lock(self)

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['tier'] = self.tier
        d['ward'] = self.ward
        d['newWardAddr'] = self.newWardAddr
        d['fsmstate'] = self.fsmstate
        d['status'] = self._status
        d['diagnosis'] = self._diagnosis
        d['treatment'] = self._treatment
        d['lastUpdateTime'] = self.lastUpdateTime
        return d

    def __setstate__(self, d):
        patches.Agent.__setstate__(self, d)
        self.tier = d['tier']
        self.ward = d['ward']
        self.newWardAddr = d['newWardAddr']
        self.fsmstate = d['fsmstate']
        self._status = d['status']
        self._diagnosis = d['diagnosis']
        self._treatment = d['treatment']
        self.lastUpdateTime = d['lastUpdateTime']


def describeSelf():
    print """This should write some documentation"""


def main():
    trace = False
    verbose = False  # @UnusedVariable
    debug = False
    deterministic = False

    for a in sys.argv[1:]:
        if a == '-v':
            verbose = True  # @UnusedVariable
        elif a == '-d':
            debug = True
        elif a == '-t':
            trace = True
        elif a == '-D':
            deterministic = True
        else:
            describeSelf()
            sys.exit('unrecognized argument %s' % a)

    comm = patches.getCommWorld()

    if deterministic:
        seed(1234 + comm.rank)  # Set the random number generator seed

    patchGroup = patches.PatchGroup(comm, trace=trace, deterministic=deterministic)
    nPatches = 2
    for j in xrange(nPatches):  # @UnusedVariable
        patch = patchGroup.addPatch(patches.Patch(patchGroup))

        hosp = Facility('Hospital_%s' % str(patch.tag), patch)
        rehab = Facility('Rehab_%s' % str(patch.tag), patch)
        assisted = Facility('Assisted_%s' % str(patch.tag), patch)
        community = Facility('Community_%s' % str(patch.tag), patch)

        allItr = [hosp.reqQueue, hosp.holdQueue,
                  rehab.reqQueue, rehab.holdQueue,
                  assisted.reqQueue, assisted.holdQueue,
                  community.reqQueue, community.holdQueue]
        allAgents = [hosp.manager, rehab.manager, assisted.manager, community.manager]

        count = 0
        for tier, nBeds in [(CareTier.HOSP, 100), (CareTier.HOSP, 50), (CareTier.ICU, 20)]:
            allItr.append(hosp.addWard(Ward('Ward_%s_%s_%d' % (hosp.name,
                                                               CareTier.names[tier],
                                                               count),
                                            patch, tier, nBeds)))
            count += 1

        count = 0
        for tier, nBeds in [(CareTier.REHAB, 200)]:
            allItr.append(rehab.addWard(Ward('Ward_%s_%s_%d' % (rehab.name,
                                                                CareTier.names[tier],
                                                                count),
                                             patch, tier, nBeds)))
            count += 1

        count = 0
        for tier, nBeds in [(CareTier.NURSING, 100), (CareTier.ASSISTED, 1000)]:
            allItr.append(assisted.addWard(Ward('Ward_%s_%s_%d' % (assisted.name,
                                                                   CareTier.names[tier],
                                                                   count),
                                                patch, tier, nBeds)))
            count += 1

        count = 0
        communityWards = []
        for tier, nBeds in [(CareTier.HOME, 5500), (CareTier.HOME, 5500)]:
            w = community.addWard(CommunityWard('Community_%s_%s_%d' % (community.name,
                                                                        CareTier.names[tier],
                                                                        count),
                                                patch, nBeds))
            communityWards.append(w)
            allItr.append(w)
            count += 1

        for i in xrange(10000):
            ward = choice(communityWards)
            a = PatientAgent('PatientAgent_%s_%d' % (str(patch.tag), i),
                             patch, choice(communityWards), debug=debug)
            ward.lock(a)
            a.ward = ward
            allAgents.append(a)

        patch.addInteractants(allItr)
        patch.addAgents(allAgents)

    allAgents[-1].printSummary()

    patchGroup.start()
    print '%s all done (from main)' % patchGroup.name


############
# Main hook
############

if __name__ == "__main__":
    main()
