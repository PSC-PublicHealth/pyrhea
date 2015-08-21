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

import math
from random import randint, random, shuffle
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
        if dT > 0:  # moving from one ward to another can trigger two updates the same day
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
                    if self.newWardAddr is None:
                        # Nowhere to go; try again tomorrow
                        print ('%s is stuck at %s; going back to sleep at %s' %
                               (self.name, self.ward, timeNow))
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
