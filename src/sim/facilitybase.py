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

from random import randint, random, shuffle
import patches
import pyrheabase
from pyrheautils import enum, namedtuple


CareTier = enum('HOME', 'NURSING', 'HOSP', 'ICU')

DiagClassA = enum('HEALTHY', 'POORHEALTH', 'REHAB', 'SICK', 'VERYSICK', 'DEATH')

DiagClassB = enum('CLEAR', 'COLONIZED', 'INFECTED')

TreatmentProtocol = enum('NORMAL', 'REHAB')  # for things like patient isolation

PatientStatus = namedtuple('PatientStatus',
                           ['diagClassA',           # one of DiagClassA
                            'startDateA',           # date diagClassA status was entered
                            'diagClassB',           # one of DiagClassB
                            'startDateB'            # date diagClassB status was entered
                            ],
                           field_types=[DiagClassA, None, DiagClassB, None])

PatientDiagnosis = namedtuple('PatientDiagnosis',
                              ['diagClassA',           # one of DiagClassA
                               'diagClassB',           # one of DiagClassB
                               ],
                              field_types=[DiagClassA, DiagClassB])


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
            return (CareTier.HOME, TreatmentProtocol.NORMAL)
        if patientDiagnosis.diagClassA == DiagClassA.POORHEALTH:
            return (CareTier.NURSING, TreatmentProtocol.NORMAL)
        if patientDiagnosis.diagClassA == DiagClassA.REHAB:
            return (CareTier.NURSING, TreatmentProtocol.REHAB)
        elif patientDiagnosis.diagClassA == DiagClassA.SICK:
            return (CareTier.HOSP, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.VERYSICK:
            return (CareTier.ICU, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.DEATH:
            return (None, TreatmentProtocol.NORMAL)
        else:
            raise RuntimeError('Unknown DiagClassA %s' % str(patientDiagnosis.diagClassA))

    def diagnosisFromCareTier(self, careTier):
        """
        If I look at a patient under a given care tier, what would I expect their diagnostic class
        to be?
        """
        return {CareTier.HOME: PatientDiagnosis(DiagClassA.HEALTHY, DiagClassB.CLEAR),
                CareTier.NURSING: PatientDiagnosis(DiagClassA.POORHEALTH, DiagClassB.CLEAR),
                CareTier.HOSP: PatientDiagnosis(DiagClassA.SICK, DiagClassB.CLEAR),
                CareTier.ICU: PatientDiagnosis(DiagClassA.VERYSICK, DiagClassB.CLEAR),
                }[careTier]

    def getStatusChangeCDF(self, patientStatus, careTier, treatment, startTime, timeNow):
        """
        Return a list of pairs of the form [(prob1, newStatus), (prob2, newStatus), ...]
        where newStatus is a full patientStatus and prob is the probability of reaching
        that status at the end of timeNow.  The sum of all the probs must equal 1.0.
        """
        return [(1.0, patientStatus)]


PatientState = enum('ATWARD', 'MOVING', 'JUSTARRIVED')


class PatientAgent(patches.Agent):

    diseaseTransitionPDFs = {}

    def __init__(self, name, patch, ward, timeNow=0, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.ward = ward
        self.tier = ward.tier
        self.newWardAddr = None
        self.fsmstate = PatientState.ATWARD
        self._diagnosis = self.ward.fac.diagnosisFromCareTier(self.ward.tier)
        self._status = PatientStatus(self._diagnosis.diagClassA, 0,
                                     self._diagnosis.diagClassB, 0)
        newTier, self._treatment = self.ward.fac.prescribe(self._diagnosis,  # @UnusedVariable
                                                           TreatmentProtocol.NORMAL)
        self.lastUpdateTime = timeNow

    def printSummary(self):
        print '%s as of %s' % (self.name, self.lastUpdateTime)
        print '  status: %s ' % str(self._status)
        print '  diagnosis: %s' % str(self._diagnosis)
        print '  treatment: %s' % TreatmentProtocol.names[self._treatment]

    def updateDiseaseState(self, treatment, facility, timeNow):
        """This should embody healing, community-acquired infection, etc."""
        dT = timeNow - self.lastUpdateTime
        if dT > 0:  # moving from one ward to another can trigger two updates the same day
            probOutcomes = self.ward.fac.getStatusChangeCDF(self._status, self.tier,
                                                            self._treatment,
                                                            self.lastUpdateTime, timeNow)
            r = random()
            for deltaProb, newStatus in probOutcomes:
                r -= deltaProb
                if r < 0.0:
                    if self._status != newStatus:
                        print '%s new status %s' % (self.name, newStatus)
                    self._status = newStatus
                    break
            else:
                raise RuntimeError('Facility %s produced an incomplete PDF' % self.ward.fac.name)

    def updateEverything(self, timeNow):
        self.updateDiseaseState(self._treatment, self.ward.fac, timeNow)
        if self._status.diagClassA == DiagClassA.DEATH:
            return None
        self._diagnosis = self.ward.fac.diagnose(self._status)
        newTier, self._treatment = self.ward.fac.prescribe(self._diagnosis, self._treatment)
        self.lastUpdateTime = timeNow
        return newTier

    def run(self, startTime):
        timeNow = startTime
        # For locations with long checkintervals, random sleep to desynchronize
        if self.fsmstate == PatientState.ATWARD and self.ward.checkInterval > 1:
            timeNow = self.sleep(randint(0, self.ward.checkInterval-1))
        while True:
            if self.fsmstate == PatientState.ATWARD:
                tier = self.updateEverything(timeNow)
                if tier is None:
                    assert self._status.diagClassA == DiagClassA.DEATH, \
                        '%s is not dead?' % self.name
                    print 'Alas poor %s! %s' % (self.name, timeNow)
                    self.patch.launch(pyrheabase.DepartureMsg(self.name + '_depMsg',
                                                              self.patch, self.ward.tier,
                                                              self.ward.getGblAddr(),
                                                              self.ward.fac.reqQueue.getGblAddr()),
                                      timeNow)
                    timeNow = self.ward.unlock(self)
                    break
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
                    self.fsmstate = PatientState.JUSTARRIVED
                    self.ward = addr
                    self.tier = self.ward.tier
                    print '%s arrived at new ward %s at %s' % (self.name, addr._name, timeNow)
                timeNow = addr.lock(self)
            elif self.fsmstate == PatientState.JUSTARRIVED:
                self.fsmstate = PatientState.ATWARD
                timeNow = self.sleep(randint(0, self.ward.checkInterval-1))

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
