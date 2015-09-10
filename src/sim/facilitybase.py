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

import types
from random import randint, random, shuffle
import patches
import pyrheabase
from pyrheautils import enum, namedtuple
from math import fabs


CareTier = enum('HOME', 'NURSING', 'HOSP', 'ICU')

PatientOverallHealth = enum('HEALTHY', 'FRAIL')

DiagClassA = enum('HEALTHY', 'NEEDSREHAB', 'SICK', 'VERYSICK', 'DEATH')

DiagClassB = enum('CLEAR', 'COLONIZED', 'INFECTED')

TreatmentProtocol = enum('NORMAL', 'REHAB')  # for things like patient isolation

PatientStatus = namedtuple('PatientStatus',
                           ['overall',              # one of patientOverallHealth
                            'diagClassA',           # one of DiagClassA
                            'startDateA',           # date diagClassA status was entered
                            'diagClassB',           # one of DiagClassB
                            'startDateB'            # date diagClassB status was entered
                            ],
                           field_types=[PatientOverallHealth, DiagClassA, None, DiagClassB, None])

PatientDiagnosis = namedtuple('PatientDiagnosis',
                              ['overall',              # one of PatientOverallHealth
                               'diagClassA',           # one of DiagClassA
                               'diagClassB',           # one of DiagClassB
                               ],
                              field_types=[PatientOverallHealth, DiagClassA, DiagClassB])


class CachedCDFGenerator:
    """
    This supports a common operation performed by Facilities- given a treatment interval (relative
    to a start date of zero), return a likelihood.  For example, this may be the likelihood that
    the patient will be discharged.  Since the interval bounds are integers, caching of the
    generated values is very effective.
    """
    def __init__(self, frozenPDF):
        self.frozenPDF = frozenPDF
        self.cache = {}

    def intervalProb(self, start, end):
        key = (start, end)
        if key in self.cache:
            return self.cache[key]
        else:
            cP = self.frozenPDF.cdf(end) - self.frozenPDF.cdf(start)
            self.cache[key] = cP
            return cP


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
        return PatientDiagnosis(patientStatus.overall,
                                patientStatus.diagClassA, patientStatus.diagClassB)

    def prescribe(self, patientDiagnosis, patientTreatment):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            return (CareTier.HOME, TreatmentProtocol.NORMAL)
        if patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
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
        to be?  This is used for patient initialization purposes.
        """
        if careTier == CareTier.HOME:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.HEALTHY, DiagClassB.CLEAR)
        elif careTier == CareTier.NURSING:
            return PatientDiagnosis(PatientOverallHealth.FRAIL,
                                    DiagClassA.HEALTHY, DiagClassB.CLEAR)
        elif careTier == CareTier.HOSP:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.SICK, DiagClassB.CLEAR)
        elif careTier == CareTier.ICU:
            return PatientDiagnosis(PatientOverallHealth.HEALTHY,
                                    DiagClassA.VERYSICK, DiagClassB.CLEAR)
        else:
            raise RuntimeError('Unknown care tier %s' % careTier)

    def getStatusChangeTree(self, patientStatus, careTier, treatment, startTime, timeNow):
        """
        Return a Bayes tree the traversal of which yields a patientStatus.

        Elements of the tree are either a terminal patientStatus value or a 3-tuple
        of the form (k opt1 opt2), where opt1 and opt2 are likewise elements of the tree
        and k is a float in the range 0 to 1. The tree is evaluated by generating a random
        number each time a 3-tuple is encountered and choosing the first or second opt
        if the random value is <= or > k respectively.  When a value is encountered which
        is not a 3-tuple, that value is returned as the result of the traversal.
        """
        return patientStatus

    @staticmethod
    def foldCDF(linearCDF):
        """
        Given a linear CDF, return an equivalent CDF in Bayes Tree form

        The input CDF is of the form:

           [(p1, r1), (p2, r2), ...]

        such that the sum of p1, p2, etc is 1.0 .  The output CDF is of the form:

            pdf = (tp1 pdf1 pdf2)

        or

            pdf = r1

        where tp1 is a float between 0 and 1, and the second form is equivalent to (1.0 r1 None)
        """
        tol = 0.00001
        lCDF = len(linearCDF)
        if lCDF == 1:
            p, r = linearCDF[0]
            assert fabs(p - 1.0) <= tol, 'CDF terms do not sum to 1 (%s)' % linearCDF
            return r
        elif lCDF == 2:
            p1, r1 = linearCDF[0]
            p2, r2 = linearCDF[1]
            assert fabs(p1 + p2 - 1.0) <= tol, 'CDF terms do not sum to 1 (%s)' % linearCDF
            if p2 >= p1:
                return [p1, r1, r2]
            else:
                return [p2, r2, r1]
        else:
            linearCDF.sort()
            part1 = linearCDF[:lCDF/2]
            part2 = linearCDF[lCDF/2:]
            pivot = part1[-1][0]
            if pivot == 0.0:
                return Facility.foldCDF(part2)
            else:
                w1 = sum([p for p, r in part1])
                w2 = sum([p for p, r in part2])
                return [pivot,
                        Facility.foldCDF([(p / w1, r) for p, r in part1]),
                        Facility.foldCDF([(p / w2, r) for p, r in part2])
                        ]

PatientState = enum('ATWARD', 'MOVING', 'JUSTARRIVED')


class PatientAgent(patches.Agent):

    def __init__(self, name, patch, ward, timeNow=0, debug=False):
        patches.Agent.__init__(self, name, patch, debug=debug)
        self.ward = ward
        self.tier = ward.tier
        self.newWardAddr = None
        self.fsmstate = PatientState.ATWARD
        self._diagnosis = self.ward.fac.diagnosisFromCareTier(self.ward.tier)
        self._status = PatientStatus(PatientOverallHealth.HEALTHY,
                                     self._diagnosis.diagClassA, 0,
                                     self._diagnosis.diagClassB, 0)
        newTier, self._treatment = self.ward.fac.prescribe(self._diagnosis,  # @UnusedVariable
                                                           TreatmentProtocol.NORMAL)
        self.lastUpdateTime = timeNow

    def printSummary(self):
        print '%s as of %s' % (self.name, self.lastUpdateTime)
        print '  status: %s ' % str(self._status)
        print '  diagnosis: %s' % str(self._diagnosis)
        print '  treatment: %s' % TreatmentProtocol.names[self._treatment]

    @staticmethod
    def _traverseBayesTree(tree):
        # print 'walk tree: %s' % str(tree)
        while True:
            if isinstance(tree, types.ListType):
                if random() <= tree[0]:
                    tree = tree[1]
                else:
                    tree = tree[2]
            else:
                # print 'walk returning %s' % str(tree)
                return tree

    def updateDiseaseState(self, treatment, facility, timeNow):
        """This should embody healing, community-acquired infection, etc."""
        dT = timeNow - self.lastUpdateTime
        if dT > 0:  # moving from one ward to another can trigger two updates the same day
            tree = self.ward.fac.getStatusChangeTree(self._status, self.tier, self._treatment,
                                                     self.lastUpdateTime, timeNow)
            self._status = self._traverseBayesTree(tree)

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
                    # print 'Alas poor %s! %s' % (self.name, timeNow)
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
                    # print ('%s wants a tier %s ward at %s' %
                    #        (self.name, CareTier.names[tier], timeNow))
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
                        # print ('%s is stuck at %s; going back to sleep at %s' %
                        #        (self.name, self.ward, timeNow))
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
                    # print '%s arrived at new ward %s at %s' % (self.name, addr._name, timeNow)
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
