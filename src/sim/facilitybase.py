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
from random import randint, random
import pyrheabase
from pyrheautils import enum, namedtuple
from math import fabs
import logging
from phacsl.utils.notes.statval import HistoVal

logger = logging.getLogger(__name__)


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


class PatientStatusSetter(object):
    def __init__(self):
        pass

    def set(self, patientStatus, timeNow):
        """This base class just returns a copy"""
        return PatientStatus._make(patientStatus)

    def __str__(self):
        return 'PatientStatusSetter()'


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
            cP = ((self.frozenPDF.cdf(end) - self.frozenPDF.cdf(start))
                  / (1.0 - self.frozenPDF.cdf(start)))
            self.cache[key] = cP
            return cP


class Ward(pyrheabase.Ward):
    def __init__(self, name, patch, tier, nBeds):
        pyrheabase.Ward.__init__(self, name, patch, tier, nBeds)
        self.checkInterval = 1  # check health daily


class FacilityManager(pyrheabase.FacilityManager):
    pass


class Facility(pyrheabase.Facility):
    class PatientRecord(object):
        def __init__(self, ward, patientID, arrivalDate):
            self.patientID = patientID
            self.arrivalDate = arrivalDate
            self.departureDate = None
            self.prevVisits = 0

    def __init__(self, name, patch):
        pyrheabase.Facility.__init__(self, name, patch, managerClass=FacilityManager)
        self.noteHolder = None
        self.idCounter = 0
        self.patientDataDict = {}

    def setNoteHolder(self, noteHolder):
        self.noteHolder = noteHolder
        self.noteHolder.addNote({'name': self.name})

    def getNoteHolder(self):
        return self.noteHolder

    def getArrivalMsgPayload(self, patientAgent):
        return patientAgent.id

    def handleWardArrival(self, ward, payload, timeNow):
        """
        The locking agent should be an arriving patient.  This allows a sort of 'check-in' for
        record-keeping purposes.
        """
        patientID = payload
        if patientID in self.patientDataDict:
            logger.warning('Patient %s has returned to %s' % (patientID, self.name))
            patientRec = self.patientDataDict[patientID]
            patientRec.prevVisits += 1
            patientRec.arrivalDate = timeNow
            patientRec.departureDate = None
        else:
            patientRec = Facility.PatientRecord(ward, patientID, timeNow)
            self.patientDataDict[patientID] = patientRec
        nh = self.getNoteHolder()
        if nh:
            nh.addNote({(CareTier.names[ward.tier] + '_arrivals'): 1})

    def getDepartureMsgPayload(self, patientAgent):
        return patientAgent.id

    def handleWardDeparture(self, ward, payload, timeNow):
        """
        This routine is called in the time slice of the Facility Manager when the manager
        learns of the patient's departure.
        """
        patientID = payload
        if patientID not in self.patientDataDict:
            logger.error('%s has no record of patient %s' % (self.name, patientID))
        patientRec = self.patientDataDict[patientID]
        patientRec.departureDate = timeNow
        lengthOfStay = timeNow - patientRec.arrivalDate
        losKey = (CareTier.names[ward.tier] + '_LOS')
        nh = self.getNoteHolder()
        if losKey not in nh:
            nh.addNote({losKey: HistoVal([])})
        nh.addNote({(CareTier.names[ward.tier] + '_departures'): 1, losKey: lengthOfStay})

    def getBedRequestPayload(self, patientAgent, desiredTier):
        return (0, desiredTier)  # number of bounces and tier

    def handleBedRequestResponse(self, ward, payload, timeNow):
        """
        This routine is called in the time slice of the Facility Manager when the manager
        responds to a request for a bed.  If the request was denied, 'ward' will be None.
        The return value of this message becomes the new payload.
        """
        nBounces, tier = payload
        return (nBounces + 1, tier)  # new number of bounces and tier

    def handleBedRequestFate(self, ward, payload, timeNow):
        """
        This routine is called in the time slice of the Facility Manager when the manager
        responds to a request for a bed.  If the search for a bed failed, 'ward' will be None.
        """
        nBounces, tier = payload
        if ward is None:
            nFail = 1
            nSuccess = 0
        else:
            nFail = 0
            nSuccess = 1
        nh = self.getNoteHolder()
        bounceKey = CareTier.names[tier] + '_bounce_histo'
        if bounceKey not in nh:
            nh.addNote({bounceKey: HistoVal([])})
        self.getNoteHolder().addNote({(CareTier.names[tier] + '_found'): nSuccess,
                                      (CareTier.names[tier] + '_notfound'): nFail,
                                      bounceKey: nBounces})

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
        return PatientStatusSetter()  # A simple setter that just makes a copy

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


class PatientAgent(pyrheabase.PatientAgent):
    idCounter = 0  # to provide a reliable identifier for each patient

    def __init__(self, name, patch, ward, timeNow=0, debug=False):
        pyrheabase.PatientAgent.__init__(self, name, patch, ward, timeNow=timeNow, debug=debug)
        self._diagnosis = self.ward.fac.diagnosisFromCareTier(self.ward.tier)
        self._status = PatientStatus(PatientOverallHealth.HEALTHY,
                                     self._diagnosis.diagClassA, 0,
                                     self._diagnosis.diagClassB, 0)
        newTier, self._treatment = self.ward.fac.prescribe(self._diagnosis,  # @UnusedVariable
                                                           TreatmentProtocol.NORMAL)
        self.lastUpdateTime = timeNow
        self.id = (patch.patchId, PatientAgent.idCounter)
        PatientAgent.idCounter += 1
        self.logger = logging.getLogger(__name__ + '.PatientAgent')

    def printSummary(self):
        print '%s as of %s' % (self.name, self.lastUpdateTime)
        print '  status: %s ' % str(self._status)
        print '  diagnosis: %s' % str(self._diagnosis)
        print '  treatment: %s' % TreatmentProtocol.names[self._treatment]

    @staticmethod
    def _printBayesTree(tree, indent=0):
        if isinstance(tree, types.ListType):
            prob, path1, path2 = tree
            print '%s(%f' % (' '*indent, prob)
            PatientAgent._printBayesTree(path1, indent+4)
            PatientAgent._printBayesTree(path2, indent+4)
            print '%s)' % (' '*indent)
        else:
            print '%s%s' % (' '*indent, tree)

    @staticmethod
    def _traverseBayesTree(tree):
        # PatientAgent._printBayesTree(tree)
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
            setter = self._traverseBayesTree(tree)
            self._status = setter.set(self._status, timeNow)

    def updateEverything(self, timeNow):
        self.updateDiseaseState(self._treatment, self.ward.fac, timeNow)
        if self._status.diagClassA == DiagClassA.DEATH:
            return None
        self._diagnosis = self.ward.fac.diagnose(self._status)
        newTier, self._treatment = self.ward.fac.prescribe(self._diagnosis, self._treatment)
        if self.debug:
            self.logger.debug('%s: status -> %s, diagnosis -> %s, treatment -> %s'
                              % (self.name, self._status, self._diagnosis, self._treatment))
        self.lastUpdateTime = timeNow
        return newTier

    def getPostArrivalPauseTime(self, timeNow):
        if self.ward.checkInterval > 1:
            return randint(0, self.ward.checkInterval-1)
        else:
            return 0

    def handleTierUpdate(self, timeNow):
        return self.updateEverything(timeNow)

    def handlePatientDeath(self, timeNow):
        assert self._status.diagClassA == DiagClassA.DEATH, '%s is not dead?' % self.name
        self.ward.fac.noteHolder.addNote({"death": 1})
        if self.debug:
            self.logger.debug('Alas poor %s! %s' % (self.name, timeNow))

    def __getstate__(self):
        d = pyrheabase.PatientAgent.__getstate__(self)
        d['status'] = self._status
        d['diagnosis'] = self._diagnosis
        d['treatment'] = self._treatment
        d['lastUpdateTime'] = self.lastUpdateTime
        d['id'] = self.id
        return d

    def __setstate__(self, d):
        pyrheabase.PatientAgent.__setstate__(self, d)
        self._status = d['status']
        self._diagnosis = d['diagnosis']
        self._treatment = d['treatment']
        self.lastUpdateTime = d['lastUpdateTime']
        self.id = d['id']
