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

import logging
import math
import types
from random import random
from collections import defaultdict
from scipy.stats import expon
import pyrheautils
import pathogenutils as pthu
from phacsl.utils.collections.phacollections import SingletonMetaClass
from stats import CachedCDFGenerator, BayesTree, fullCRVFromPDFModel
from facilitybase import CareTier, PatientOverallHealth, DiagClassA, TreatmentProtocol
from pathogenbase import Pathogen, PthStatus, defaultPthStatus
from facilitybase import PatientStatusSetter, PthStatusSetter
from genericCommunity import FreezerError

"""
This pathogen supports colonization but is never transmitted and has no consequences.  It
is useful to fill the requirements for a pathogen without having the impact of an actual
pathogen on patient flow.
"""


pathogenName = 'DUMMYPATHOGEN'
_constants_values = '$(CONSTANTS)/dummy_pathogen_constants.yaml'
_constants_schema = 'dummy_pathogen_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)

def _parseFracByTierByFacilityByCategory(fieldStr):
    return pthu.parseFracByTierByFacilityByCategory(fieldStr, _constants)

def _parseFracByTierByCategory(fieldStr, eltKey='frac'):
    return pthu.parseFracByTierByCategory(fieldStr, _constants, eltKey='frac')

def _parseScaleByTierByCategory(fieldStr):
    return pthu.parseFracByTierByCategory(fieldStr, _constants, eltKey='scale')

def _parseTierTierScaleList(fieldStr):
    return pthu.parseTierTierScaleList(fieldStr, _constants)


class DUMMYPathCore(object):
    """This is where we put things that are best shared across all CRE instances"""
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        self.initialFracColonizedTbl = _parseFracByTierByFacilityByCategory('initialFractionColonized')
        self.categoryInitialFracColonizedTbl = \
            _parseFracByTierByCategory('categoryInitialFractionColonized')
        self.colDischDelayTbl = pthu.parseValByTier('colonizedDischargeDelayTime', _constants)


    def _getInitialFracColonized(self, abbrev, category, tier):
        tierStr = CareTier.names[tier]
        if (category in self.initialFracColonizedTbl
                and abbrev in self.initialFracColonizedTbl[category]
                and tierStr in self.initialFracColonizedTbl[category][abbrev]):
            initialFracColonized = (self.initialFracColonizedTbl[category]
                                         [abbrev][tierStr])
        elif (category in self.categoryInitialFracColonizedTbl
              and tierStr in self.categoryInitialFracColonizedTbl[category]):
            initialFracColonized = (self.categoryInitialFracColonizedTbl[category][tierStr])
        else:
            raise RuntimeError('No way to set initialFracColonized for %s tier %s' %
                               (abbrev, CareTier.names[tier]))

        return initialFracColonized

    def flushCaches(self):
        """Force cached info to be regenerated"""
        pass

class DUMMYPATHOGEN(Pathogen):
    def __init__(self, ward, implCategory):
        """
        implCategory will typically be the same as that listed in the facility
        description file for the ward, but it may differ if category mapping has
        occurred.  For example, ward.fac.category might be 'SNF' while implCategory
        is 'NURSINGHOME' because that is the category definition in the class which
        actually implements 'SNF'.  This is intended to support category name mapping
        between facility descriptions and implementations.
        """
        super(DUMMYPATHOGEN, self).__init__(ward, implCategory)
        self.core = DUMMYPathCore()
        self.propogationInfo = {}
        self.propogationInfoKey = None
        self.propogationInfoTime = None
        self.treatmentProbModifierDict = None

        self.clearColonizedStatusProb = 0.0
        self.initialFracColonized = self._getInitialFracColonized(ward.fac.abbrev,
                                                                  ward.fac.category,
                                                                  ward.tier)

        self.colDischDelayTime = pthu.getValByTier(self.core.colDischDelayTbl,
                                                   'colonizedDischargeDelayTime',
                                                   ward, default=0.0)
        self.infDischDelayTime = 0.0  # the default

    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        pass

    def _getInitialFracColonized(self, abbrev, category, tier):
        return self.core._getInitialFracColonized(abbrev, category, tier)

    @classmethod
    def getCore(cls):
        if DUMMYPathCore in SingletonMetaClass._instances:
            return SingletonMetaClass._instances[DUMMYPathCore]
        else:
            dummyCore = DUMMYPathCore()  # Which will persist because it is a singleton
            return dummyCore

    @classmethod
    def getRelativeProb(cls, pthStatus, fromTier, toTier):
        """
        If the probability of transfer from fromTier to toTier of a patient at
        PthStatus.CLEAR is P, and the probability for a patient at the given PthStatus
        is kP, this routine returns the value of k.  Note that kP must still be less than 1.0,
        so there is an implied upper bound of k of 1.0/P.
        """
        return 1.0


    @classmethod
    def getEstimatedPrevalence(cls, pthStatus, abbrev, category, tier):
        """
        The return value provides an estimate prevalence (0.0 <= return val <= 1.0)
        of the pathogen for the given pathogen status at the facility named by abbrev,
        of the given category, at the given care tier.  One use of this value is to
        help keep patient flows in the correct range while rescaling the flow of a
        particular category of patient in response to colonization, etc.
        """
        if pthStatus == PthStatus.CLEAR:
            return 1.0 - cls.getCore()._getInitialFracColonized(abbrev, category, tier)
        elif pthStatus == PthStatus.COLONIZED:
            return cls.getCore()._getInitialFracColonized(abbrev, category, tier)
        else:
            return 0.0

    def _emptyPatientPth(self):
        return {k:0 for k in PthStatus.names.keys()}

    def __str__(self):
        return '<%s>' % pathogenName

    def initializePatientState(self, patient):
        """
        This method assigns the patient a Status appropriate to time 0- that is,
        it implements the initial seeding of the patient population with pathogen.
        """
        pthStatus = (PthStatus.COLONIZED if (random() <= self.initialFracColonized)
                     else defaultPthStatus)
        canClear = False
        patient._status = patient._status._replace(pthStatus=pthStatus)._replace(canClear=canClear)

    def getPatientStateKey(self, status, treatment):
        """ treatment is the patient's current treatment status """
        stateL = ['-' if status.pthStatus in [PthStatus.CLEAR, PthStatus.RECOVERED] else '+']
        for tP in self.ward.fac.treatmentPolicies:
            if hasattr(type(tP), 'treatmentKey'):
                treatmentKey = getattr(type(tP), 'treatmentKey')
                stateL.append('+' if treatment._asdict()[treatmentKey] else '-')
        return ''.join(stateL)

    def getPropogationInfo(self, timeNow):
        """Maintain and return a cached table of patient counts by status and treatment"""
        if self.propogationInfoTime != timeNow:
            pI = defaultdict(lambda: 0)
            try:
                for pt in self.ward.getPatientList():
                    pI[self.getPatientStateKey(pt.getStatus(), pt.getTreatmentProtocol())] += 1
            except FreezerError:
                pass  # We're going to have to ignore freeze-dried patients
            pI = {key: val for key, val in pI.items()
                  if not key.startswith('-')} # only sick people spread
            lst = pI.items()[:]
            lst.sort()
            self.propogationInfo = pI
            self.propogationInfoKey = tuple(lst)
            self.propogationInfoTime = timeNow
        return self.propogationInfo, self.propogationInfoKey

    def getTreatmentProbModifierDict(self):
        """Maintain and return a cached table of efficacies of combinations of treatments"""
        if not self.treatmentProbModifierDict:
            dct = {'+': (1.0, 1.0),
                   '-': (0.0, 1.0)}  # Store the two coefficients in from-to order
            for tP in self.ward.fac.treatmentPolicies:
                if hasattr(type(tP), 'treatmentKey'):
                    treatmentKey = getattr(type(tP), 'treatmentKey')
                    fromFac = tP.getTransmissionFromMultiplier(self.ward.tier,
                                                               **{treatmentKey: True})
                    toFac = tP.getTransmissionToMultiplier(self.ward.tier,
                                                           **{treatmentKey: True})
                    newDct = {}
                    for kStr, (fromProb, toProb) in dct.items():
                        newDct[kStr + '-'] = (fromProb, toProb)
                        newDct[kStr + '+'] = (fromFac * fromProb, toFac*toProb)
                    dct = newDct
            self.treatmentProbModifierDict = dct
        return self.treatmentProbModifierDict

    def getStatusChangeTree(self, patientAgent, modifierDct, startTime, timeNow):
        return BayesTree(PatientStatusSetter())  # Nothing changes.

    def filterStatusChangeTrees(self, treeList, patientAgent, startTime, timeNow):
        # Find and edit any trees containing the 'LOS' tag
        patientStatus = patientAgent.getStatus()
        newTreeList = []
        if patientStatus.pthStatus in [PthStatus.CLEAR, PthStatus.UNDETCOLONIZED]:
            key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        elif patientStatus.pthStatus in [PthStatus.COLONIZED, PthStatus.CHRONIC]:
            delayedS = max(startTime - (patientStatus.startDateA + self.colDischDelayTime), 0.0)
            delayedE = delayedS + (timeNow - startTime)
            key = (delayedS, delayedE)
        elif patientStatus.pthStatus == PthStatus.INFECTED:
            delayedS = max(startTime - (patientStatus.startDateA + self.infDischDelayTime), 0.0)
            delayedE = delayedS + (timeNow - startTime)
            key = (delayedS, delayedE)
        else:
            raise RuntimeError('patient has unexpected PthStatus %s'
                               % PthStatus.names[patientStatus.pthStatus])

        for tree in treeList:
            losTree = tree.findTag('LOS')
            if losTree:
                lTree, rTree, notReallyProb, topTag = losTree.getParts()  # @UnusedVariable
                if isinstance(notReallyProb, (types.FunctionType, types.MethodType)):
                    prob = notReallyProb(*key)
                    assert isinstance(prob, types.FloatType), 'LOS function did not return a float'
                    replacementTree = BayesTree(lTree, rTree, prob, tag="DELAYEDLOS")
                    newTreeList.append(tree.replaceSubtree('LOS', replacementTree))
                else:
                    if self.colDischDelayTime != 0.0:
                        logger.warn('%s filterStatusChangeTrees found LOS that was already evaluated',
                                    pathogenName)
                    newTreeList.append(tree)
            else:
                newTreeList.append(tree)
        return newTreeList


def getPathogenClass():
    return DUMMYPATHOGEN


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(_constants_values,
                                         _constants_schema)
