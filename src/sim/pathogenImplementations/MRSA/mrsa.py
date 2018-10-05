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
from random import random
import types
from collections import defaultdict
import pyrheautils
import pathogenutils as pthu
from phacsl.utils.collections.phacollections import SingletonMetaClass
from stats import CachedCDFGenerator, BayesTree, fullCRVFromPDFModel
from facilitybase import CareTier
from facilitybase import PatientStatusSetter, PthStatus, PthStatusSetter
from pathogenbase import Pathogen
from genericCommunity import FreezerError

pathogenName = 'MRSA'
_constants_values = '$(CONSTANTS)/mrsa_constants.yaml'
_constants_schema = 'mrsa_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)


def _parseValByTier(fieldStr):
    return pthu.parseValByTier(fieldStr, _constants)

def _parseFracByTierByFacilityByCategory(fieldStr):
    return pthu.parseFracByTierByFacilityByCategory(fieldStr, _constants)

def _parseFracByTierByCategory(fieldStr, eltKey='frac'):
    return pthu.parseFracByTierByCategory(fieldStr, _constants, eltKey='frac')

def _parseScaleByTierByCategory(fieldStr):
    return pthu.parseFracByTierByCategory(fieldStr, _constants, eltKey='scale')

def _parseValByTierByCategory(fieldStr):
    return pthu.parseFracByTierByCategory(fieldStr, _constants, eltKey='value')

def _parseTierTierScaleList(fieldStr):
    return pthu.parseTierTierScaleList(fieldStr, _constants)

class MRSACore(object):
    """This is where we put things that are best shared across all MRSA instances"""
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        infCRV = fullCRVFromPDFModel(_constants['infectionDurationPDF'])
        self.infCachedCDF = CachedCDFGenerator(infCRV)
        reboundCRV = fullCRVFromPDFModel(_constants['reboundPDF'])
        self.reboundCachedCDF = CachedCDFGenerator(reboundCRV)
        sponLossCRV = fullCRVFromPDFModel(_constants['spontaneousLossPDF'])
        self.sponLossCachedCDF = CachedCDFGenerator(sponLossCRV)
        colToInfectCRV = fullCRVFromPDFModel(_constants['colonizationToInfectionPDF'])
        self.colToInfectCachedCDF = CachedCDFGenerator(colToInfectCRV)

        self.exposureProbCache = {}  # so simple we don't need a CachedCDFGenerator

        self.initialFracColonizedTbl = \
            _parseFracByTierByFacilityByCategory('initialFractionColonized')
        self.categoryInitialFracColonizedTbl = \
            _parseFracByTierByCategory('categoryInitialFractionColonized')
        self.initialFracInfectedTbl = \
            _parseFracByTierByFacilityByCategory('initialFractionInfected')
        self.categoryInitialFracInfectedTbl = \
            _parseFracByTierByCategory('categoryInitialFractionInfected')
        self.initialFracUndetColTbl = \
            _parseFracByTierByFacilityByCategory('initialFractionUndetectablyColonized')
        self.categoryInitialFracUndetColTbl = \
            _parseFracByTierByCategory('categoryInitialFractionUndetectablyColonized')
        self.initialFracChronicTbl = \
            _parseFracByTierByFacilityByCategory('initialFractionChronicallyColonized')
        self.categoryInitialFracChronicTbl = \
            _parseFracByTierByCategory('categoryInitialFractionChronicallyColonized')

        self.tauTbl = _parseFracByTierByCategory('tau')
        self.tauOverrideTbl = _parseFracByTierByFacilityByCategory('facilityTau')
        self.exposureCutoffTbl = _parseScaleByTierByCategory('exposureCutoff')
        self.colDischDelayTbl = _parseValByTier('colonizedDischargeDelayTime')
        self.infDischDelayTbl = _parseValByTier('infectedDischargeDelayTime')

        self.probNewExposuresUndet = _constants['probNewExposuresAreUndetectable']['value']
        self.fracPermanentlyColonized = _constants['fracPermanentlyColonized']['value']
        self.colTransferProbScaleDict = _parseTierTierScaleList('colonizedTransferProbScale')
        self.infTransferProbScaleDict = _parseTierTierScaleList('infectedTransferProbScale')


    def _getInitialPthFrac(self, pthStatus, abbrev, category, tier):
        tierStr = CareTier.names[tier]
        if pthStatus == PthStatus.COLONIZED:
            tbl1, tbl2 = self.initialFracColonizedTbl, self.categoryInitialFracColonizedTbl
        elif pthStatus == PthStatus.INFECTED:
            tbl1, tbl2 = self.initialFracInfectedTbl, self.categoryInitialFracInfectedTbl
        elif pthStatus == PthStatus.UNDETCOLONIZED:
            tbl1, tbl2 = self.initialFracUndetColTbl, self.categoryInitialFracUndetColTbl
        elif pthStatus == PthStatus.CHRONIC:
            tbl1, tbl2 = self.initialFracChronicTbl, self.categoryInitialFracChronicTbl
        else:
            raise RuntimeError("No implementation for initial population in pathogen state %s"
                               % PthStatus.names[pthStatus])
        if category in tbl1 and abbrev in tbl1[category] and tierStr in tbl1[category][abbrev]:
            initialFrac = tbl1[category][abbrev][tierStr]
        elif category in tbl2 and tierStr in tbl2[category]:
            initialFrac = tbl2[category][tierStr]
        else:
            raise RuntimeError('No way to set initial fraction %s for %s tier %s' %
                               (PthStatus.names[pthStatus], abbrev, CareTier.names[tier]))
        return initialFrac

    def flushCaches(self):
        """Force cached info to be regenerated"""
#        self.spontaneousLossTreeCache = {}
        self.exposureProbCache = {}


class MRSA(Pathogen):
    def __init__(self, ward, implCategory):
        """
        implCategory will typically be the same as that listed in the facility
        description file for the ward, but it may differ if category mapping has
        occurred.  For example, ward.category might be 'SNF' while implCategory
        is 'NURSINGHOME'.  This is intended to support category name mapping between
        facility descriptions and implementations.
        """
        super(MRSA, self).__init__(ward, implCategory)
        self.core = MRSACore()
        self.propogationInfo = {}
        self.propogationInfoKey = None
        self.propogationInfoTime = None
        self.treatmentProbModifierDict = None

        initialFractionColonized = self.core._getInitialPthFrac(PthStatus.COLONIZED,
                                                                ward.fac.abbrev,
                                                                ward.fac.category, ward.tier)
        self.initialFracColonized = initialFractionColonized
        initialFractionInfected = self.core._getInitialPthFrac(PthStatus.INFECTED,
                                                               ward.fac.abbrev,
                                                               ward.fac.category, ward.tier)
        self.initialFracInfected = initialFractionInfected
        initialFractionUndetCol = self.core._getInitialPthFrac(PthStatus.UNDETCOLONIZED,
                                                               ward.fac.abbrev,
                                                               ward.fac.category, ward.tier)
        self.initialFracUndetCol = initialFractionUndetCol
        initialFractionChronic = self.core._getInitialPthFrac(PthStatus.CHRONIC,
                                                              ward.fac.abbrev,
                                                              ward.fac.category, ward.tier)
        self.initialFracChronic = initialFractionChronic

        fracClear = 1.0 - (initialFractionColonized + initialFractionInfected
                           + initialFractionUndetCol + initialFractionChronic)
        self.initializationBayesTree = BayesTree.fromLinearCDF([(initialFractionColonized,
                                                                 PthStatusSetter(PthStatus.COLONIZED)),
                                                                (initialFractionInfected,
                                                                 PthStatusSetter(PthStatus.INFECTED)),
                                                                (initialFractionUndetCol,
                                                                 PthStatusSetter(PthStatus.UNDETCOLONIZED)),
                                                                (initialFractionChronic,
                                                                 PthStatusSetter(PthStatus.CHRONIC)),
                                                                (fracClear,
                                                                 PatientStatusSetter())])
        self.tau = pthu.getValByTierByCategory(self.core.tauTbl, 'tau', ward, ward.fac.category,
                                               overrideTbl=self.core.tauOverrideTbl)
        self.exposureCutoff = pthu.getValByTierByCategory(self.core.exposureCutoffTbl,
                                                          'exposureCutoff',
                                                          ward, ward.fac.category)
        self.colDischDelayTime = pthu.getValByTier(self.core.colDischDelayTbl,
                                                   'colonizedDischargeDelayTime',
                                                   ward, default=0.0)
        self.infDischDelayTime = pthu.getValByTier(self.core.colDischDelayTbl,
                                                   'infectedDischargeDelayTime',
                                                   ward, default=0.0)

    def __str__(self):
        return '<%s>' % pathogenName

    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        self.core.flushCaches()
        self.treatmentProbModifierDict = None
        self.propogationInfoTime = None

    @classmethod
    def getCore(cls):
        """Access the singleton core based only on class"""
        if MRSACore in SingletonMetaClass._instances:
            return SingletonMetaClass._instances[MRSACore]
        else:
            dummyCore = MRSACore()  # Which will persist because it is a singleton
            return dummyCore

    @classmethod
    def getRelativeProb(cls, pthStatus, fromTier, toTier):
        """
        If the probability of transfer from fromTier to toTier of a patient at
        PthStatus.CLEAR is P, and the probability for a patient at the given PthStatus
        is kP, this routine returns the value of k.  Note that kP must still be less than 1.0,
        so there is an implied upper bound of k of 1.0/P.
        """
        key = (fromTier, toTier)
        if pthStatus in [PthStatus.COLONIZED, PthStatus.CHRONIC]:
            tD = cls.getCore().colTransferProbScaleDict
        elif pthStatus == PthStatus.INFECTED:
            tD = cls.getCore().infTransferProbScaleDict
        else:
            # All other states presumed to look 'healthy'
            tD = {}
        if key in tD:
            return tD[key]
        else:
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
        if pthStatus in [PthStatus.COLONIZED, PthStatus.INFECTED, PthStatus.UNDETCOLONIZED,
                         PthStatus.CHRONIC]:
            return cls.getCore()._getInitialPthFrac(pthStatus, abbrev, category, tier)
        elif pthStatus == PthStatus.CLEAR:
            return 1.0 - sum([cls.getCore()._getInitialPthFrac(pthS, abbrev, category, tier)
                              for pthS in [PthStatus.COLONIZED, PthStatus.INFECTED,
                                           PthStatus.UNDETCOLONIZED]])
        else:
            return 0.0

    def initializePatientState(self, patient):
        """
        This method assigns the patient a Status appropriate to time 0- that is,
        it implements the initial seeding of the patient population with pathogen.
        """
        patient._status = self.initializationBayesTree.traverse().set(patient.getStatus(), 0)
        canClear = (random() > self.core.fracPermanentlyColonized)
        patient._status = patient.getStatus()._replace(canClear=canClear)
        

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
                  if not key.startswith('-')} # because only sick people spread
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
        patientStatus = patientAgent.getStatus()
        ward = patientAgent.ward
        treatment = patientAgent.getTreatmentProtocol()
        if patientStatus.pthStatus in [PthStatus.CLEAR, PthStatus.UNDETCOLONIZED]:
            pI, pIKey = self.getPropogationInfo(timeNow)
            dT = timeNow - startTime
            #
            # Note that this caching scheme assumes all wards with the same category and tier
            # have the same infectivity constant(s)
            #
            key = (self.ward.fac.category, self.ward.tier, treatment, pIKey, dT, self.tau)
            if key not in self.core.exposureProbCache:

                patientKey = self.getPatientStateKey(patientStatus, treatment)
                tPMD = self.getTreatmentProbModifierDict()
                selfProt = tPMD[patientKey][1]
                logPSafe = 0.0
                totPop = float(sum(pI.values()))
                expScale = 1.0
                if self.exposureCutoff < totPop:
                    expScale = self.exposureCutoff / totPop
                for k2, ct in pI.items():
                    logPSafe += math.log(1.0 - (tPMD[k2][0] * selfProt * self.tau)) * ct
                pNotExposed = math.exp(dT * expScale * logPSafe)
                self.core.exposureProbCache[key] = pNotExposed

            # For CLEAR individuals, the only possible state change is to colonized via exposure.
            # for UNDETCOLONIZED individuals, there's a possible random transition as well.
            if patientStatus.pthStatus == PthStatus.CLEAR:
                pNotExposed = self.core.exposureProbCache[key]
                return BayesTree(PatientStatusSetter(),
                                 BayesTree(PthStatusSetter(PthStatus.UNDETCOLONIZED),
                                           PthStatusSetter(PthStatus.COLONIZED),
                                           self.core.probNewExposuresUndet),
                                 pNotExposed)
            else:
                # status is PthStatus.UNDETCOLONIZED
                timeKey = (startTime - patientStatus.startDatePth,
                           timeNow - patientStatus.startDatePth)
                pRebound = self.core.reboundCachedCDF.intervalProb(*timeKey)
                pSponLoss = self.core.sponLossCachedCDF.intervalProb(*timeKey)
                pNotExposed = self.core.exposureProbCache[key]
                pTransition = pSponLoss + pRebound - (pSponLoss * pRebound)
                # The following works if intervals are short enough to consider the pdfs constant
                pRatio = pSponLoss/(pSponLoss + pRebound)
                innerTree = BayesTree(BayesTree(PthStatusSetter(PthStatus.CLEAR),
                                                PthStatusSetter(PthStatus.COLONIZED),
                                                pRatio),
                                      PatientStatusSetter(),
                                      pTransition)
                return BayesTree(innerTree,
                                 BayesTree(PthStatusSetter(PthStatus.UNDETCOLONIZED),
                                           PthStatusSetter(PthStatus.COLONIZED),
                                           self.core.probNewExposuresUndet),
                                 pNotExposed)
        elif patientStatus.pthStatus == PthStatus.COLONIZED:
            # Infection has not been fully implemented- TODO add the colonized-to-infected transition
            timeKey = (startTime - patientStatus.startDatePth, timeNow - patientStatus.startDatePth)
            if patientStatus.canClear:
                pSponLoss = self.core.sponLossCachedCDF.intervalProb(*timeKey)
                return BayesTree(PthStatusSetter(PthStatus.CLEAR),
                                 PatientStatusSetter(),
                                 pSponLoss)
            else:
                return BayesTree(PatientStatusSetter())
        elif patientStatus.pthStatus == PthStatus.CHRONIC:
            # Chronic infection is forever
            return BayesTree(PatientStatusSetter())
        elif patientStatus.pthStatus == PthStatus.INFECTED:
            # Infection has not yet been fully implemented
            timeKey = (startTime - patientStatus.startDatePth, timeNow - patientStatus.startDatePth)
            pRecovery = self.core.infCachedCDF.intervalProb(*timeKey)
            return BayesTree(PthStatusSetter(PthStatus.CLEAR),
                             PatientStatusSetter(),
                             pRecovery)
        else:
            raise RuntimeError('patient has unexpected PthStatus %s' %
                               PthStatus.names[patientStatus.pthStatus])

    def filterStatusChangeTrees(self, treeList, patientAgent, startTime, timeNow):
        # Find and edit any trees containing the 'LOS' tag
        patientStatus = patientAgent.getStatus()
        ward = patientAgent.ward
        treatment = patientAgent.getTreatmentProtocol()
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
                        logger.warn('%s filterStatusChangeTrees: LOS was already evaluated!',
                                    pathogenName)
                    newTreeList.append(tree)
            else:
                newTreeList.append(tree)
        return newTreeList


def getPathogenClass():
    """What is the pathogen class this module provides?"""
    return MRSA


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(_constants_values,
                                         _constants_schema)
