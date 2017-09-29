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

import os.path
import logging
import math
import types
from collections import defaultdict
import pyrheautils
from phacsl.utils.collections.phacollections import SingletonMetaClass
from stats import CachedCDFGenerator, BayesTree, fullCRVFromPDFModel
from facilitybase import CareTier, PatientOverallHealth, DiagClassA, TreatmentProtocol
from facilitybase import PatientStatusSetter, PthStatusSetter
from pathogenbase import Pathogen, PthStatus, defaultPthStatus

pathogenName = 'MRSA'
_constants_values = 'mrsa_constants.yaml'
_constants_schema = 'mrsa_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)


def _valFromCategoryEntry(key, ctg, constantsJSON):
    if key not in constantsJSON:
        raise RuntimeError('Constant list for %s was not found' % key)
    for v in constantsJSON[key]:
        if v['category'] == ctg:
            return v['frac']['value']
    raise RuntimeError('Constant entry for category %s was not found' % ctg)

def _parseValByTier(fieldStr):
    topD = {}
    for elt in _constants[fieldStr]:
        topD[elt['tier']] = elt['value']
    return topD

def _parseFracByTierByFacilityByCategory(fieldStr):
    topD = {}
    for elt in _constants[fieldStr]:
        cat = elt['category']
        facilities = elt['facilities']
        if cat not in topD:
            topD[cat] = {}
        for facElt in facilities:
            abbrev = facElt['abbrev']
            tiers = facElt['tiers']
            if abbrev not in topD[cat]:
                topD[cat][abbrev] = {}
            for tierElt in tiers:
                tier = tierElt['tier']
                val = tierElt['frac']['value']
                assert tier not in topD[cat][abbrev], ('Redundant %s for %s %s' %
                                                      (fieldStr, abbrev, CareTier.names[tier]))
                topD[cat][abbrev][tier] = val
    return topD


def _parseFracByTierByCategory(fieldStr, eltKey='frac'):
    topD = {}
    for elt in _constants[fieldStr]:
        cat = elt['category']
        tiers = elt['tiers']
        if cat not in topD:
            topD[cat] = {}
        for tierElt in tiers:
            tier = tierElt['tier']
            val = tierElt[eltKey]['value']
            assert tier not in topD[cat], ('Redundant %s for %s %s' %
                                           (fieldStr, cat, CareTier.names[tier]))
            topD[cat][tier] = val
    return topD

def _parseScaleByTierByCategory(fieldStr):
    return _parseFracByTierByCategory(fieldStr, eltKey='scale')

def _parseValByTierByCategory(fieldStr):
    return _parseFracByTierByCategory(fieldStr, eltKey='value')

def _parseTierTierScaleList(fieldStr):
    result = {}
    for elt in _constants[fieldStr]:
        fmTier = CareTier.__dict__[elt['tierFrom']]
        toTier = CareTier.__dict__[elt['tierTo']]
        value = elt['scale']['value']
        result[(fmTier, toTier)] = value
    return result

def _getValByTierByCategory(tbl, tblNameForErr, ward, wardCategory, overrideTbl=None,
                            default=None):
    tierStr = CareTier.names[ward.tier]
    if overrideTbl:
        # Potential override values are stored by [category][abbrev][tierName]
        abbrev = ward.fac.abbrev
        if (wardCategory in overrideTbl
            and abbrev in overrideTbl[wardCategory]
            and tierStr in overrideTbl[wardCategory][abbrev]):
            return overrideTbl[wardCategory][abbrev][tierStr]
    if wardCategory in tbl and tierStr in tbl[wardCategory]:
        return tbl[wardCategory][tierStr]
    else:
        if default is not None:
            return default
        else:
            raise RuntimeError('No way to set %s for %s tier %s' %
                               (tblNameForErr, ward.fac.abbrev, CareTier.names[ward.tier]))


class MRSACore(object):
    """This is where we put things that are best shared across all MRSA instances"""
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        colCRV = fullCRVFromPDFModel(_constants['colonizationDurationPDF'])
        self.colCachedCDF = CachedCDFGenerator(colCRV)
        infCRV = fullCRVFromPDFModel(_constants['infectionDurationPDF'])
        self.infCachedCDF = CachedCDFGenerator(infCRV)
        reboundCRV = fullCRVFromPDFModel(_constants['reboundPDF'])
        self.reboundCachedCDF = CachedCDFGenerator(reboundCRV)
        sponLossCRV = fullCRVFromPDFModel(_constants['spontaneousLossPDF'])
        self.sponLossCachedCDF = CachedCDFGenerator(sponLossCRV)

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
        self.patientPth = self._emptyPatientPth()
        self.patientPthTime = None
        self.propogationInfo = {}
        self.propogationInfoKey = None
        self.propogationInfoTime = None
        self.treatmentProbModifierDict = None

        initialFractionColonized = self.core._getInitialPthFrac(PthStatus.COLONIZED,
                                                                ward.fac.abbrev,
                                                                ward.fac.category, ward.tier)
        initialFractionInfected = self.core._getInitialPthFrac(PthStatus.INFECTED,
                                                               ward.fac.abbrev,
                                                               ward.fac.category, ward.tier)
        initialFractionUndetCol = self.core._getInitialPthFrac(PthStatus.UNDETCOLONIZED,
                                                               ward.fac.abbrev,
                                                               ward.fac.category, ward.tier)
        initialFractionChronic = self.core._getInitialPthFrac(PthStatus.CHRONIC,
                                                              ward.fac.abbrev,
                                                              ward.fac.category, ward.tier)
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
        self.tau = _getValByTierByCategory(self.core.tauTbl, 'tau', ward, ward.fac.category,
                                           overrideTbl=self.core.tauOverrideTbl)
        self.exposureCutoff = _getValByTierByCategory(self.core.exposureCutoffTbl,
                                                      'exposureCutoff',
                                                      ward, ward.fac.category)
        self.colDischDelayTime = _getValByTierByCategory(self.core.colDischDelayTbl,
                                                         'colonizedDischargeDelayTime',
                                                         ward, ward.fac.category, default=0.0)
        self.infDischDelayTime = _getValByTierByCategory(self.core.colDischDelayTbl,
                                                         'infectedDischargeDelayTime',
                                                         ward, ward.fac.category, default=0.0)


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

    def _emptyPatientPth(self):
        return {k:0 for k in PthStatus.names.keys()}

    def __str__(self):
        return '<%s>' % pathogenName

    def initializePatientState(self, patient):
        """
        This method assigns the patient a Status appropriate to time 0- that is,
        it implements the initial seeding of the patient population with pathogen.
        """
        patient._status = self.initializationBayesTree.traverse().set(patient._status, 0)

    def getPatientPthCounts(self, timeNow):
        """
        Returns a dict of the form {PthStatus.CLEAR : nClearPatients, ...}
        """
        if self.patientPthTime != timeNow:
            dct = self._emptyPatientPth()
            for pt in self.ward.getPatientList():
                dct[pt._status.pthStatus] += 1
            self.patientPth = dct
            self.patientPthTime = timeNow
        return self.patientPth

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
            for pt in self.ward.getPatientList():
                pI[self.getPatientStateKey(pt._status, pt._treatment)] += 1
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

    def getStatusChangeTree(self, patientStatus, ward, treatment, startTime, timeNow):
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
                return BayesTree(PthStatusSetter(PthStatus.COLONIZED),
                                 BayesTree(BayesTree(PthStatusSetter(PthStatus.CLEAR),
                                                     PatientStatusSetter(),
                                                     pSponLoss),
                                           PthStatusSetter(PthStatus.COLONIZED),
                                           pNotExposed),
                                 pRebound)
        elif patientStatus.pthStatus == PthStatus.COLONIZED:
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

    def filterStatusChangeTrees(self, treeList, patientStatus, careTier, treatment,
                                startTime, timeNow):
        # Find and edit any trees containing the 'LOS' tag
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
_constants = pyrheautils.importConstants(os.path.join(os.path.dirname(__file__),
                                                      _constants_values),
                                         _constants_schema)
