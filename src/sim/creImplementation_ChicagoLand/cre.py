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
from phacsl.utils.collections.phacollections import SingletonMetaClass
from stats import CachedCDFGenerator, BayesTree, fullCRVFromPDFModel
from facilitybase import CareTier, PatientOverallHealth, DiagClassA, TreatmentProtocol
from pathogenbase import Pathogen, PthStatus, defaultPthStatus
from facilitybase import PatientStatusSetter, PthStatusSetter

pathogenName = 'CRE'
_constants_values = '$(CONSTANTS)/cre_constants.yaml'
_constants_schema = 'cre_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)

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
            assert tier not in topD[cat], ('Redundant categoryInitialFracColonized for %s %s' %
                                           (cat, CareTier.names[tier]))
            topD[cat][tier] = val
    return topD

def _parseScaleByTierByCategory(fieldStr):
    return _parseFracByTierByCategory(fieldStr, eltKey='scale')

def _parseTierTierScaleList(fieldStr):
    result = {}
    for elt in _constants[fieldStr]:
        fmTier = CareTier.__dict__[elt['tierFrom']]
        toTier = CareTier.__dict__[elt['tierTo']]
        value = elt['scale']['value']
        result[(fmTier, toTier)] = value
    return result

def _getValByTierByCategory(tbl, tblNameForErr, ward, implCategory, overrideTbl=None):
    wardCat = ward.fac.category
    tierStr = CareTier.names[ward.tier]
    if overrideTbl:
        # Potential override values are stored by [category][abbrev][tierName]
        abbrev = ward.fac.abbrev
        if (wardCat in overrideTbl
            and abbrev in overrideTbl[wardCat]
            and tierStr in overrideTbl[wardCat][abbrev]):
            return overrideTbl[wardCat][abbrev][tierStr]
    if (wardCat in tbl and tierStr in tbl[wardCat]):
        return tbl[wardCat][tierStr]
    else:
        raise RuntimeError('No way to set %s for %s tier %s' %
                           (tblNameForErr, ward.fac.abbrev, CareTier.names[ward.tier]))


class CRECore(object):
    """This is where we put things that are best shared across all CRE instances"""
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        self.initialFracColonizedTbl = _parseFracByTierByFacilityByCategory('initialFractionColonized')
        self.categoryInitialFracColonizedTbl = \
            _parseFracByTierByCategory('categoryInitialFractionColonized')
        self.tauTbl = _parseFracByTierByCategory('tau')
        self.tauOverrideTbl = _parseFracByTierByFacilityByCategory('facilityTau')
        self.exposureCutoffTbl = _parseScaleByTierByCategory('exposureCutoff')
        self.fracPermanentlyColonized = _constants['fracPermanentlyColonized']['value']
        self.spontaneousLossTimeConstant = _constants['spontaneousLossTimeConstant']['value']
        self.transferProbScaleDict = _parseTierTierScaleList('colonizedTransferProbScale')
        self.exposureTreeCache = {}
        self.spontaneousLossTreeCache = {}
        self.spontaneousLossCachedCDF = CachedCDFGenerator(expon(scale=self.spontaneousLossTimeConstant))
        
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
        self.spontaneousLossTreeCache = {}
        self.exposureTreeCache = {}

class CRE(Pathogen):
    def __init__(self, ward, implCategory):
        """
        implCategory will typically be the same as that listed in the facility
        description file for the ward, but it may differ if category mapping has
        occurred.  For example, ward.fac.category might be 'SNF' while implCategory
        is 'NURSINGHOME' because that is the category definition in the class which
        actually implements 'SNF'.  This is intended to support category name mapping
        between facility descriptions and implementations.
        """
        super(CRE, self).__init__(ward, implCategory)
        self.core = CRECore()
        self.patientPth = self._emptyPatientPth()
        self.patientPthTime = None
        self.propogationInfo = {}
        self.propogationInfoKey = None
        self.propogationInfoTime = None
        self.treatmentProbModifierDict = None

        self.initialFracColonized = self._getInitialFracColonized(ward.fac.abbrev,
                                                                  ward.fac.category,
                                                                  ward.tier)

        self.tau = _getValByTierByCategory(self.core.tauTbl, 'tau', ward, implCategory,
                                           overrideTbl=self.core.tauOverrideTbl)
#         if 'PRES_100_L' in ward._name:
#             print '########## tau is %s' % self.tau
        self.exposureCutoff = _getValByTierByCategory(self.core.exposureCutoffTbl,
                                                      'exposureCutoff',
                                                      ward, implCategory)

        tierName = CareTier.names[self.ward.tier]
        self.colDischDelayTime = 0.0  # the default
        for ent in _constants['colonizedDischargeDelayTime']:
            if ent['tier'] == tierName:
                self.colDischDelayTime = ent['value']
                break

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

    def _getInitialFracColonized(self, abbrev, category, tier):
        return self.core._getInitialFracColonized(abbrev, category, tier)

    @classmethod
    def getCore(cls):
        if CRECore in SingletonMetaClass._instances:
            return SingletonMetaClass._instances[CRECore]
        else:
            dummyCore = CRECore()  # Which will persist because it is a singleton
            return dummyCore

    @classmethod
    def getRelativeProb(cls, pthStatus, fromTier, toTier):
        """
        If the probability of transfer from fromTier to toTier of a patient at
        PthStatus.CLEAR is P, and the probability for a patient at the given PthStatus
        is kP, this routine returns the value of k.  Note that kP must still be less than 1.0,
        so there is an implied upper bound of k of 1.0/P.
        """
        if pthStatus == PthStatus.COLONIZED:
            key = (fromTier, toTier)
            tD = cls.getCore().transferProbScaleDict
            if key in tD:
                return tD[key]
            else:
                return 1.0
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
        canClear = (random() > self.core.fracPermanentlyColonized)
        patient._status = patient._status._replace(pthStatus=pthStatus)._replace(canClear=canClear)

    def getPatientPthCounts(self, timeNow):
        """
        Returns a dict of the form {PthStatus.CLEAR : nClearPatients, ...}
        """
        if self.patientPthTime != timeNow:
            d = self._emptyPatientPth()
            for p in self.ward.getPatientList():
                d[p._status.pthStatus] += 1
            self.patientPth = d
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

    def getStatusChangeTree(self, patientStatus, ward, treatment, startTime, timeNow):
        if patientStatus.pthStatus == PthStatus.CLEAR:
            pI, pIKey = self.getPropogationInfo(timeNow)
            dT = timeNow - startTime
            #
            # Note that this caching scheme assumes all wards with the same category and tier
            # have the same infectivity constant(s)
            #
            key = (self.ward.fac.category, self.ward.tier, treatment, pIKey, dT, self.tau)
#             if 'THC_4058_L' in ward._name:
#                 print '%s %s %s %s' % (ward._name, patientStatus.justArrived, pI, treatment)
            if key not in self.core.exposureTreeCache:

                patientKey = self.getPatientStateKey(patientStatus, treatment)
                tPMD = self.getTreatmentProbModifierDict()
                selfProt = tPMD[patientKey][1]
                logPSafe = 0.0
                totPop = float(sum(pI.values()))
                expScale = 1.0
                if self.exposureCutoff < totPop:
                    expScale = (self.exposureCutoff) / (totPop)
                for k2, ct in pI.items():
                    logPSafe += math.log(1.0 - (tPMD[k2][0] * selfProt * self.tau)) * ct
#                 if 'PRES_100_L' in ward._name:
#                     print '%s %s %s %s %s %s' % (ward._name, patientKey, logPSafe, tPMD,
#                                                  selfProt, self.tau)
                pSafe = math.exp(dT * expScale * logPSafe)
                tree = BayesTree(PatientStatusSetter(),
                                 PthStatusSetter(PthStatus.COLONIZED),                                 
                                 pSafe)
                self.core.exposureTreeCache[key] = tree
            return self.core.exposureTreeCache[key]
        else:
            assert patientStatus.pthStatus == PthStatus.COLONIZED, ('patient has unexpected PthStatus %s' %
                                                                    PthStatus.names[patientStatus.pthStatus])
            key = (startTime - patientStatus.startDatePth, timeNow - patientStatus.startDatePth)
            if patientStatus.canClear:
                if key not in self.core.spontaneousLossTreeCache:
                    changeProb = self.core.spontaneousLossCachedCDF.intervalProb(*key)
                    tree = BayesTree(PthStatusSetter(PthStatus.CLEAR),
                                     PatientStatusSetter(),
                                     changeProb)
                    self.core.spontaneousLossTreeCache[key] = tree
                return self.core.spontaneousLossTreeCache[key]
            else:
                return BayesTree(PatientStatusSetter())

    def filterStatusChangeTrees(self, treeList, patientStatus, careTier, treatment, startTime, timeNow):
        # Find and edit any trees containing the 'LOS' tag
        newTreeList = []
        if patientStatus.pthStatus == PthStatus.CLEAR:
            key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        elif patientStatus.pthStatus == PthStatus.COLONIZED:
            delayedS = max(startTime - (patientStatus.startDateA + self.colDischDelayTime), 0.0)
            delayedE = delayedS + (timeNow - startTime)
            key = (delayedS, delayedE)
        else:
            raise RuntimeError('patient has unexpected PthStatus %s' % PthStatus.names[patientStatus.pthStatus])
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
    return CRE


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(_constants_values,
                                         _constants_schema)
