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

import os.path
import jsonschema
import yaml
import math
from scipy.stats import lognorm
import logging
from random import shuffle

import pyrheabase
import pyrheautils
import schemautils
from facilitybase import CareTier, PthStatus, DiagClassA
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent
from facilitybase import PatientStatusSetter, LTACQueue
from stats import CachedCDFGenerator, BayesTree
from hospital import estimateWork as hospitalEstimateWork
from hospital import buildChangeTree, biasTransfers

category = 'LTAC'
_schema = 'hospitalfacts_schema.yaml'
_constants_values = '$(MODELDIR)/constants/ltac_constants.yaml'
_constants_schema = 'ltac_constants_schema.yaml'
_validator = None
_constants = None

logger = logging.getLogger(__name__)

class LTAC(Facility):
    def __init__(self, descr, patch, policyClasses=None, categoryNameMapper=None):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr,
                          descr, patch,
                          reqQueueClasses=[LTACQueue],
                          policyClasses=policyClasses,
                          categoryNameMapper=categoryNameMapper)
        descr = self.mapDescrFields(descr)
        bedsPerWard = _constants['bedsPerWard']['value']

        _c = _constants
        totDsch = float(descr['totalDischarges']['value'])
        totTO = sum([elt['count']['value'] for elt in descr['totalTransfersOut']])
        ttoD = {elt['category']: elt['count']['value'] for elt in descr['totalTransfersOut']}
        fracNotTransferred = (float(totDsch) - float(totTO)) / float(totDsch)
        assert fracNotTransferred >= 0.0, '%s has more transfers than discharges' % self.name

        self.lclRates = {}
        self.lclRates['death'] = _c['dischargeViaDeathFrac']['value']
        assert fracNotTransferred >= self.lclRates['death'], '%s has more deaths than non-transfers'
        self.lclRates['home'] = fracNotTransferred - self.lclRates['death']

        expectedKeys = ['HOSPITAL', 'LTAC', 'NURSINGHOME', 'VSNF']
        if totTO:
            for key in ttoD.keys():
                assert key in expectedKeys, \
                    '%s records outgoing transfers to unexpected category %s' % (self.name, key)
            for key in expectedKeys:
                nTransfers = ttoD[key] if key in ttoD else 0
                self.lclRates[key.lower()] = float(nTransfers) / float(totDsch)
        else:
            for key in expectedKeys:
                self.lclRates[key.lower()] = 0.0
            logger.warning('%s has no transfers out', self.name)

        for subCat, cat, key in [('icu', 'hospital', 'hospTransferToICURate'),
                                 ('vent', 'vsnf', 'vsnfTransferToVentRate'),
                                 ('skilnrs', 'vsnf', 'vsnfTransferToSkilNrsRate')]:
            self.lclRates[subCat] = _c[key]['value'] * self.lclRates[cat]
            self.lclRates[cat] -= self.lclRates[subCat]

        # VSNF-specific care tiers have been separated out, so any remaining vsnf fraction
        # is actually just nursing.
        self.lclRates['nursinghome'] += self.lclRates['vsnf']
        self.lclRates['vsnf'] = 0.0

        # pthRates is the equivalent of lclRates but for pathogen carriers.  It must
        # be initialized lazily because the pathogen has not yet been defined.
        self.pthRates = None

        if 'nBeds' in descr:
            nBeds = descr['nBeds']['value']
            nWards = int(float(nBeds)/bedsPerWard) + 1
        else:
            meanPop = descr['meanPop']['value']
            nWards = int(math.ceil(meanPop / bedsPerWard))

        nBeds = nWards * bedsPerWard

        # Add a way to dial the LOS up and down
        scaledLOSParms = self.scaleLOS(descr['losModel'],
                                       (descr['scaleLengthOfStay']['value'] if
                                        'scaleLengthOfStay' in descr else None))
        self.cachedCDF = CachedCDFGenerator(lognorm(scaledLOSParms[1],
                                                    scale=math.exp(scaledLOSParms[0])))
        self.treeCache = {}

        for i in xrange(nWards):
            self.addWard(Ward(('%s_%s_%s_%s_%d' %
                               (category, patch.name, descr['abbrev'], 'LTAC', i)),
                              patch, CareTier.LTAC, bedsPerWard))

    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        self.treeCache = {}

    def getStatusChangeTree(self, patientAgent, startTime, timeNow):
        patientStatus = patientAgent.getStatus()
        ward = patientAgent.ward
        treatment = patientAgent.getTreatmentProtocol()
        if not self.pthRates:
            # Lazy initialization
            infectiousAgent = self.getWards(CareTier.LTAC)[0].iA
            self.lclRates, self.pthRates = biasTransfers(self.lclRates,
                                                         infectiousAgent,
                                                         CareTier.ICU,
                                                         self.abbrev, self.category,
                                                         CareTier.LTAC)


        if ward.tier == CareTier.LTAC:
            if patientAgent.getDiagnosis().diagClassA == DiagClassA.NEEDSLTAC:
                biasFlag = (patientStatus.pthStatus not in [PthStatus.CLEAR, PthStatus.RECOVERED])
                key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA,
                       biasFlag)
                if key in self.treeCache:
                    return self.treeCache[key]
                else:
                    if biasFlag:
                        changeTree = buildChangeTree(self.pthRates,
                                                     forceRelocateDiag=DiagClassA.NEEDSLTAC)
                    else:
                        changeTree = buildChangeTree(self.lclRates,
                                                     forceRelocateDiag=DiagClassA.NEEDSLTAC)

                    tree = BayesTree(changeTree,
                                     PatientStatusSetter(),
                                     self.cachedCDF.intervalProb, tag='LOS')
                    self.treeCache[key] = tree
                    return tree
            else:
                # This patient doesn't belong in this ward
                logger.warning('fac %s patient: %s careTier %s with status %s startTime: %s: '
                               'this patient should be gone by now'
                               % (self.name, patientAgent.name, CareTier.names[ward.tier],
                                  DiagClassA.names[patientStatus.diagClassA], startTime))
                return BayesTree(PatientStatusSetter())

        else:
            raise RuntimeError('LTACs do not provide care tier %s' % ward.tier)

    def getInitialOverallHealth(self, ward, timeNow):  # @UnusedVariable
        # LTAC patients have by definition been in a hospital in the last year.
        return PatientOverallHealth.UNHEALTHY


def _populate(fac, descr, patch):
    assert 'meanPop' in descr, \
        "LTAC description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = float(descr['meanPop']['value'])
    if 'nBeds' in descr and meanPop > descr['nBeds']['value']:
        logger.warning('LTAC %s meanPop %s > nBeds %s'
                       % (descr['abbrev'], meanPop, descr['nBeds']['value']))
        meanPop = descr['nBeds']['value']
    agentList = []
    for i in xrange(int(round(meanPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.LTAC)
        assert ward is not None, 'Ran out of LTAC beds populating %(abbrev)s!' % descr
        a = PatientAgent('PatientAgent_LTAC_%s_%d' % (ward._name, i), patch, ward)
        ward.lock(a)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              None)
        agentList.append(a)
    return agentList


def generateFull(facilityDescr, patch, policyClasses=None, categoryNameMapper=None):
    fac = LTAC(facilityDescr, patch, policyClasses=policyClasses,
               categoryNameMapper=categoryNameMapper)
    return [fac], fac.getWards(), _populate(fac, facilityDescr, patch)


def estimateWork(descr):
    return hospitalEstimateWork(descr)


def checkSchema(facilityDescr):
    global _validator
    if _validator is None:
        _validator = schemautils.getValidator(_schema)
    nErrors = sum([1 for e in _validator.iter_errors(facilityDescr)])  # @UnusedVariable
    return nErrors


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(_constants_values,
                                         _constants_schema)
