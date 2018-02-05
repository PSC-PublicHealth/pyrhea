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
from facilitybase import DiagClassA, CareTier, TreatmentProtocol
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent
from facilitybase import PatientStatusSetter, LTACQueue, tierToQueueMap
from stats import CachedCDFGenerator, BayesTree
from hospital import ClassASetter
from hospital import estimateWork as hospitalEstimateWork

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
        self.dischargeViaDeathFrac = _constants['dischargeViaDeathFrac']['value']
        nTotalTransfersOut = sum([v['count']['value'] for v in descr['totalTransfersOut']])
        self.totalTransferFrac = (float(nTotalTransfersOut)
                                  / float(descr['totalDischarges']['value']))
        self.fracDischargeHealthy = 1.0 - (self.totalTransferFrac
                                           + self.dischargeViaDeathFrac)
        transferFrac = {}
        for v in descr['totalTransfersOut']:
            implCat = v['category']
            if implCat not in transferFrac:
                transferFrac[implCat] = 0.0
            transferFrac[implCat] += (self.totalTransferFrac * float(v['count']['value'])
                                      / float(nTotalTransfersOut))
        self.fracTransferHosp = (transferFrac['HOSPITAL'] if 'HOSPITAL' in transferFrac else 0.0)
        self.fracTransferLTAC = (transferFrac['LTAC'] if 'LTAC' in transferFrac else 0.0)
        self.fracTransferNH = (transferFrac['NURSINGHOME'] if 'NURSINGHOME' in transferFrac
                               else 0.0)
        if 'nBeds' in descr:
            nBeds = descr['nBeds']['value']
            nWards = int(float(nBeds)/bedsPerWard) + 1
        else:
            meanPop = descr['meanPop']['value']
            nWards = int(math.ceil(meanPop / bedsPerWard))

        nBeds = nWards * bedsPerWard

        assert descr['losModel']['pdf'] == 'lognorm(mu=$0,sigma=$1)', \
            'LTAC %(abbrev)s LOS PDF is not lognorm(mu=$0,sigma=$1)' % descr
        self.cachedCDF = CachedCDFGenerator(lognorm(descr['losModel']['parms'][1],
                                                    scale=math.exp(descr['losModel']
                                                                   ['parms'][0])))
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

    def getStatusChangeTree(self, patientStatus, ward, treatment, startTime, timeNow):
        assert not treatment.rehab, \
            "LTACs do not offer rehab treatment; found %s" % treatment
        careTier = ward.tier
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        _c = _constants
        if careTier == CareTier.LTAC:
            if key in self.treeCache:
                return self.treeCache[key]
            else:
                changeProb = self.cachedCDF.intervalProb(*key)
                rehabFrac = _c['fracOfDischargesRequiringRehab']['value']
                changeTree = BayesTree.fromLinearCDF([(self.dischargeViaDeathFrac,
                                                       ClassASetter(DiagClassA.DEATH)),
                                                      (self.fracTransferHosp,
                                                       ClassASetter(DiagClassA.SICK)),
                                                      (self.fracTransferLTAC,
                                                       ClassASetter(DiagClassA.NEEDSLTAC)),
                                                      (((self.fracTransferNH
                                                         + self.fracDischargeHealthy)
                                                        * rehabFrac),
                                                       ClassASetter(DiagClassA.NEEDSREHAB)),
                                                      (((self.fracTransferNH
                                                         + self.fracDischargeHealthy)
                                                        * (1.0 - rehabFrac)),
                                                       ClassASetter(DiagClassA.WELL))])
                tree = BayesTree(changeTree,
                                 PatientStatusSetter(),
                                 changeProb)
                self.treeCache[key] = tree
                return tree
        else:
            raise RuntimeError('LTACs do not provide care tier %s' % careTier)

    def prescribe(self, ward, patientId, patientDiagnosis, patientTreatment,
                  timeNow=None):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.WELL:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
                return (CareTier.HOME, patientTreatment._replace(rehab=False))
            else:
                return (CareTier.NURSING, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
            return (CareTier.NURSING, patientTreatment._replace(rehab=True))
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSLTAC:
            return (CareTier.LTAC, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.SICK:
            return (CareTier.HOSP, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.VERYSICK:
            return (CareTier.ICU, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.DEATH:
            return (None, patientTreatment._replace(rehab=False))
        else:
            raise RuntimeError('Unknown DiagClassA %s' % str(patientDiagnosis.diagClassA))


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
