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
from facilitybase import DiagClassA, CareTier, TreatmentProtocol
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent
from facilitybase import PatientStatusSetter, LTACQueue, tierToQueueMap
from stats import CachedCDFGenerator, BayesTree
from hospital import ClassASetter
from hospital import checkSchema as hospitalCheckSchema, estimateWork as hospitalEstimateWork

category = 'LTAC'
_schema = 'facilityfacts_schema.yaml'
_constants_values = 'ltac_constants.yaml'
_constants_schema = 'ltac_constants_schema.yaml'
_validator = None
_constants = None

logger = logging.getLogger(__name__)


class LTAC(Facility):
    def __init__(self, descr, patch):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr, patch,
                          reqQueueClasses=[LTACQueue])
        bedsPerWard = _constants['bedsPerWard']['value']
        self.dischargeViaDeathFrac = _constants['dischargeViaDeathFrac']['value']
        nTotalTransfersOut = sum([v['count']['value'] for v in descr['totalTransfersOut']])
        self.totalTransferFrac = (float(nTotalTransfersOut)
                                  / float(descr['totalDischarges']['value']))
        self.fracDischargeHealthy = 1.0 - (self.totalTransferFrac
                                           + self.dischargeViaDeathFrac)
        transferFrac = {v['category']: (self.totalTransferFrac * float(v['count']['value'])
                                        / float(nTotalTransfersOut))
                        for v in descr['totalTransfersOut']}
        self.fracTransferHosp = (transferFrac['HOSPITAL'] if 'HOSPITAL' in transferFrac else 0.0)
        self.fracTransferLTAC = (transferFrac['LTAC'] if 'LTAC' in transferFrac else 0.0)
        self.fracTransferNH = (transferFrac['NURSINGHOME'] if 'NURSINGHOME' in transferFrac
                               else 0.0)
        if 'nBeds' in descr:
            nBeds = descr['nBeds']
            nWards = int(float(nBeds)/bedsPerWard) + 1
        else:
            meanPop = descr['meanPop']
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

    def getOrderedCandidateFacList(self, oldTier, newTier, timeNow):
        queueClass = tierToQueueMap[newTier]
        facAddrList = [tpl[1] for tpl in self.manager.patch.serviceLookup(queueClass.__name__)]
        shuffle(facAddrList)
        return facAddrList

    def getStatusChangeTree(self, patientStatus, ward, treatment, startTime, timeNow):
        assert treatment == TreatmentProtocol.NORMAL, \
            "LTACs only offer 'NORMAL' treatment; found %s" % treatment
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
                                                       ClassASetter(DiagClassA.HEALTHY))])
                tree = BayesTree(changeTree,
                                 PatientStatusSetter(),
                                 changeProb)
                self.treeCache[key] = tree
                return tree
        else:
            raise RuntimeError('LTACs do not provide care tier %s' % careTier)

    def prescribe(self, patientDiagnosis, patientTreatment):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
                return (CareTier.HOME, TreatmentProtocol.NORMAL)
            else:
                return (CareTier.NURSING, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
            return (CareTier.NURSING, TreatmentProtocol.REHAB)
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSLTAC:
            return (CareTier.LTAC, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.SICK:
            return (CareTier.HOSP, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.VERYSICK:
            return (CareTier.ICU, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.DEATH:
            return (None, TreatmentProtocol.NORMAL)
        else:
            raise RuntimeError('Unknown DiagClassA %s' % str(patientDiagnosis.diagClassA))


def _populate(fac, descr, patch):
    assert 'meanPop' in descr, \
        "LTAC description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = float(descr['meanPop'])
    agentList = []
    for i in xrange(int(round(meanPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.LTAC)
        assert ward is not None, 'Ran out of LTAC beds populating %(abbrev)s!' % descr
        a = PatientAgent('PatientAgent_LTAC_%s_%d' % (ward._name, i), patch, ward)
        ward.lock(a)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              0)
        agentList.append(a)
    return agentList


def generateFull(facilityDescr, patch):
    fac = LTAC(facilityDescr, patch)
    return [fac], fac.getWards(), _populate(fac, facilityDescr, patch)


def estimateWork(descr):
    return hospitalEstimateWork(descr)


def checkSchema(facilityDescr):
    return hospitalCheckSchema(facilityDescr)


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(os.path.join(os.path.dirname(__file__),
                                                      _constants_values),
                                         os.path.join(os.path.dirname(__file__),
                                                      _constants_schema))
