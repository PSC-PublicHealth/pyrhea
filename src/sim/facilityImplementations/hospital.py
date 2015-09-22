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
from itertools import cycle
import jsonschema
import yaml
import math
from scipy.stats import lognorm
import logging

import pyrheautils
from facilitybase import DiagClassA, PatientStatus, CareTier, TreatmentProtocol
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent, CachedCDFGenerator

category = 'HOSPITAL'
_schema = 'facilityfacts_schema.yaml'
_constants_values = 'hospital_constants.yaml'
_constants_schema = 'hospital_constants_schema.yaml'
_validator = None
_constants = None

logger = logging.getLogger(__name__)


class Hospital(Facility):
    def __init__(self, descr, patch):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr, patch)
        bedsPerWard = _constants['bedsPerWard']['value']
        bedsPerICUWard = _constants['bedsPerWard']['value']
        self.hospDischargeViaDeathFrac = _constants['hospDischargeViaDeathFrac']['value']
        nTotalTransfersOut = sum([v['count']['value'] for v in descr['totalTransfersOut']])
        self.hospTotalTransferFrac = (float(nTotalTransfersOut)
                                      / float(descr['totalDischarges']['value']))
        self.fracDischargeHealthy = 1.0 - (self.hospTotalTransferFrac
                                           + self.hospDischargeViaDeathFrac)
        transferFrac = {v['category']: (self.hospTotalTransferFrac * float(v['count']['value'])
                                        / float(nTotalTransfersOut))
                        for v in descr['totalTransfersOut']}
        self.fracTransferHosp = (transferFrac['HOSPITAL'] if 'HOSPITAL' in transferFrac else 0.0)
        self.fracTransferLTAC = (transferFrac['LTAC'] if 'LTAC' in transferFrac else 0.0)
        self.fracTransferNH = (transferFrac['NURSINGHOME'] if 'NURSINGHOME' in transferFrac
                               else 0.0)
        self.icuDischargeViaDeathFrac = _constants['icuDischargeViaDeathFrac']['value']
        if 'nBeds' in descr:
            icuBeds = descr['fracAdultPatientDaysICU'] * descr['nBeds']
            icuWards = int(icuBeds/bedsPerICUWard) + 1
            nonICUBeds = max(descr['nBeds'] - icuBeds, 0)
            nonICUWards = int(float(nonICUBeds)/bedsPerWard) + 1
        else:
            meanPop = descr['meanPop']
            meanICUPop = meanPop * descr['fracAdultPatientDaysICU']
            meanNonICUPop = meanPop - meanICUPop
            icuWards = int(math.ceil(meanICUPop / bedsPerICUWard))
            nonICUWards = int(math.ceil(meanNonICUPop / bedsPerWard))

        icuBeds = icuWards * bedsPerICUWard
        nonICUBeds = nonICUWards * bedsPerWard

        assert descr['losModel']['pdf'] == 'lognorm(mu=$0,sigma=$1)', \
            'Hospital %(abbrev)s LOS PDF is not lognorm(mu=$0,sigma=$1)' % descr
        self.hospCachedCDF = CachedCDFGenerator(lognorm(descr['losModel']['parms'][0],
                                                        scale=math.exp(descr['losModel']
                                                                       ['parms'][1])))
        self.hospTreeCache = {}
        if descr['meanLOSICU']['value'] == 0.0:
            self.icuCachedCDF = None
        else:
            mean = descr['meanLOSICU']['value']
            sigma = _constants['icuLOSLogNormSigma']['value']
            mu = math.log(mean) - (0.5 * sigma * sigma)
            self.icuCachedCDF = CachedCDFGenerator(lognorm(sigma, scale=math.exp(mu)))
        self.icuTreeCache = {}

        for i in xrange(icuWards):
            self.addWard(Ward(('%s_%s_%s_%s_%d' %
                               (category, patch.name, descr['abbrev'], 'ICU', i)),
                              patch, CareTier.ICU, bedsPerICUWard))
        for i in xrange(nonICUWards):
            self.addWard(Ward(('%s_%s_%s_%s_%d' %
                               (category, patch.name, descr['abbrev'], 'HOSP', i)),
                              patch, CareTier.HOSP, bedsPerWard))

    def getStatusChangeTree(self, patientStatus, careTier, treatment, startTime, timeNow):
        assert treatment == TreatmentProtocol.NORMAL, \
            "Hospitals only offer 'NORMAL' treatment; found %s" % treatment
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        if careTier == CareTier.HOSP:
            if key in self.hospTreeCache:
                return self.hospTreeCache[key]
            else:
                changeProb = self.hospCachedCDF.intervalProb(*key)
                tree = [changeProb,
                        Hospital.foldCDF([(self.hospDischargeViaDeathFrac,
                                           patientStatus._replace(diagClassA=DiagClassA.DEATH,
                                                                  startDateA=timeNow)),
                                          ((self.fracTransferHosp + self.fracTransferLTAC),
                                           patientStatus._replace(diagClassA=DiagClassA.SICK,
                                                                  startDateA=timeNow)),
                                          (((self.fracTransferNH + self.fracDischargeHealthy)
                                            * _constants['fracOfDischargesRequiringRehab']['value']),
                                           patientStatus._replace(diagClassA=DiagClassA.NEEDSREHAB,
                                                                  startDateA=timeNow)),
                                          (((self.fracTransferNH + self.fracDischargeHealthy)
                                            * (1.0 - _constants['fracOfDischargesRequiringRehab']['value'])),
                                           patientStatus._replace(diagClassA=DiagClassA.HEALTHY,
                                                                  startDateA=timeNow))]),
                        patientStatus]
                self.hospTreeCache[key] = tree
                return tree
        elif careTier == CareTier.ICU:
            if key in self.icuTreeCache:
                return self.icuTreeCache[key]
            else:
                changeProb = self.icuCachedCDF.intervalProb(*key)
#                 changeProb = self._icuChangeProb(startTime - patientStatus.startDateA,
#                                                  timeNow - patientStatus.startDateA)
                tree = [changeProb,
                        Hospital.foldCDF([(self.icuDischargeViaDeathFrac,
                                           patientStatus._replace(diagClassA=DiagClassA.DEATH,
                                                                  startDateA=timeNow)),
                                          ((1.0 - self.icuDischargeViaDeathFrac),
                                           patientStatus._replace(diagClassA=DiagClassA.SICK,
                                                                  startDateA=timeNow))]),
                        patientStatus]
                self.icuTreeCache[key] = tree
                return tree
        else:
            raise RuntimeError('Hospitals do not provide care tier %s' % careTier)

    def prescribe(self, patientDiagnosis, patientTreatment):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
                return (CareTier.HOME, TreatmentProtocol.NORMAL)
            else:
                return (CareTier.NURSING, TreatmentProtocol.NORMAL)
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
            return (CareTier.NURSING, TreatmentProtocol.REHAB)
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
        "Hospital description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = float(descr['meanPop'])
    meanICUPop = meanPop * descr['fracAdultPatientDaysICU']
    meanHospPop = meanPop - meanICUPop
    icuWards = fac.getWards(CareTier.ICU)
    hospWards = fac.getWards(CareTier.HOSP)
    agentList = []
    for i, ward in zip(xrange(int(round(meanICUPop))), cycle(icuWards)):
        a = PatientAgent('PatientAgent_ICU_%s_%d' % (ward._name, i),
                         patch, ward)
        ward.lock(a)
        agentList.append(a)
    for i, ward in zip(xrange(int(round(meanHospPop))), cycle(hospWards)):
        a = PatientAgent('PatientAgent_HOSP_%s_%d' % (ward._name, i),
                         patch, ward)
        ward.lock(a)
        agentList.append(a)
    return agentList


def generateFull(facilityDescr, patch):
    fac = Hospital(facilityDescr, patch)
    return [fac], fac.getWards(), _populate(fac, facilityDescr, patch)


def estimateWork(facRec):
    if 'meanPop' in facRec:
        return facRec['meanPop']
    elif 'nBeds' in facRec:
        return facRec['nBeds']
    else:
        logger.warn('Cannot estimate work for %(abbrev)s' % facRec)
        return 0


def checkSchema(facilityDescr):
    global _validator
    if _validator is None:
        with open(os.path.join(os.path.dirname(__file__), _schema), 'rU') as f:
            schemaJSON = yaml.safe_load(f)
        _validator = jsonschema.validators.validator_for(schemaJSON)(schema=schemaJSON)
    nErrors = sum([1 for e in _validator.iter_errors(facilityDescr)])  # @UnusedVariable
    return nErrors


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(os.path.join(os.path.dirname(__file__),
                                                      _constants_values),
                                         os.path.join(os.path.dirname(__file__),
                                                      _constants_schema))
