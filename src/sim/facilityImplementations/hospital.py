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

import facilitybase
from facilitybase import DiagClassA, PatientStatus, CareTier
from facilitybase import Facility, PatientAgent

category = 'HOSPITAL'
_schema = 'facilityfacts_schema.yaml'
_constants_values = 'hospital_constants.yaml'
_constants_schema = 'hospital_constants_schema.yaml'
_validator = None
_constants = None


class Hospital(Facility):
    def __init__(self, descr, patch):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr, patch)
        bedsPerWard = _constants['bedsPerWard']['value']
        bedsPerICUWard = _constants['bedsPerWard']['value']
        self.hospDischargeViaDeathFrac = _constants['hospDischargeViaDeathFrac']['value']
        nTotalTransfersOut = sum([v['count']['value'] for v in descr['totalTransfersOut']])
        hospTotalTransferFrac = (float(nTotalTransfersOut)
                                 / float(descr['totalDischarges']['value']))
        self.fracDischargeHealthy = 1.0 - (hospTotalTransferFrac + self.hospDischargeViaDeathFrac)
        transferFrac = {v['category']: hospTotalTransferFrac * (float(v['count']['value'])
                                                                / float(nTotalTransfersOut))
                        for v in descr['totalTransfersOut']}
        self.fracTransferHosp = (transferFrac['HOSPITAL'] if 'HOSPITAL' in transferFrac else 0.0)
        self.fracTransferLTAC = (transferFrac['LTAC'] if 'LTAC' in transferFrac else 0.0)
        self.fracTransferNH = (transferFrac['NURSINGHOME'] if 'NURSINGHOME' in transferFrac
                               else 0.0)
        self.icuDischargeViaDeathFrac = _constants['icuDischargeViaDeathFrac']['value']
        self.hospProbCache = {}
        self.icuProbCache = {}
        if 'nBeds' in descr:
            icuBeds = descr['fracAdultPatientDaysICU'] * descr['nBeds']
            icuWards = int(icuBeds/bedsPerICUWard) + 1
            nonICUBeds = max(descr['nBeds'] - icuBeds, 0)
            nonICUWards = int(float(nonICUBeds)/bedsPerWard) + 1
#             baseStr = 'nBeds = %s' % descr['nBeds']
        else:
            meanPop = descr['meanPop']
            meanICUPop = meanPop * descr['fracAdultPatientDaysICU']
            meanNonICUPop = meanPop - meanICUPop
            icuWards = int(math.ceil(meanICUPop / bedsPerICUWard))
            nonICUWards = int(math.ceil(meanNonICUPop / bedsPerWard))
#             baseStr = 'meanPop = %s' % descr['meanPop']
        icuBeds = icuWards * bedsPerICUWard
        nonICUBeds = nonICUWards * bedsPerWard
#         dumpFName = 'dump_%s.csv' % patch.name
#         with open(dumpFName,'a') as f:
#             f.write('%s, %f, %f, %f, %f, %f\n' %
#                     (descr['abbrev'], icuBeds, nonICUBeds, meanICUPop, meanNonICUPop,
#                      (icuBeds + nonICUBeds)/(meanICUPop+meanNonICUPop)))
#         print '%s: %d + %d = %d vs. %s' % (descr['abbrev'], icuBeds, nonICUBeds,
#                                            icuBeds+nonICUBeds, baseStr)

        assert descr['losModel']['pdf'] == 'lognorm(mu=$0,sigma=$1)', \
            'Hospital %(abbrev)s LOS PDF is not lognorm(mu=$0,sigma=$1)' % descr
        self.frozenHospPDF = self._calcHospPDF(descr['losModel']['parms'][0],
                                               descr['losModel']['parms'][1])

        self.frozenICUPDF = self._calcICUPDF(descr['meanLOSICU']['value'],
                                             _constants['icuLOSLogNormSigma']['value'])

        for i in xrange(icuWards):
            self.addWard(facilitybase.Ward(('%s_%s_%s_%s_%d' %
                                            (category, patch.name, descr['abbrev'], 'ICU', i)),
                                           patch, CareTier.ICU, bedsPerICUWard))
        for i in xrange(nonICUWards):
            self.addWard(facilitybase.Ward(('%s_%s_%s_%s_%d' %
                                            (category, patch.name, descr['abbrev'], 'HOSP', i)),
                                           patch, CareTier.HOSP, bedsPerWard))

    def _calcHospPDF(self, mu, sigma):
        return lognorm(sigma, scale=math.exp(mu))

    def _calcICUPDF(self, mean, sigma):
        if mean > 0.0:
            mu = math.log(mean) - (0.5 * sigma * sigma)
            return lognorm(sigma, scale=math.exp(mu))
        else:
            return None

    def _hospChangeProb(self, startTime, timeNow):
        key = (startTime, timeNow - startTime)
        if key in self.hospProbCache:
            return self.hospProbCache[key]
        else:
            changeProb = self.frozenHospPDF.cdf(timeNow) - self.frozenHospPDF.cdf(startTime)
            self.hospProbCache[key] = changeProb
            return changeProb

    def _icuChangeProb(self, startTime, timeNow):
        key = (startTime, timeNow - startTime)
        if key in self.icuProbCache:
            return self.icuProbCache[key]
        else:
            if not self.frozenICUPDF:
                raise RuntimeError('ICU LOS distribution for %s is needed but unknown!'
                                   % self.name)
            changeProb = self.frozenICUPDF.cdf(timeNow) - self.frozenICUPDF.cdf(startTime)
            self.icuProbCache[key] = changeProb
            return changeProb

    def getStatusChangeCDF(self, patientStatus, careTier, treatment, startTime, timeNow):
        """
        Return a list of pairs of the form [(prob1, newStatus), (prob2, newStatus), ...]
        where newStatus is a full patientStatus and prob is the probability of reaching
        that status at the end of timeNow.  The sum of all the probs must equal 1.0.
        """
        assert treatment == facilitybase.TreatmentProtocol.NORMAL, \
            "Hospitals only offer 'NORMAL' treatment; found %s" % treatment
        if careTier == CareTier.HOSP:
            changeProb = self._hospChangeProb(startTime - patientStatus.startDateA,
                                              timeNow - patientStatus.startDateA)
            return [(changeProb * self.hospDischargeViaDeathFrac,
                     patientStatus._replace(diagClassA=DiagClassA.DEATH, startDateA=timeNow)),
                    (changeProb * (self.fracTransferHosp + self.fracTransferLTAC),
                     patientStatus._replace(diagClassA=DiagClassA.SICK, startDateA=timeNow)),
                    (changeProb * self.fracTransferNH,
                     patientStatus._replace(diagClassA=DiagClassA.REHAB, startDateA=timeNow)),
                    (changeProb * self.fracDischargeHealthy,
                     patientStatus._replace(diagClassA=DiagClassA.HEALTHY, startDateA=timeNow)),
                    (1.0 - changeProb, patientStatus)]
        elif careTier == CareTier.ICU:
            changeProb = self._icuChangeProb(startTime - patientStatus.startDateA,
                                             timeNow - patientStatus.startDateA)
            return [(changeProb * self.icuDischargeViaDeathFrac,
                     patientStatus._replace(diagClassA=DiagClassA.DEATH, startDateA=timeNow)),
                    (changeProb * (1.0 - self.icuDischargeViaDeathFrac),
                     patientStatus._replace(diagClassA=DiagClassA.SICK, startDateA=timeNow)),
                    (1.0 - changeProb, patientStatus)]
        else:
            raise RuntimeError('Hospitals do not provide care tier %s' % careTier)


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
        print '####### Cannot estimate work for %(abbrev)s' % facRec
        return 0


def checkSchema(facilityDescr):
    global _validator
    if _validator is None:
        with open(os.path.join(os.path.dirname(__file__), _schema), 'rU') as f:
            schemaJSON = yaml.safe_load(f)
        _validator = jsonschema.validators.validator_for(schemaJSON)(schema=schemaJSON)
    nErrors = sum([1 for e in _validator.iter_errors(facilityDescr)])  # @UnusedVariable
    return nErrors


def _importConstants():
    global _constants
    if _constants is None:
        with open(os.path.join(os.path.dirname(__file__), _constants_values), 'rU') as f:
            cJSON = yaml.safe_load(f)
        with open(os.path.join(os.path.dirname(__file__), _constants_schema), 'rU') as f:
            schemaJSON = yaml.safe_load(f)
        validator = jsonschema.validators.validator_for(schemaJSON)(schema=schemaJSON)
        nErrors = sum([1 for e in validator.iter_errors(cJSON)])  # @UnusedVariable
        if nErrors:
            for e in validator.iter_errors(cJSON):
                print ('Schema violation: %s: %s' %
                       (' '.join([str(word) for word in e.path]), e.message))
            raise RuntimeError('%s does not satisfy the schema %s' %
                               (_constants_values, _constants_schema))
        _constants = cJSON

###########
# Initialize the module
###########
_importConstants()
