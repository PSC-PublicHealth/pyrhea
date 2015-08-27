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

category = 'HOSPITAL'
schema = 'facilityfacts_schema.yaml'
constants_values = 'hospital_constants.yaml'
constants_schema = 'hospital_constants_schema.yaml'

validator = None
constants = None

# """An empirical estimate"""
# bedsPerWard = 20
# """Another empirical estimate"""
# bedsPerICUWard = 12

# class CommunityWard(Ward):
#     """This 'ward' type represents being out in the community"""
#
#     def __init__(self, name, patch, nBeds):
#         Ward.__init__(self, name, patch, CareTier.HOME, nBeds)
#         self.checkInterval = 30  # check health monthly
#         #self.checkInterval = 1  # check health monthly


class Hospital(facilitybase.Facility):
    def __init__(self, descr, patch):
        facilitybase.Facility.__init__(self, '%(category)s_%(abbrev)s' % descr, patch)
        bedsPerWard = constants['bedsPerWard']['value']
        bedsPerICUWard = constants['bedsPerWard']['value']
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
        mu = descr['losModel']['parms'][0]
        sigma = descr['losModel']['parms'][1]
        self.frozenPDF = lognorm(sigma, scale=math.exp(mu))
        for i in xrange(icuWards):
            self.addWard(facilitybase.Ward(('%s_%s_%s_%s_%d' %
                                            (category, patch.name, descr['abbrev'], 'ICU', i)),
                                           patch, facilitybase.CareTier.ICU, bedsPerICUWard))
        for i in xrange(nonICUWards):
            self.addWard(facilitybase.Ward(('%s_%s_%s_%s_%d' %
                                            (category, patch.name, descr['abbrev'], 'HOSP', i)),
                                           patch, facilitybase.CareTier.HOSP, bedsPerWard))

    def getStatusChangeCDF(self, patientStatus, careTier, treatment, startTime, timeNow):
        """
        Return a list of pairs of the form [(prob1, newStatus), (prob2, newStatus), ...]
        where newStatus is a full patientStatus and prob is the probability of reaching
        that status at the end of timeNow.  The sum of all the probs must equal 1.0.
        """
        if careTier == facilitybase.CareTier.HOSP:
            changeProb = self.frozenPDF.cdf(timeNow) - self.frozenPDF.cdf(startTime)
            fracDischargesViaDeath = 0.0081
            return [(changeProb * fracDischargesViaDeath,
                     facilitybase.PatientStatus(facilitybase.DiagClassA.DEATH, timeNow,
                                                patientStatus.diagClassB,
                                                patientStatus.startDateB)),
                    (changeProb * (1.0-fracDischargesViaDeath),
                     facilitybase.PatientStatus(facilitybase.DiagClassA.HEALTHY, timeNow,
                                                patientStatus.diagClassB,
                                                patientStatus.startDateB)),
                    (1.0 - changeProb, patientStatus)]
        elif careTier == facilitybase.CareTier.ICU:
            return [(1.0, patientStatus)]
        else:
            raise RuntimeError('Hospitals do not provide care tier %s' % careTier)


def _populate(fac, descr, patch):
    assert 'meanPop' in descr, \
        "Hospital description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = float(descr['meanPop'])
    meanICUPop = meanPop * descr['fracAdultPatientDaysICU']
    meanHospPop = meanPop - meanICUPop
    icuWards = fac.getWards(facilitybase.CareTier.ICU)
    hospWards = fac.getWards(facilitybase.CareTier.HOSP)
    agentList = []
    for i, ward in zip(xrange(int(round(meanICUPop))), cycle(icuWards)):
        a = facilitybase.PatientAgent('PatientAgent_ICU_%s_%d' % (ward._name, i),
                                      patch, ward)
        ward.lock(a)
        agentList.append(a)
    for i, ward in zip(xrange(int(round(meanHospPop))), cycle(hospWards)):
        a = facilitybase.PatientAgent('PatientAgent_HOSP_%s_%d' % (ward._name, i),
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
    global validator
    if validator is None:
        with open(os.path.join(os.path.dirname(__file__), schema), 'rU') as f:
            schemaJSON = yaml.safe_load(f)
        validator = jsonschema.validators.validator_for(schemaJSON)(schema=schemaJSON)
    nErrors = sum([1 for e in validator.iter_errors(facilityDescr)])  # @UnusedVariable
    return nErrors


def _importConstants():
    global constants
    if constants is None:
        with open(os.path.join(os.path.dirname(__file__), constants_values), 'rU') as f:
            cJSON = yaml.safe_load(f)
        with open(os.path.join(os.path.dirname(__file__), constants_schema), 'rU') as f:
            schemaJSON = yaml.safe_load(f)
        validator = jsonschema.validators.validator_for(schemaJSON)(schema=schemaJSON)
        nErrors = sum([1 for e in validator.iter_errors(cJSON)])  # @UnusedVariable
        if nErrors:
            for e in validator.iter_errors(cJSON):
                print ('Schema violation: %s: %s' %
                       (' '.join([str(word) for word in e.path]), e.message))
            raise RuntimeError('%s does not satisfy the schema %s' %
                               (constants_values, constants_schema))
        constants = cJSON

###########
# Initialize the module
###########
_importConstants()
