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

from itertools import cycle
import os.path
import yaml
import jsonschema
from random import random
import math
from scipy.stats import lognorm, expon
import logging

import pyrheabase
import pyrheautils
from facilitybase import DiagClassA, CareTier, TreatmentProtocol, NURSINGQueue
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent, CachedCDFGenerator
from hospital import checkSchema as hospitalCheckSchema, estimateWork as hospitalEstimateWork
from hospital import createClassASetter, createOverallHealthSetter, createCopier

logger = logging.getLogger(__name__)

category = 'NURSINGHOME'
_constants_values = 'nursinghome_constants.yaml'
_constants_schema = 'nursinghome_constants_schema.yaml'
_constants = None


class NursingHome(Facility):
    def __init__(self, descr, patch):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr, patch,
                          reqQueueClasses=[NURSINGQueue])
        assert 'nBeds' in descr, 'Nursing home %(abbrev) description is missing nBeds' % descr
        nBeds = int(descr['nBeds'])
        if 'losModel' in descr:
            losModel = descr['losModel']
        else:
            losModel = _constants['nhLOSModel']
        assert losModel['pdf'] == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)', \
            "Unexpected losModel form %s for %s!" % (losModel['pdf'], descr['abbrev'])

        self.rehabCachedCDF = CachedCDFGenerator(lognorm(losModel['parms'][2],
                                                         scale=math.exp(losModel['parms'][1])))
        self.rehabTreeCache = {}
        self.residentCachedCDF = CachedCDFGenerator(expon(scale=1.0/losModel['parms'][3]))
        self.residentTreeCache = {}

        self.addWard(Ward('%s_%s_%s' % (category, patch.name, descr['abbrev']),
                          patch, CareTier.NURSING, nBeds))

    def getStatusChangeTree(self, patientStatus, careTier, treatment, startTime, timeNow):
        assert careTier == CareTier.NURSING, \
            "Nursing homes only offer CareTier 'NURSING'; found %s" % careTier
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        _c = _constants
        if treatment == TreatmentProtocol.REHAB:
            if key in self.rehabTreeCache:
                return self.rehabTreeCache[key]
            else:
                changeProb = self.rehabCachedCDF.intervalProb(*key)
                adverseProb = (_c['rehabDeathRate']['value']
                               + _c['rehabSickRate']['value']
                               + _c['rehabVerySickRate']['value'])
                tree = [changeProb,
                        [adverseProb,
                         Facility.foldCDF([(_c['rehabDeathRate']['value'] / adverseProb,
                                            createClassASetter(DiagClassA.DEATH)),
                                          (_c['rehabSickRate']['value'] / adverseProb,
                                           createClassASetter(DiagClassA.SICK)),
                                          (_c['rehabVerySickRate']['value'] / adverseProb,
                                           createClassASetter(DiagClassA.VERYSICK)),
                                           ]),
                         createClassASetter(DiagClassA.HEALTHY)],
                        createCopier()]
                self.rehabTreeCache[key] = tree
                return tree
        elif treatment == TreatmentProtocol.NORMAL:
            if key in self.residentTreeCache:
                return self.residentTreeCache[key]
            else:
                changeProb = self.residentCachedCDF.intervalProb(*key)
                totRate = (_c['residentDeathRate']['value']
                           + _c['residentSickRate']['value']
                           + _c['residentVerySickRate']['value']
                           + _c['residentReturnToCommunityRate']['value'])
                tree = [changeProb,
                        Facility.foldCDF(
                            [(_c['residentDeathRate']['value']/totRate,
                              createClassASetter(DiagClassA.DEATH)),
                             (_c['residentSickRate']['value']/totRate,
                              createClassASetter(DiagClassA.SICK)),
                             (_c['residentVerySickRate']['value']/totRate,
                              createClassASetter(DiagClassA.VERYSICK)),
                             (_c['residentReturnToCommunityRate']['value']/totRate,
                              createOverallHealthSetter(PatientOverallHealth.HEALTHY)),
                             ]),
                        createCopier()]
                self.residentTreeCache[key] = tree
                return tree
        else:
            raise RuntimeError('Nursing homes do not provide treatment protocol %s' % treatment)

    def prescribe(self, patientDiagnosis, patientTreatment):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
                return (CareTier.HOME, TreatmentProtocol.NORMAL)
            elif patientDiagnosis.overall == PatientOverallHealth.FRAIL:
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
        "Nursing home description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = descr['meanPop']
    if meanPop > descr['nBeds']:
        logger.warning('Nursing Home %(abbrev)s meanPop %(meanPop)s > nBeds %(nBeds)s'
                       % descr)
        meanPop = descr['nBeds']
    if 'losModel' in descr:
        losModel = descr['losModel']
    else:
        losModel = _constants['nhLOSModel']
    assert losModel['pdf'] == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)', \
        "Unexpected losModel form %s for %s!" % (losModel['pdf'], descr['abbrev'])
    # The following is approximate, but adequate...
    residentFrac = (1.0 - losModel['parms'][0])
    agentList = []
    for i in xrange(int(round(meanPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.NURSING)
        assert ward is not None, 'Ran out of beds populating %(abbrev)s!' % descr
        a = PatientAgent('PatientAgent_NURSING_%s_%d' % (ward._name, i), patch, ward)
        if random() <= residentFrac:
            a._status = a._status._replace(overall=PatientOverallHealth.FRAIL)
        else:  # They must be here for rehab
            a._status = a._status._replace(diagClassA=DiagClassA.NEEDSREHAB)
            a._treatment = TreatmentProtocol.REHAB
        ward.lock(a)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              0)
        agentList.append(a)
    return agentList


def generateFull(facilityDescr, patch):
    fac = NursingHome(facilityDescr, patch)
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
