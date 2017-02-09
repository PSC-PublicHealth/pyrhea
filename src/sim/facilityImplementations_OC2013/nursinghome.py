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
from random import random
import math
from scipy.stats import lognorm, expon
import logging

import pyrheabase
import pyrheautils
import schemautils
from stats import CachedCDFGenerator, lognormplusexp, BayesTree
from facilitybase import DiagClassA, CareTier, TreatmentProtocol, NURSINGQueue
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent
from facilitybase import PatientStatusSetter
from hospital import estimateWork as hospitalEstimateWork
from hospital import ClassASetter, OverallHealthSetter

logger = logging.getLogger(__name__)

category = 'NURSINGHOME'
_schema = 'nursinghomefacts_schema.yaml'
_constants_values = '$(MODELDIR)/constants/nursinghome_constants.yaml'
_constants_schema = 'nursinghome_constants_schema.yaml'
_constants = None
_validator = None


class NursingHome(Facility):
    def __init__(self, descr, patch, policyClasses=None, categoryNameMapper=None):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr,
                          descr, patch,
                          reqQueueClasses=[NURSINGQueue],
                          policyClasses=policyClasses,
                          categoryNameMapper=categoryNameMapper)
        descr = self.mapDescrFields(descr)
        _c = _constants
        if 'nBeds' in descr:
            nBeds = int(descr['nBeds']['value'])
        else:
            assert 'meanPop' in descr, (('Nursing home %(abbrev)s description has neither'
                                         ' nBeds nor meanPop') % descr)
            nBeds = int(math.ceil(descr['meanPop']['value']))
            logger.warning('Nursing Home %s has no nBeds; using ceil(meanPop) = %d'
                           % (descr['abbrev'], nBeds))
        if 'losModel' in descr:
            losModel = descr['losModel']
        else:
            losModel = _constants['nhLOSModel']
        assert losModel['pdf'] == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)', \
            "Unexpected losModel form %s for %s!" % (losModel['pdf'], descr['abbrev'])

        self.rehabCachedCDF = CachedCDFGenerator(lognorm(losModel['parms'][2],
                                                         scale=math.exp(losModel['parms'][1])))
        self.rehabTreeCache = {}
        lMP = losModel['parms']
#         self.frailCachedCDF = CachedCDFGenerator(lognormplusexp(scale=math.exp(lMP[1]),
#                                                                 s=lMP[2],
#                                                                 k=lMP[0],
#                                                                 lmda=lMP[3]))
        self.frailCachedCDF = CachedCDFGenerator(expon(scale=1.0/lMP[3]))
        self.frailTreeCache = {}
        self.frailRehabTreeCache = {}
        self.addWard(Ward('%s_%s_%s' % (category, patch.name, descr['abbrev']),
                          patch, CareTier.NURSING, nBeds))

    def getStatusChangeTree(self, patientStatus, ward, treatment, startTime, timeNow):
        careTier = ward.tier
        assert careTier == CareTier.NURSING, \
            "Nursing homes only offer CareTier 'NURSING'; found %s" % careTier
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        _c = _constants
        if patientStatus.overall == PatientOverallHealth.FRAIL:
            if treatment.rehab:
                if key in self.frailRehabTreeCache:
                    return self.frailRehabTreeCache[key]
                else:
                    # If no major changes happen, allow for the patient's completing rehab
                    innerTree = BayesTree(ClassASetter(PatientOverallHealth.HEALTHY),
                                          PatientStatusSetter(),
                                          self.rehabCachedCDF.intervalProb(*key))
            else:
                if key in self.frailTreeCache:
                    return self.frailTreeCache[key]
                else:
                    innerTree = PatientStatusSetter()  # No change of state
            changeProb = self.frailCachedCDF.intervalProb(*key)
            totRate = (_c['residentDeathRate']['value']
                       + _c['residentSickRate']['value']
                       + _c['residentVerySickRate']['value']
                       + _c['residentNeedsLTACRate']['value']
                       + _c['residentReturnToCommunityRate']['value'])
            healthySetter = OverallHealthSetter(PatientOverallHealth.HEALTHY)
            changeTree = BayesTree.fromLinearCDF([(_c['residentDeathRate']['value']/totRate,
                                                   ClassASetter(DiagClassA.DEATH)),
                                                  (_c['residentNeedsLTACRate']['value']/totRate,
                                                   ClassASetter(DiagClassA.NEEDSLTAC)),
                                                  (_c['residentSickRate']['value']/totRate,
                                                   ClassASetter(DiagClassA.SICK)),
                                                  (_c['residentVerySickRate']['value']/totRate,
                                                   ClassASetter(DiagClassA.VERYSICK)),
                                                  (_c['residentReturnToCommunityRate']['value']
                                                   / totRate,
                                                   healthySetter)])
            tree = BayesTree(changeTree,
                             innerTree,
                             changeProb)
            if treatment.rehab:
                self.frailRehabTreeCache[key] = tree
            else:
                self.frailTreeCache[key] = tree
            return tree
        else:
            if treatment.rehab:
                if key in self.rehabTreeCache:
                    return self.rehabTreeCache[key]
                else:
                    changeProb = self.rehabCachedCDF.intervalProb(*key)
                    adverseProb = (_c['rehabDeathRate']['value']
                                   + _c['rehabSickRate']['value']
                                   + _c['rehabVerySickRate']['value']
                                   + _c['rehabNeedsLTACRate']['value'])
                    adverseTree = BayesTree.fromLinearCDF([(_c['rehabDeathRate']['value']
                                                            / adverseProb,
                                                            ClassASetter(DiagClassA.DEATH)),
                                                           (_c['rehabNeedsLTACRate']['value']
                                                            / adverseProb,
                                                            ClassASetter(DiagClassA.NEEDSLTAC)),
                                                           (_c['rehabSickRate']['value']
                                                            / adverseProb,
                                                            ClassASetter(DiagClassA.SICK)),
                                                           (_c['rehabVerySickRate']['value']
                                                            / adverseProb,
                                                            ClassASetter(DiagClassA.VERYSICK))])
                    tree = BayesTree(BayesTree(adverseTree,
                                               ClassASetter(DiagClassA.HEALTHY),
                                               adverseProb),
                                     PatientStatusSetter(),
                                     changeProb)
                    self.rehabTreeCache[key] = tree
                    return tree
            elif not treatment.rehab:
                if patientStatus.diagClassA in ([DiagClassA.NEEDSLTAC, DiagClassA.SICK,
                                                 DiagClassA.VERYSICK]):
                    logger.warning('fac %s status: %s careTier %s startTime: %s: '
                                   'this patient should be gone by now'
                                   % (self.name, str(patientStatus), CareTier.names[careTier],
                                      startTime))
                    return BayesTree(PatientStatusSetter())
                else:
                    raise RuntimeError('Patients with NORMAL overall health should only be'
                                       ' in NURSING care for rehab')
            else:
                raise RuntimeError('Nursing homes do not provide treatment %s'
                                   % treatment)

    def prescribe(self, patientDiagnosis, patientTreatment):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.HEALTHY:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
                return (CareTier.HOME, patientTreatment._replace(rehab=False))
            elif patientDiagnosis.overall == PatientOverallHealth.FRAIL:
                return (CareTier.NURSING, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSREHAB:
            return (CareTier.NURSING, patientTreatment._replace(rehab=True))
        elif patientDiagnosis.diagClassA == DiagClassA.SICK:
            return (CareTier.HOSP, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.NEEDSLTAC:
            return (CareTier.LTAC, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.VERYSICK:
            return (CareTier.ICU, patientTreatment._replace(rehab=False))
        elif patientDiagnosis.diagClassA == DiagClassA.DEATH:
            return (None, patientTreatment._replace(rehab=False))
        else:
            raise RuntimeError('Unknown DiagClassA %s' % str(patientDiagnosis.diagClassA))


def _populate(fac, descr, patch):
    assert 'meanPop' in descr, \
        "Nursing home description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = descr['meanPop']['value']
    if 'nBeds' in descr and meanPop > descr['nBeds']['value']:
        logger.warning('Nursing Home %s meanPop %s > nBeds %s'
                       % (descr['abbrev'], meanPop, descr['nBeds']['value']))
        meanPop = descr['nBeds']['value']
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
            a._treatment = a._treatment._replace(rehab=True)
        ward.lock(a)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              0)
        agentList.append(a)
    return agentList


def generateFull(facilityDescr, patch, policyClasses=None, categoryNameMapper=None):
    fac = NursingHome(facilityDescr, patch, policyClasses=policyClasses,
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
