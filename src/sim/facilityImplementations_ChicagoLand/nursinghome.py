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
_schema = 'nursinghomefacts_ChicagoLand_schema.yaml'
_constants_values = '$(MODELDIR)/constants/nursinghome_constants.yaml'
_constants_schema = 'nursinghome_ChicagoLand_constants_schema.yaml'
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
        nBeds = int(descr['nBeds']['value'])
        losModel = descr['losModel']
        assert losModel['pdf'] == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)', \
            "Unexpected losModel form %s for %s!" % (losModel['pdf'], descr['abbrev'])

        totDsch = float(descr['totalDischarges']['value'])
        totTO = sum([elt['count']['value'] for elt in descr['totalTransfersOut']])
        ttoD = {elt['category']: elt['count']['value'] for elt in descr['totalTransfersOut']}
        fracNotTransferred = (float(totDsch) - float(totTO)) / float(totDsch)
        assert fracNotTransferred >= 0.0, '%s has more transfers than discharges' % self.name

        self.lclRates = {}
        self.lclRates['death'] = _c['deathRate']['value']
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

        self.lclRates['icu'] = _c['hospTransferToICURate']['value'] * self.lclRates['hospital']
        self.lclRates['hospital'] = ((1.0 -_c['hospTransferToICURate']['value'])
                                     * self.lclRates['hospital'])

        lMP = losModel['parms']
        self.rehabCachedCDF = CachedCDFGenerator(lognorm(lMP[2],
                                                         scale=math.exp(lMP[1])))
        self.rehabTreeCache = {}
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
            if treatment == TreatmentProtocol.REHAB:
                if key in self.frailRehabTreeCache:
                    return self.frailRehabTreeCache[key]
                else:
                    # If no major changes happen, allow for the patient's completing rehab
                    innerTree = BayesTree(ClassASetter(DiagClassA.HEALTHY),
                                          PatientStatusSetter(),
                                          self.rehabCachedCDF.intervalProb(*key))
            else:
                if key in self.frailTreeCache:
                    return self.frailTreeCache[key]
                else:
                    innerTree = PatientStatusSetter()  # No change of state
            changeProb = self.frailCachedCDF.intervalProb(*key)
            healthySetter = OverallHealthSetter(PatientOverallHealth.HEALTHY)
            changeTree = BayesTree.fromLinearCDF([(self.lclRates['death'],
                                                   ClassASetter(DiagClassA.DEATH)),
                                                  (self.lclRates['nursinghome'],
                                                   ClassASetter(DiagClassA.NEEDSREHAB)),
                                                  (self.lclRates['vsnf'],
                                                   ClassASetter(DiagClassA.NEEDSSKILNRS)),
                                                  (self.lclRates['ltac'],
                                                   ClassASetter(DiagClassA.NEEDSLTAC)),
                                                  (self.lclRates['hospital'],
                                                   ClassASetter(DiagClassA.SICK)),
                                                  (self.lclRates['icu'],
                                                   ClassASetter(DiagClassA.VERYSICK)),
                                                  (self.lclRates['home'],
                                                   healthySetter)])

            tree = BayesTree(changeTree,
                             innerTree,
                             changeProb)
            if treatment == TreatmentProtocol.REHAB:
                self.frailRehabTreeCache[key] = tree
            else:
                self.frailTreeCache[key] = tree
            return tree
        else:
            if treatment == TreatmentProtocol.REHAB:
                if key in self.rehabTreeCache:
                    return self.rehabTreeCache[key]
                else:
                    changeProb = self.rehabCachedCDF.intervalProb(*key)
                    adverseProb = (self.lclRates['death']
                                   + self.lclRates['hospital']
                                   + self.lclRates['icu']
                                   + self.lclRates['ltac']
                                   + self.lclRates['nursinghome']
                                   + self.lclRates['vsnf'])
                    if adverseProb > 0.0:
                        adverseTree = BayesTree.fromLinearCDF([(self.lclRates['death']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.DEATH)),
                                                               (self.lclRates['nursinghome']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.NEEDSREHAB)),
                                                               (self.lclRates['vsnf']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.NEEDSSKILNRS)),
                                                               (self.lclRates['ltac']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.NEEDSLTAC)),
                                                               (self.lclRates['hospital']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.SICK)),
                                                               (self.lclRates['icu']
                                                                / adverseProb,
                                                                ClassASetter(DiagClassA.VERYSICK))])
                        tree = BayesTree(BayesTree(adverseTree,
                                                   ClassASetter(DiagClassA.HEALTHY),
                                                   adverseProb),
                                         PatientStatusSetter(),
                                         changeProb)
                    else:
                        tree = BayesTree(ClassASetter(DiagClassA.HEALTHY),
                                         PatientStatusSetter(),
                                         changeProb)
                    self.rehabTreeCache[key] = tree
                    return tree
            elif treatment == TreatmentProtocol.NORMAL:
                if patientStatus.diagClassA in ([DiagClassA.NEEDSLTAC, DiagClassA.SICK,
                                                 DiagClassA.VERYSICK, DiagClassA.NEEDSVENT,
                                                 DiagClassA.NEEDSSKILNRS]):
                    logger.warning('fac %s status: %s careTier %s startTime: %s: '
                                   'this patient should be gone by now'
                                   % (self.name, str(patientStatus), CareTier.names[careTier],
                                      startTime))
                    return BayesTree(PatientStatusSetter())
                else:
                    raise RuntimeError('Patients with NORMAL overall health should only be'
                                       ' in NURSING care for rehab')
            else:
                raise RuntimeError('Nursing homes do not provide treatment protocol %s'
                                   % treatment)


def _populate(fac, descr, patch):
    assert 'meanPop' in descr, \
        "Nursing home description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = descr['meanPop']['value']
    if 'nBeds' in descr and meanPop > descr['nBeds']['value']:
        logger.warning('Nursing Home %s meanPop %s > nBeds %s'
                       % (descr['abbrev'], meanPop, descr['nBeds']['value']))
        meanPop = descr['nBeds']['value']
    losModel = descr['losModel']
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
