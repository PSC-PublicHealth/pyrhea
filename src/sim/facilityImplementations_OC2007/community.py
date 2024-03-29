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
from scipy.stats import expon
import logging

import pyrheabase
import pyrheautils
from facilitybase import DiagClassA, CareTier, TreatmentProtocol, BirthQueue, HOMEQueue
from facilitybase import PatientOverallHealth, Facility, Ward, PatientAgent, PatientStatusSetter
from stats import CachedCDFGenerator, BayesTree
from hospital import ClassASetter
import schemautils

logger = logging.getLogger(__name__)

category = 'COMMUNITY'
_schema = 'communityfacts_schema.yaml'
_validator = None
_constants_values = 'community_constants.yaml'
_constants_schema = 'community_constants_schema.yaml'
_constants = None


class CommunityWard(Ward):
    """This 'ward' type represents being out in the community"""

    def __init__(self, name, patch, nBeds):
        Ward.__init__(self, name, patch, CareTier.HOME, nBeds)
        self.checkInterval = _constants['communityPatientCheckInterval']['value']


class Community(Facility):
    def __init__(self, descr, patch, policyClasses=None):
        Facility.__init__(self, '%(category)s_%(abbrev)s' % descr,
                          descr, patch,
                          reqQueueClasses=[pyrheabase.FacRequestQueue, BirthQueue, HOMEQueue],
                          policyClasses=policyClasses)
        meanPop = descr['meanPop']
        nBeds = int(round(3.0*meanPop))
        losModel = _constants['communityLOSModel']
        assert losModel['pdf'] == 'expon(lambda=$0)', \
            "Unexpected losModel form %s for %s!" % (losModel['pdf'], descr['abbrev'])
        self.addWard(CommunityWard('%s_%s_%s' % (category, patch.name, descr['abbrev']),
                                   patch, nBeds))
        self.cachedCDF = CachedCDFGenerator(expon(scale=1.0/losModel['parms'][0]))
        self.treeCache = {}

    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        self.hospTreeCache = {}

    def getStatusChangeTree(self, patientAgent, startTime, timeNow):
        patientStatus = patientAgent.getStatus()
        ward = patientAgent.ward
        treatment = patientAgent.getTreatmentProtocol()
        careTier = ward.tier
        assert careTier == CareTier.HOME, \
            "The community only offers CareTier 'HOME'; found %s" % careTier
        assert not any(treatment), \
            "The community offers no treatment; found %s" % treatment
        key = (startTime - patientStatus.startDateA, timeNow - patientStatus.startDateA)
        if key in self.treeCache:
            return self.treeCache[key]
        else:
            changeProb = self.cachedCDF.intervalProb(*key)
            deathRate = _constants['communityDeathRate']['value']
            verySickRate = _constants['communityVerySickRate']['value']
            needsLTACRate = _constants['communityNeedsLTACRate']['value']
            sickRate = 1.0 - (deathRate + verySickRate + needsLTACRate)
            tree = BayesTree(BayesTree.fromLinearCDF([(deathRate,
                                                       ClassASetter(DiagClassA.DEATH)),
                                                      (sickRate,
                                                       ClassASetter(DiagClassA.SICK)),
                                                      (verySickRate,
                                                       ClassASetter(DiagClassA.VERYSICK)),
                                                      (needsLTACRate,
                                                       ClassASetter(DiagClassA.NEEDSLTAC)),
                                                      ]),
                             PatientStatusSetter(),
                             changeProb)
            self.treeCache[key] = tree
            return tree

    def prescribe(self, ward, patientId, patientDiagnosis, patientTreatment,
                  timeNow=None):
        """This returns a tuple (careTier, patientTreatment)"""
        if patientDiagnosis.diagClassA == DiagClassA.WELL:
            if patientDiagnosis.overall == PatientOverallHealth.HEALTHY:
                return (CareTier.HOME, patientTreatment._replace(rehab=False))
            elif patientDiagnosis.overall == PatientOverallHealth.FRAIL:
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
        "Hospital description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = float(descr['meanPop'])
    agentList = []
    for i in xrange(int(round(meanPop))):
        ward = fac.manager.allocateAvailableBed(CareTier.HOME)
        assert ward is not None, 'Ran out of beds populating %(abbrev)s!' % descr
        a = PatientAgent('PatientAgent_HOME_%s_%d' % (ward._name, i), patch, ward)
        ward.lock(a)
        fac.handleIncomingMsg(pyrheabase.ArrivalMsg,
                              fac.getMsgPayload(pyrheabase.ArrivalMsg, a),
                              None)
        agentList.append(a)
    return agentList


def generateFull(facilityDescr, patch, policyClasses=None):
    fac = Community(facilityDescr, patch, policyClasses=policyClasses)
    return [fac], fac.getWards(), _populate(fac, facilityDescr, patch)


def estimateWork(facRec):
    if 'meanPop' in facRec:
        return facRec['meanPop'] / _constants['communityPatientCheckInterval']['value']
    elif 'nBeds' in facRec:
        return facRec['nBeds'] / _constants['communityPatientCheckInterval']['value']
    else:
        logger.warning('Cannot estimate work for %(abbrev)s' % facRec)
        return 0


def checkSchema(facilityDescr):
    global _validator
    if _validator is None:
        _validator = schemautils.getValidator(_schema)
    nErrors = sum([1 for e in _validator.iter_errors(facilityDescr)])  # @UnusedVariable
    return nErrors


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(os.path.join(os.path.dirname(__file__),
                                                      _constants_values),
                                         _constants_schema)
