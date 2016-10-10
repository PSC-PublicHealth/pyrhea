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
import logging
import math
import pyrheautils
from phacsl.utils.collections.phacollections import SingletonMetaClass
from stats import CachedCDFGenerator, BayesTree, fullCRVFromPDFModel
from facilitybase import CareTier, PatientOverallHealth, DiagClassA, TreatmentProtocol
from facilitybase import PatientStatusSetter, Pathogen, PthStatus, PthStatusSetter

pathogenName = 'MRSA'
_constants_values = 'mrsa_constants.yaml'
_constants_schema = 'mrsa_constants_schema.yaml'
_constants = None

logger = logging.getLogger(__name__)


def _valFromCategoryEntry(key, ctg, constantsJSON):
    if key not in constantsJSON:
        raise RuntimeError('Constant list for %s was not found' % key)
    for v in constantsJSON[key]:
        if v['category'] == ctg:
            return v['frac']['value']
    raise RuntimeError('Constant entry for category %s was not found' % ctg)


class MRSACore(object):
    """This is where we put things that are best shared across all MRSA instances"""
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        colCRV = fullCRVFromPDFModel(_constants['colonizationDurationPDF'])
        self.colCachedCDF = CachedCDFGenerator(colCRV)
        self.colTreeCache = {}
        infCRV = fullCRVFromPDFModel(_constants['infectionDurationPDF'])
        self.infCachedCDF = CachedCDFGenerator(infCRV)
        self.infTreeCache = {}

        self.exposureTreeCache = {}  # so simple we don't need a CachedCDFGenerator


class MRSA(Pathogen):
    def __init__(self, ward):
        super(MRSA, self).__init__(ward)
        self.core = MRSACore()
        self.patientPth = self._emptyPatientPth()
        self.patientPthTime = None
        initialFractionColonized = _valFromCategoryEntry('initialFractionColonized',
                                                         self.ward.fac.category,
                                                         _constants)
        initialFractionInfected = _valFromCategoryEntry('initialFractionInfected',
                                                        self.ward.fac.category,
                                                        _constants)
        colTree = BayesTree(PthStatusSetter(PthStatus.CHRONIC),
                            PthStatusSetter(PthStatus.COLONIZED),
                            _constants['fractionColonizedChronic']['value'])
        fracClear = 1.0 - (initialFractionColonized + initialFractionInfected)
        self.initializationBayesTree = BayesTree.fromLinearCDF([(initialFractionColonized,
                                                                 colTree),
                                                                (initialFractionInfected,
                                                                 PthStatusSetter(PthStatus.INFECTED)),
                                                                (fracClear,
                                                                 PatientStatusSetter())])

    def _emptyPatientPth(self):
        return {k:0 for k in PthStatus.names.keys()}

    def __str__(self):
        return '<%s>' % pathogenName

    def initializePatientState(self, patient):
        """
        This method assigns the patient a Status appropriate to time 0- that is,
        it implements the initial seeding of the patient population with pathogen.
        """
        patient._status = self.initializationBayesTree.traverse().set(patient._status, 0)

    def getPatientPthCounts(self, timeNow):
        """
        Returns a dict of the form {PthStatus.CLEAR : nClearPatients, ...}
        """
        if self.patientPthTime != timeNow:
            d = self._emptyPatientPth()
            for p in self.ward.getPatientList():
                d[p._status.pthStatus] += 1
            self.patientPth = d
            self.patientPthTime = timeNow
        return self.patientPth

    def getStatusChangeTree(self, patientStatus, careTier, treatment, startTime, timeNow):
        if patientStatus.pthStatus == PthStatus.CLEAR:
            pP = self.getPatientPthCounts(timeNow)
            nExposures = pP[PthStatus.COLONIZED] + pP[PthStatus.INFECTED] + pP[PthStatus.CHRONIC]
            nTot = sum(pP.values())
            dT = timeNow - startTime
            key = (nExposures, nTot, dT)
            if key in self.core.exposureTreeCache:
                return self.core.exposureTreeCache[key]
            else:
                pSafe = math.pow((1.0 - _constants['beta']['value']), dT)
                tree = BayesTree(PatientStatusSetter(),
                                 BayesTree(PthStatusSetter(PthStatus.CHRONIC),
                                           PthStatusSetter(PthStatus.COLONIZED),
                                           _constants['fractionColonizedChronic']['value']),
                                 pSafe)
                self.core.exposureTreeCache[key] = tree
                return tree
        else:
            key = (startTime - patientStatus.startDatePth, timeNow - patientStatus.startDatePth)
            if patientStatus.pthStatus == PthStatus.COLONIZED:
                if key in self.core.colTreeCache:
                    return self.core.colTreeCache[key]
                else:
                    changeProb = self.core.colCachedCDF.intervalProb(*key)
                    infFrac = _constants['colonizationToInfectionFrac']['value']
                    tree = BayesTree(BayesTree(PthStatusSetter(PthStatus.INFECTED),
                                               PthStatusSetter(PthStatus.CLEAR),
                                               infFrac),
                                     PatientStatusSetter(),
                                     changeProb)
                    self.core.colTreeCache[key] = tree
                    return tree
            elif patientStatus.pthStatus == PthStatus.INFECTED:
                if key in self.core.infTreeCache:
                    return self.core.infTreeCache[key]
                else:
                    changeProb = self.core.infCachedCDF.intervalProb(*key)
                    tree = BayesTree(PthStatusSetter(PthStatus.COLONIZED),
                                     PatientStatusSetter(),
                                     changeProb)
                    self.core.infTreeCache[key] = tree
                    return tree
            elif patientStatus.pthStatus == PthStatus.CHRONIC:
                    return BayesTree(PatientStatusSetter())
            else:
                raise RuntimeError('Unexpected %s status %s at %s' %
                                   (self, patientStatus.pthStatus, self.ward._name))


def getPathogenClass():
    return MRSA


###########
# Initialize the module
###########
_constants = pyrheautils.importConstants(os.path.join(os.path.dirname(__file__),
                                                      _constants_values),
                                         _constants_schema)
