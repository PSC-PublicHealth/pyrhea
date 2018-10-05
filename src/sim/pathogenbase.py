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

from phacsl.utils.collections.phacollections import enum, namedtuple
from freezerbase import FreezerError

PthStatus = enum('CLEAR', 'COLONIZED', 'CHRONIC', 'INFECTED', 'RECOVERED', 'UNDETCOLONIZED')
defaultPthStatus = PthStatus.CLEAR

class Pathogen(object):
    def __init__(self, ward, implCategory):
        """
        implCategory will typically be the same as that listed in the facility
        description file for the ward, but it may differ if category mapping has
        occurred.  For example, ward.fac.category might be 'SNF' while implCategory
        is 'NURSINGHOME' because that is the category definition in the class which
        actually implements 'SNF'.  This is intended to support category name mapping
        between facility descriptions and implementations.
        """
        self.ward = ward
        self.patientPth = self._emptyPatientPth()
        self.patientPthTime = None

    def flushCaches(self):
        """
        Derived classes often cache things like BayesTrees, but the items in the cache
        can become invalid when a new scenario starts and the odds of transitions are
        changed.  This method is called when the environment wants to trigger a cache
        flush.
        """
        pass

    def _emptyPatientPth(self):
        return {k:0 for k in PthStatus.names.keys()}

    def getStatusChangeTree(self, patientAgent, modifierDct, startTime, timeNow):
        raise RuntimeError('Pathogen base class getStatusChangeTree was called.')

    def filterStatusChangeTrees(self, treeList, patientAgent, startTime, timeNow):
        return treeList[:]

    def initializePatientState(self, patient):
        """
        This method assigns the patient a Status appropriate to time 0- that is,
        it implements the initial seeding of the patient population with pathogen.
        """
        raise RuntimeError('Pathogen base class initializePatientState was called')

    def getPatientPthCounts(self, timeNow):
        """
        Returns a dict of the form {PthStatus.CLEAR : nClearPatients, ...}
        """
        if self.patientPthTime != timeNow:
            dct = self._emptyPatientPth()
            try:
                for pt in self.ward.getPatientList():
                    dct[pt._status.pthStatus] += 1
            except FreezerError:
                pass  # We're going to have to ignore freeze-dried patients
            self.patientPth = dct
            self.patientPthTime = timeNow
        return self.patientPth
