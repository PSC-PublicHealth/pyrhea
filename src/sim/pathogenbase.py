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

from phacsl.utils.collections.phacollections import enum, namedtuple

PthStatus = enum('CLEAR', 'COLONIZED', 'CHRONIC', 'INFECTED', 'RECOVERED')
defaultPthStatus = PthStatus.CLEAR

class Pathogen(object):
    def __init__(self, ward, useWardCategory):
        """
        useWardCategory will typically be the same as that listed in the facility
        description file for the ward, but it may differ if category mapping has
        occurred.  For example, ward.fac.category might be 'SNF' while useWardCategory
        is 'NURSINGHOME' because that is the category definition in the class which
        actually implements 'SNF'.  This is intended to support category name mapping
        between facility descriptions and implementations.
        """
        self.ward = ward

    def getStatusChangeTree(self, patientStatus, careTier, treatment, startTime, timeNow):
        raise RuntimeError('Pathogen base class getStatusChangeTree was called.')

    def filterStatusChangeTrees(self, treeList, patientStatus, careTier, treatment, startTime, timeNow):
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
        raise RuntimeError('Pathogen base class getPatientPthCounts was called')
