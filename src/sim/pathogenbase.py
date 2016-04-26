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

from phacsl.utils.collections.phacollections import enum

Status = enum('CLEAR', 'COLONIZED', 'CHRONIC', 'INFECTED', 'RECOVERED')
defaultStatus = Status.CLEAR


class Pathogen(object):
    def __init__(self, ward):
        self.ward = ward

    def getStatusChangeTree(self, patientStatus, careTier, treatment, startTime, timeNow):
        raise RuntimeError('Pathogen base class getStatusChangeTree was called.')

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
