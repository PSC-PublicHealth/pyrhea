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

import hospital

category = 'LTAC'
schema = 'facilityfacts_schema.yaml'


class LTAC(hospital.Hospital):
    pass


def generateFull(facilityDescr, patch):
    fac = LTAC(facilityDescr, patch)
    return [fac], fac.getWards(), hospital._populate(fac, facilityDescr, patch)


def estimateWork(descr):
    return hospital.estimateWork(descr)


def checkSchema(facilityDescr):
    return hospital.checkSchema(facilityDescr)
