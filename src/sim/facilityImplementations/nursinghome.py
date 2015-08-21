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

import facilitybase
from hospital import checkSchema as hospitalCheckSchema, estimateWork as hospitalEstimateWork

category = 'NURSINGHOME'
schema = 'facilityfacts_schema.yaml'


# class CommunityWard(Ward):
#     """This 'ward' type represents being out in the community"""
#
#     def __init__(self, name, patch, nBeds):
#         Ward.__init__(self, name, patch, CareTier.HOME, nBeds)
#         self.checkInterval = 30  # check health monthly
#         #self.checkInterval = 1  # check health monthly


class NursingHome(facilitybase.Facility):
    def __init__(self, descr, patch):
        facilitybase.Facility.__init__(self, '%(category)s_%(abbrev)s' % descr, patch)
        assert 'nBeds' in descr, 'Nursing home %(abbrev) description is missing nBeds' % descr
        nBeds = int(descr['nBeds'])
        self.addWard(facilitybase.Ward('%s_%s_%s' % (category, patch.name, descr['abbrev']),
                                       patch, facilitybase.CareTier.NURSING, nBeds))


def _populate(fac, descr, patch):
    assert 'meanPop' in descr, \
        "Nursing home description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = descr['meanPop']
    if meanPop > descr['nBeds']:
        print ('####### Nursing Home %(abbrev)s meanPop %(meanPop)s > nBeds %(nBeds)s #########'
               % descr)
        meanPop = descr['nBeds']
    wards = fac.getWards()
    agentList = []
    for i, ward in zip(xrange(int(round(meanPop))), cycle(wards)):
        a = facilitybase.PatientAgent('PatientAgent_NURSING_%s_%d' % (ward._name, i), patch, ward)
        ward.lock(a)
        agentList.append(a)
    return agentList


def generateFull(facilityDescr, patch):
    fac = NursingHome(facilityDescr, patch)
    return [fac], fac.getWards(), _populate(fac, facilityDescr, patch)


def estimateWork(descr):
    return hospitalEstimateWork(descr)


def checkSchema(facilityDescr):
    return hospitalCheckSchema(facilityDescr)
