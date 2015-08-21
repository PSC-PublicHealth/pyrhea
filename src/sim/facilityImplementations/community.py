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
from itertools import cycle
import jsonschema
import yaml
import facilitybase

category = 'COMMUNITY'
schema = 'communityfacts_schema.yaml'
validator = None


class CommunityWard(facilitybase.Ward):
    """This 'ward' type represents being out in the community"""

    def __init__(self, name, patch, nBeds):
        facilitybase.Ward.__init__(self, name, patch, facilitybase.CareTier.HOME, nBeds)
        self.checkInterval = 30  # check health monthly
        #self.checkInterval = 1  # check health monthly


class Community(facilitybase.Facility):
    def __init__(self, descr, patch):
        facilitybase.Facility.__init__(self, '%(category)s_%(abbrev)s' % descr, patch)
        meanPop = descr['meanPop']
        nBeds = int(round(meanPop))
        self.addWard(CommunityWard('%s_%s_%s' % (category, patch.name, descr['abbrev']),
                                   patch, nBeds))


def _populate(fac, descr, patch):
    assert 'meanPop' in descr, \
        "Hospital description %(abbrev)s is missing the expected field 'meanPop'" % descr
    meanPop = float(descr['meanPop'])
    wards = fac.getWards()
    agentList = []
    for i, ward in zip(xrange(int(round(meanPop))), cycle(wards)):
        a = facilitybase.PatientAgent('PatientAgent_HOME_%s_%d' % (ward._name, i),
                                      patch, ward)
        ward.lock(a)
        agentList.append(a)
    return agentList


def generateFull(facilityDescr, patch):
    fac = Community(facilityDescr, patch)
    return [fac], fac.getWards(), _populate(fac, facilityDescr, patch)


def estimateWork(facRec):
    if 'meanPop' in facRec:
        return facRec['meanPop']
    elif 'nBeds' in facRec:
        return facRec['nBeds']
    else:
        print '####### Cannot estimate work for %(abbrev)s' % facRec
        return 0


def checkSchema(facilityDescr):
    global validator
    if validator is None:
        with open(os.path.join(os.path.dirname(__file__), schema), 'rU') as f:
            schemaJSON = yaml.safe_load(f)
        validator = jsonschema.validators.validator_for(schemaJSON)(schema=schemaJSON)
    nErrors = sum([1 for e in validator.iter_errors(facilityDescr)])  # @UnusedVariable
    return nErrors
