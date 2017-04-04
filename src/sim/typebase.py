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
from quilt.netinterface import GblAddr
from pathogenbase import PthStatus

CareTier = enum('HOME', 'NURSING', 'LTAC', 'HOSP', 'ICU', 'VENT', 'SKILNRS')

PatientOverallHealth = enum('HEALTHY', 'FRAIL')

DiagClassA = enum('HEALTHY', 'NEEDSREHAB', 'NEEDSLTAC', 'SICK', 'VERYSICK', 'DEATH',
                  'NEEDSVENT', 'NEEDSSKILNRS')

TreatmentProtocol = namedtuple('TreatmentProtocol',
                               ['rehab',
                                'contactPrecautions',
                                'creBundle'
                                ],
                               field_types=[bool, bool, bool])

TREATMENT_DEFAULT = TreatmentProtocol(rehab=False, contactPrecautions=False, creBundle=False)

PatientStatus = namedtuple('PatientStatus',
                           ['overall',              # one of PatientOverallHealth
                            'diagClassA',           # one of DiagClassA
                            'startDateA',           # date diagClassA status was entered
                            'pthStatus',            # one of PthStatus
                            'startDatePth',         # date PthStatus status was entered
                            'relocateFlag',         # true if patient needs relocation
                            'justArrived',          # true on patient's first day in new location
                            'canClear',             # true if patient can spontaneously clear infection
                            'homeAddr'              # GblAddr of patient's home tract or NH
                            ],
                           field_types=[PatientOverallHealth, DiagClassA, None, PthStatus, None,
                                        bool, bool, bool, GblAddr])

PatientDiagnosis = namedtuple('PatientDiagnosis',
                              ['overall',              # one of PatientOverallHealth
                               'diagClassA',           # one of DiagClassA
                               'startDateA',           # Date diagClassA was entered
                               'pthStatus',            # one of PthStatus
                               'relocateFlag'          # true if patient needs relocation
                               ],
                              field_types=[PatientOverallHealth, DiagClassA, None, PthStatus, bool])


