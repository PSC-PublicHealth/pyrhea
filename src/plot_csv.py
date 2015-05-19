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

_rhea_svn_id_="$Id$"

import csv
import os.path
import pyplot

import csv_tools


knownAbbrevs = set()
hospitalAbbrevs = set()
nhAbbrevs = set()

danteTestDir = '/home/welling/git/rhea-dante/test/nursing_home_CI_decolonization_2014'
lclDir = '/home/welling/workspace/pyRHEA/src'

"""
This one appears to contain an NxN matrix of transfer counts, where N is the number 
of facilities
"""
fName = '2007_OC_Nursing_Homes_and_Hospitals_Patient_Flow_UPS_RHEA_full-direct.csv'
with open(os.path.join(danteTestDir, fName), 'rU') as f:
    rdr = csv.reader(f)
    for idx, row in enumerate(rdr):
        print '%d: %d elements' % (idx, len(row))
#         if idx == 99:
#             print row

""" 
This one appears to contain an NxN matrix of transfer counts, where N is the number of facilities.
But this one is 'readmissions' rather than direct.  I don't understand why a transfer to a nursing 
home and then back doesn't just show up as two entries in the full-direct table.

The last two rows are bogus; the csv contains a stray entry which appears to be a column sum.
"""
fName = '2007_OC_Nursing_Homes_and_Hospitals_Patient_Flow_UPS_RHEA_full-readmissions.csv'
with open(os.path.join(danteTestDir, fName), 'rU') as f:
    rdr = csv.reader(f)
    for idx, row in enumerate(rdr):
        print '%d: %d elements' % (idx, len(row))
        if idx == 100:
            print 'should be empty'
            print row
        if idx == 101:
            print 'should have one entry, which is a misplaced column sum'
            print row

""" This contains facility information for all facilities """
fName = 'RHEA_Hosp-NH_Inputs_ADULT_2007_v03_06APR2012_properties_MRSA-STRAT-LOS_tweaked.csv'
with open(os.path.join(lclDir, fName), 'rU') as f:
    keys, recs = csv_tools.parseCSV(f)
abbrevKeys = []
for k in keys:
    if k.lower().endswith('_abbrev'):
        abbrevKeys.append(k)
infoRecDict = {}
for rec in recs:
    l = []
    for k in abbrevKeys:
        l.append(rec[k])
    firstAbbrev = l[0]
    assert all([a == firstAbbrev for a in l]), 'Mismatched abbreviation in %s' % l
    knownAbbrevs.add(firstAbbrev)
    infoRecDict[firstAbbrev] = rec
#print len(knownAbbrevs)

"""
This contains average daily census for each facility.  There are two columns;
the first is the place abbreviation and the second is the value as an integer.
Why are they integers and not floats?
"""
fName = 'Average_Daily_Census_2007_for_OC_Nursing_Homes_10-31-11.csv'
aveCensusDict = {}
with open(os.path.join(danteTestDir, fName), 'rU') as f:
    rdr = csv.reader(f)
    for row in rdr:
        if row[0].startswith('#'):
            continue
        assert len(row) == 2
        k = row[0].strip()
        assert k not in aveCensusDict, 'duplicate row for %s in %s' % (k, fName)
        aveCensusDict[k] = int(row[1])
# print aveCensusDict
for k in aveCensusDict.keys():
    if k not in infoRecDict:
        print 'no infoRec for nursing home <%s> with census %s' % (k, aveCensusDict[k])
    if k not in knownAbbrevs:
        print 'found an abbrev late: <%s>' % k
        knownAbbrevs.add(k)
    nhAbbrevs.add(k)

"""
Is this one actually input data? Were these entries derived from some larger collection of
patient LOS information?  It's a CSV format with row descriptions in comments.

The beta value info in the bottom two sections is not used.  The county value (ORNG) is
essentially the defaults; BWH just has an anomalous amount of info because they are a
research hospital.
"""
fName = 'ICU_Estimated_LOS_7-22-09.txt'

"""
This one is a two-column CSV with the first column being an abbrev and the second an integer,
presumably the length of stay of one individual in days.
"""
fName = 'Length_of_Stay_2007_to_2009_OC_Nursing_Homes-12-18-11_SMB_with_abbrev_RHEA.csv'
nhLOSDict = {}
with open(os.path.join(danteTestDir, fName), 'rU') as f:
    rdr = csv.reader(f)
    for row in rdr:
        if row[0].startswith('#'):
            continue
        assert len(row) == 2
        k = row[0].strip()
        if k not in nhLOSDict:
            nhLOSDict[k] = []
        nhLOSDict[k].append(int(row[1]))
# for k, v in nhLOSDict.items():
#     print '%s : %d entries' % (k, len(v))
for k in nhLOSDict.keys():
    if k not in nhAbbrevs:
        print 'Found a nursing home abbrev late: <%s>' % k
        nhAbbrevs.add(k)
for k in nhAbbrevs:
    if k not in nhLOSDict:
        print 'No LOS records for nursing home <%s>' % k

"""
If you get sent to a hospital from an assisted living nursing home, they hold your bed.
This is a histogram of lengths of time people were away from their nursing homes.
Any stay longer than 14 days appears only as a surplus in the 'All LOS' bin.
"""
fName = 'Mean_Temporary_Discharge_LOS_for_2007_to_2009_2-2-12.csv'
nhLOSHistoDict = {}
with open(os.path.join(danteTestDir, fName), 'rU') as f:
    keys, recs = csv_tools.parseCSV(f)
    assert set(keys) == set([u'# LTC', u'LOS = 0d', u'LOS = 1d', u'LOS = 2d',
                             u'LOS = 3d', u'LOS = 4d', u'LOS = 5d', u'LOS = 6d',
                             u'LOS = 7', u'LOS = 8d', u'LOS = 9d', u'LOS = 10d',
                             u'LOS = 11d', u'LOS = 12d', u'LOS = 13d', u'LOS = 14d',
                             u'All LOS (N)'])
    for rec in recs:
        assert rec[u'# LTC'] in nhAbbrevs
        hR = {}
        for i in xrange(15):
            try:
                hR[i] = rec[u'LOS = %dd' % i]
            except KeyError:
                hR[i] = rec[u'LOS = %d' % i]
        hR['all'] = rec[u'All LOS (N)']
        nhLOSHistoDict[rec[u'# LTC']] = hR
# for k,v in nhLOSHistoDict.items():
#     print '%s: %s' % (k, v)
for k in nhAbbrevs:
    if k not in nhLOSHistoDict:
        print 'No histogram for nursing home <%s>' % k

"""
This is a 3 column table for nursing homes giving the number of patients re-entering after
hospital visits.  Column 0 is abbreviation, 1 is three-year sum, 2 is rounded average.
"""
fName = 'Re-entries_2007-2009_OC_Nursing_Homes_1-4-12_RHEA_CRM.csv'
