#! /usr/bin/env python

"""
Build a table for use by policyImplementations/transferbycapacity which where the capacity is
based on the difference between total discharges and total outgoing transfers.
"""

import json
import yaml
import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.formats.csv_tools as csv_tools
import sys

allKeySet, recs = yaml_tools.parse_all(sys.argv[1])
facDict = {r['abbrev']:r for r in recs}

newRecs = []
for rec in recs:
    try:
        entry = {'abbrev': rec['abbrev'], 'category': rec['category']}
        for key in ['totalDischarges', 'totalTransfersOut']:
            assert key in rec, "%s has no %s entry" % (rec['abbrev'], rec[key])
        tTO = sum([ent['count']['value'] for ent in rec['totalTransfersOut']])
        nCap = rec['totalDischarges']['value'] - tTO
        assert nCap > 0, '%s has deficit: %s vs %s' % (rec['abbrev'], rec['totalDischarges']['value'], tTO)
        if rec['category'] == 'HOSPITAL':
            if 'meanPop' in rec and 'meanPopICU' in rec:
                icuCap = int(round((float(rec['meanPopICU']['value']) / float(rec['meanPop']['value'])) * nCap))
            elif 'nBeds' in rec and 'nBedsICU' in rec:
                icuCap = int(round((float(rec['nBedsICU']['value']) / float(rec['nBeds']['value'])) * nCap))
            else:
                raise RuntimeError('Cannot estimate ICU capacity for %s' % rec['abbrev'])                            
            nCap -= icuCap
            entry['capacity'] = nCap
            entry['icuCapacity'] = icuCap
        else:
            entry['capacity'] = nCap
        newRecs.append(entry)
    except Exception, e:
        print 'Exception %s on %s' % (e, rec['abbrev'])

with open('discharges_which_are_not_transfers.yaml', 'w') as f:
    yaml.safe_dump(newRecs, f, default_flow_style=True, indent=4,
                   encoding='utf-8', width=130, explicit_start=True)
  