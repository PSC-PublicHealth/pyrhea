#! /usr/bin/env python

"""
facilityfacts mods:  These are location-specific tweaks
- SJMC doesn't actually need nBeds
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
        if 'nBeds' in rec:
            cap = rec['nBeds']['value']
        else:
            cap = rec['meanPop']['value']
        assert cap > 0.0, 'cap for %s is non-positive' % entry['abbrev']
        baseCap = cap
        for oKey, nKey in [('fracAdultPatientDaysICU', 'icuCapacity'),
                           ('fracVentBeds', 'ventCapacity'),
                           ('fracSkilledBeds', 'skilNrsCapacity')]:
            if oKey in rec and rec[oKey]['value'] > 0.0:
                delta = rec[oKey]['value'] * baseCap
                cap -= delta
                entry[nKey] = delta
        assert cap >= 0.0, ('%s: Not enough capacity to cover bed categories'
                            % entry['abbrev'])
        entry['capacity'] = cap
        newRecs.append(entry)
    except Exception, e:
        print 'Exception %s on %s' % (e, rec['abbrev'])

with open('bedcapacity.yaml', 'w') as f:
    yaml.safe_dump(newRecs, f, default_flow_style=True, indent=4,
                   encoding='utf-8', width=130, explicit_start=True)
  