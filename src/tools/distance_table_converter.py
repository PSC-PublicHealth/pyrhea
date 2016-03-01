#! /usr/bin/env python

import yaml
import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.formats.csv_tools as csv_tools

with open('../../models/OrangeCounty2013/transitmatrix.csv', 'rU') as f:
    keys, recs = csv_tools.parseCSV(f)

facKeys, facRecs = yaml_tools.parse_all('../../models/OrangeCounty2013/facilityfactsCurrent2013')
facDict = {r['abbrev']: r for r in facRecs}

srcDir = {}
for rec in recs:
    src = rec['From']
    dest = rec['To']
    secs = rec['Seconds']
    meters = rec['Meters']
    cat = facDict[dest]['category']
    if src not in srcDir:
        srcDir[src] = {}
    if cat not in srcDir[src]:
        srcDir[src][cat] = {}
    assert dest not in srcDir[src][cat], 'Double entry for %s %s?' % (src, dest)
    srcDir[src][cat][dest] = {'seconds': secs, 'meters': meters}

with open('../../models/OrangeCounty2013/transitmatrix.yaml', 'w') as f:
    yaml.safe_dump(srcDir, f,
                   default_flow_style=True, indent=4,
                   encoding='utf-8', width=130, explicit_start=True)
