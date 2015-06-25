#! /usr/bin/env python

import yaml
import os.path
from collections import OrderedDict

import yaml_ordered
yaml_ordered.install()


def unicode_safe_constructor(loader, node):
    return node.value

yaml.SafeLoader.add_constructor("tag:yaml.org,2002:python/unicode",
                                unicode_safe_constructor)


def parse_all(dirName):
    recs = []
    allKeys = set()
    for nm in os.listdir(dirName):
        if nm.endswith('.yaml'):
            with open(os.path.join(dirName, nm), 'r') as f:
                newD = yaml.safe_load(f)
                allKeys.update(newD.keys())
                recs.append(newD)
    return allKeys, recs


def save_all(dirName, recList):
    # Re-copy everything, but with ordered keys
    theseComeFirst = ['name', 'abbrev', 'category']
    theseComeFirstSet = set(theseComeFirst)
    newRecs = []
    for rec in recList:
        newRec = OrderedDict()
        for k in theseComeFirst:
            newRec[k] = rec[k]
        for k in rec:
            if k not in theseComeFirstSet:
                newRec[k] = rec[k]
        newRecs.append(newRec)

    # Write output
    noAbbrevCtr = 0
    for rec in newRecs:
        if 'abbrev' in rec and len(rec['abbrev']) > 0:
            ofname = rec['abbrev'] + '.yaml'
        else:
            ofname = 'none%d.yaml' % noAbbrevCtr
            noAbbrevCtr += 1
        with open(os.path.join(dirName, ofname), 'w') as f:
            yaml.safe_dump(rec, f,
                           default_flow_style=False, indent=4,
                           encoding='utf-8', width=130, explicit_start=True)
