#! /usr/bin/env python

import yaml
import os.path
from collections import OrderedDict
import types

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


def _simplify(entry):
    if isinstance(entry, types.DictType):
        pairs = []
        for k, v in entry.items():
            if k.find('_prov') >= 0:
                continue
            elif k == 'prov':
                continue
            elif k == 'value':
                return v
            else:
                simpV = _simplify(v)
                if isinstance(simpV, types.ListType):
                    newList = []
                    for item in simpV:
                        if (isinstance(item, types.DictType)
                                and 'category' in item):
                            if 'count' in item:
                                pairs.append((k + item['category'], item['count']))
                            elif 'value' in item:
                                pairs.append((k + item['category'], item['value']))
                        else:
                            newList.append(item)
                    if newList:
                        pairs.append((k, newList))
                else:
                    pairs.append((k, simpV))
        return dict(pairs)
    elif isinstance(entry, types.ListType):
        return [_simplify(e) for e in entry]
    else:
        return entry


def parse_all_simplified(dirName):
    allKeys, rawRecs = parse_all(dirName)
    return allKeys, [_simplify(r) for r in rawRecs]


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
