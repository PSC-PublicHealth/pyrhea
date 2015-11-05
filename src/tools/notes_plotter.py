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

# import csv_tools
# import yaml_tools
# import noteholder
# from util import HistoVal
import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.notes.noteholder as noteholder
from phacsl.utils.notes.statval import HistoVal
import sys
import math
import pickle
import types

import numpy as np
import matplotlib.pyplot as plt

from scipy.cluster.vq import whiten, kmeans, vq


def importNotes(fname):
    with open(fname, 'r') as f:
        stuff = pickle.load(f)
    return stuff

if len(sys.argv) > 1:
    notesFName = sys.argv[1]
else:
    notesFName = '/home/welling/workspace/pyRHEA/src/sim/notes.pkl'

notesDict = importNotes(notesFName)
categoryDict = {}
specialDict = {}
for nm, d in notesDict.items():
    try:
        category, abbrev = tuple(nm.split('_'))
        abbrev = abbrev.lower()
        if category not in categoryDict:
            categoryDict[category] = {}
        categoryDict[category][abbrev] = d
    except:
        specialDict[nm] = d

noteHolderGroup = noteholder.NoteHolderGroup()
allOfCategoryDict = {}
for category in categoryDict.keys():
    allOfCategoryDict[category] = noteHolderGroup.createNoteHolder()
    for abbrev, nhDict in categoryDict[category].items():
        allOfCategoryDict[category].addNote({k: v for k, v in nhDict.items() if k != 'name'})

figs1, axes1 = plt.subplots(nrows=len(allOfCategoryDict), ncols=1)

catNames = allOfCategoryDict.keys()[:]
nameMap = {'NURSINGHOME': 'NURSING',
           'HOSPITAL': 'HOSP',
           'COMMUNITY': 'HOME',
           'LTAC': 'HOSP'}

careTiers = ['HOME', 'NURSING', 'HOSP', 'ICU']

for offset, cat in enumerate(catNames):
    # print '%s: %s' % (cat, allOfCategoryDict[cat].keys())
    if 'death' in allOfCategoryDict[cat]:
        nDeaths = allOfCategoryDict[cat]['death']
    else:
        nDeaths = 0
    if 'births' in allOfCategoryDict[cat]:
        nBirths = allOfCategoryDict[cat]['births']
    else:
        nBirths = 0
    print '%s: %d births, %d deaths' % (cat, nBirths, nDeaths)
    bins = []
    counts = []
    for tier in careTiers:
        try:
            for k, v in allOfCategoryDict[cat][tier + '_LOS'].histogram().items():
                bins.append(k)
                counts.append(v)
        except Exception, e:
            pass

    rects = axes1[offset].bar(bins, counts, color='r')
    axes1[offset].set_ylabel('Counts')
    axes1[offset].set_xlabel('Days')
    axes1[offset].set_title('LOS for ' + cat)

figs1.tight_layout()
figs1.canvas.set_window_title("LOS Histograms By Category")

figs2, axes2 = plt.subplots(nrows=len(careTiers), ncols=1)
for offset, tier in enumerate(careTiers):
    key = '%s_bounce_histo' % tier
    hv = HistoVal([])
    for cat in catNames:
        if key in allOfCategoryDict[cat]:
            hv += allOfCategoryDict[cat][key]
    bins = []
    counts = []
    for k, v in hv.histogram().items():
        bins.append(k)
        counts.append(v)
    rects = axes2[offset].bar(bins, counts, color='r')
    print bins
    print counts
    axes2[offset].set_ylabel('Counts')
    axes2[offset].set_xlabel('Bounces')
    axes2[offset].set_title('Bounces for Bed Requests for ' + tier)
#figs2.tight_layout()
figs2.canvas.set_window_title("Bed Request Bounce Histograms")

figs3, axes3 = plt.subplots(nrows=1, ncols=len(careTiers))
for offset, tier in enumerate(careTiers):
    nArrive = 0
    nDepart = 0
    for cat in catNames:
        try:
            nArrive += allOfCategoryDict[cat]['%s_arrivals' % tier]
        except:
            pass
        try:
            nDepart += allOfCategoryDict[cat]['%s_departures' % tier]
        except:
            pass
    rects = axes3[offset].bar([0, 1], [nArrive, nDepart], color='r')
    axes3[offset].set_ylabel('Counts')
    axes3[offset].set_title(tier)
    axes3[offset].set_xticks([0.5, 1.5])
    axes3[offset].set_xticklabels(['in', 'out'])
figs3.tight_layout()
figs3.canvas.set_window_title("Patient Flow By Tier")

figs4, axes4 = plt.subplots(nrows=1, ncols=len(catNames))
for offset, cat in enumerate(catNames):
    keys = ['death'] + ['%s_found' % tier for tier in careTiers]
    labels = ['death'] + careTiers
    counts = []
    for k in keys:
        if k in allOfCategoryDict[cat]:
            counts.append(allOfCategoryDict[cat][k])
        else:
            counts.append(0)
    ct2 = []
    lbl2 = []
    for ct, lbl in zip(counts, labels):
        if ct != 0:
            ct2.append(ct)
            lbl2.append(lbl)
    axes4[offset].pie(ct2, labels=lbl2, autopct='%1.1f%%', startangle=90)
    axes4[offset].axis('equal')
    axes4[offset].set_title(cat)
figs4.tight_layout()
figs4.canvas.set_window_title("Patient Fates")

figs5, axes5 = plt.subplots(nrows=1, ncols=len(specialDict))
if len(specialDict) == 1:
    axes5 = [axes5]
for offset, (patchName, data) in enumerate(specialDict.items()):
    try:
        occDataList = data['occupancy']
        assert isinstance(occDataList, types.ListType), 'Special data %s is not a list' % patchName
        fields = {}
        for d in occDataList:
            for k, v in d.items():
                if k not in fields:
                    fields[k] = []
                fields[k].append(v)
        assert 'day' in fields, 'Date field is missing for special data %s' % patchName
        dayList = fields['day']
        del fields['day']
        keys = fields.keys()
        keys.sort()
        for k in keys:
            l = fields[k]
            assert len(l) == len(dayList), ('field %s is the wrong length in special data %s'
                                            % patchName)
            axes5[offset].plot(dayList, l, label=k)
        axes5[offset].set_xlabel('Days')
        axes5[offset].set_ylabel('Occupancy')
        axes5[offset].legend()
    except Exception, e:
        print e
    axes5[offset].set_title(patchName)
figs5.tight_layout()
figs5.canvas.set_window_title("Time History of Occupancy")

plt.show()
