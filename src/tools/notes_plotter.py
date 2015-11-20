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

import yaml
import os.path
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

from scipy.stats import lognorm, expon


nameMap = {'NURSINGHOME': 'NURSING',
           'HOSPITAL': 'HOSP',
           'COMMUNITY': 'HOME',
           'LTAC': 'HOSP'}

careTiers = ['HOME', 'NURSING', 'HOSP', 'ICU']


defaultNotesPath = '/home/welling/workspace/pyRHEA/src/sim/notes.pkl'
defaultDescrDirs = ['/home/welling/workspace/pyRHEA/models/OrangeCounty/facilityfactsCurrent',
                    '/home/welling/workspace/pyRHEA/models/OrangeCounty/communitiesCurrent']
defaultConstDir = '/home/welling/workspace/pyRHEA/src/sim/facilityImplementations'


def fullPDFFromMeanLOS(fitParms):
    def losFun(xVec):
        xArr = np.asarray(xVec, np.float64)
        mean = fitParms[0]
        sigma = fitParms[1]
        mu = math.log(mean) - (0.5 * sigma * sigma)
        print "fullPDF: sigma= %s, mu= %s" % (sigma, mu)
        pVec = lognorm.pdf(xArr, sigma, scale=math.exp(mu), loc=0.0)
        return pVec
    return losFun


def fullPDFFromLOSModel(losModel):
    if losModel['pdf'] == 'lognorm(mu=$0,sigma=$1)':
        def losFun(xVec):
            xArr = np.asarray(xVec, np.float64)
            mu, sigma = losModel['parms']
            print "fullPDF: sigma= %s, mu= %s" % (sigma, mu)
            pVec = lognorm.pdf(xArr, sigma, scale=math.exp(mu), loc=0.0)
            return pVec
    elif losModel['pdf'] == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)':
        def losFun(xVec):
            xArr = np.asarray(xVec, np.float64)
            k, mu, sigma, lmda = losModel['parms']
            print "fullPDF: k= %s, sigma= %s, mu= %s, lmda= %s" % (k, sigma, mu, lmda)
            pVec = ((k * lognorm.pdf(xArr, sigma, scale=math.exp(mu), loc=0.0)
                     + ((1.0-k) * expon.pdf(xArr, scale=1.0/lmda))))
            return pVec
    elif losModel['pdf'] == 'expon(lambda=$0)':
        def losFun(xVec):
            xArr = np.asarray(xVec, np.float64)
            lmda = losModel['parms'][0]
            print "fullPDF: lmda= %s" % lmda
            pVec = expon.pdf(xArr, scale=1.0/lmda)
            return pVec
    else:
        raise RuntimeError('Unknown LOS model %s' % losModel['pdf'])
    return losFun


class LOSPlotter(object):
    def __init__(self, descr, constants):
        self.fullPDFs = {}
        if descr['category'] == 'HOSPITAL':
            self.fullPDFs['ICU'] = fullPDFFromMeanLOS([descr['meanLOSICU'],
                                                       constants['icuLOSLogNormSigma']])
            self.fullPDFs['HOSP'] = fullPDFFromLOSModel(descr['losModel'])
        elif descr['category'] == 'LTAC':
            self.fullPDFs['HOSP'] = fullPDFFromLOSModel(descr['losModel'])
        elif descr['category'] == 'NURSINGHOME':
            self.fullPDFs['NURSING'] = fullPDFFromLOSModel(descr['losModel'])
        elif descr['category'] == 'COMMUNITY':
            self.fullPDFs['HOME'] = fullPDFFromLOSModel(constants['communityLOSModel'])
        else:
            raise RuntimeError("%s has unknown category %s - cannot plot LOS"
                               % (abbrev, descr['category']))

    def plot(self, tier, axes, nBins, rMin, rMax, scale, pattern='r-'):
        if tier in self.fullPDFs:
            curveX = np.linspace(rMin, rMax, nBins)
#             curveY = self.fullPDFs[tier](curveX) * scale * ((rMax-rMin)/nBins)
            curveY = self.fullPDFs[tier](curveX) * scale
            axes.plot(curveX, curveY, pattern, lw=2, alpha=0.6)


def loadFacilityDescription(abbrev):
    for descDir in defaultDescrDirs:
        fname = os.path.join(descDir, abbrev + '.yaml')
        if os.path.exists(fname):
            with open(fname, 'rU') as f:
                return yaml_tools._simplify(yaml.safe_load(f))
    raise RuntimeError('No description file found for %s' % abbrev)


def loadFacilityTypeConstants(category):
    with open(os.path.join(defaultConstDir, category.lower() + '_constants.yaml'), 'rU') as f:
        return yaml_tools._simplify(yaml.safe_load(f))


def importNotes(fname):
    with open(fname, 'r') as f:
        stuff = pickle.load(f)
    return stuff


def overallLOSFig(catNames, allOfCategoryDict):
    figs1, axes1 = plt.subplots(nrows=len(allOfCategoryDict), ncols=1)
    for offset, cat in enumerate(catNames):
        bins = []
        counts = []
        for tier in careTiers:
            try:
                for k, v in allOfCategoryDict[cat][tier + '_LOS'].histogram().items():
                    bins.append(k)
                    counts.append(v)
            except:
                pass

        rects = axes1[offset].bar(bins, counts, color='r')  # @UnusedVariable
        axes1[offset].set_ylabel('Counts')
        axes1[offset].set_xlabel('Days')
        axes1[offset].set_title('LOS for ' + cat)
    figs1.tight_layout()
    figs1.canvas.set_window_title("LOS Histograms By Category")


def singleLOSFig(abbrev, notesDict):
    figs1b, axes1b = plt.subplots(nrows=len(careTiers), ncols=1)
    losDescr = loadFacilityDescription(abbrev)
    constants = loadFacilityTypeConstants(losDescr['category'])
    losPlotter = LOSPlotter(losDescr, constants)
    for k in notesDict.keys():
        if k.endswith(abbrev):
            for idx, tier in enumerate(careTiers):
                bins = []
                counts = []
                tK = tier + '_LOS'
                if tK in notesDict[k]:
                    print 'hit %s for %s' % (tier, abbrev)
                    print 'dict: %s' % notesDict[k][tK].d
                    pairList = notesDict[k][tK].histogram().items()
                    pairList.sort()
                    print 'items: %s' % {k:v for k,v in pairList}
                    quantum = notesDict[k][tK].d['quantum']
                    barWidth = 1.6 * quantum
                    for b, c in pairList:
                        bins.append(b - 0.5*(quantum + barWidth))
                        counts.append(c)
                    rects = axes1b[idx].bar(bins, counts, width=barWidth,  # @UnusedVariable
                                            color='r')
                else:
                    rects = axes1b[idx].bar([], [], color='r')  # @UnusedVariable
                axes1b[idx].set_ylabel('Counts')
                axes1b[idx].set_xlabel('Days')
                axes1b[idx].set_title('LOS for ' + tier)
                if bins:
                    losPlotter.plot(tier, axes1b[idx], 300, 0.0, max(bins),
                                    sum(counts))
            break
    figs1b.tight_layout()
    figs1b.canvas.set_window_title("%s LOS Histograms By Category" % abbrev)


def bedBounceFig(allOfCategoryDict):
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
        rects = axes2[offset].bar(bins, counts, color='r')  # @UnusedVariable
        axes2[offset].set_ylabel('Counts')
        axes2[offset].set_xlabel('Bounces')
        axes2[offset].set_title('Bounces for Bed Requests for ' + tier)
    # figs2.tight_layout()
    figs2.canvas.set_window_title("Bed Request Bounce Histograms")


def patientFlowFig(allOfCategoryDict):
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
        rects = axes3[offset].bar([0, 1], [nArrive, nDepart], color='r')  # @UnusedVariable
        axes3[offset].set_ylabel('Counts')
        axes3[offset].set_title(tier)
        axes3[offset].set_xticks([0.5, 1.5])
        axes3[offset].set_xticklabels(['in', 'out'])
    figs3.tight_layout()
    figs3.canvas.set_window_title("Patient Flow By Tier")


def patientFateFig(catNames, allOfCategoryDict):
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


def occupancyTimeFig(specialDict):
    figs5, axes5 = plt.subplots(nrows=1, ncols=len(specialDict))
    if len(specialDict) == 1:
        axes5 = [axes5]
    for offset, (patchName, data) in enumerate(specialDict.items()):
        try:
            occDataList = data['occupancy']
            assert isinstance(occDataList, types.ListType), \
                'Special data %s is not a list' % patchName
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


def countBirthsDeaths(catNames, allOfCategoryDict):
    totBirths = 0
    totDeaths = 0
    for offset, cat in enumerate(catNames):  # @UnusedVariable
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
        totBirths += nBirths
        totDeaths += nDeaths
    print 'totals: %d births, %d deaths' % (totBirths, totDeaths)


if len(sys.argv) > 1:
    notesFName = sys.argv[1]
else:
    notesFName = defaultNotesPath

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

catNames = allOfCategoryDict.keys()[:]

countBirthsDeaths(catNames, allOfCategoryDict)

overallLOSFig(catNames, allOfCategoryDict)

singleLOSFig('SJUD', notesDict)
singleLOSFig('WAEC', notesDict)
singleLOSFig('CM69', notesDict)
# bedBounceFig(allOfCategoryDict)
patientFlowFig(allOfCategoryDict)
patientFateFig(catNames, allOfCategoryDict)
occupancyTimeFig(specialDict)

plt.show()
