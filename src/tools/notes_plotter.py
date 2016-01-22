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
import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.notes.noteholder as noteholder
from phacsl.utils.notes.statval import HistoVal
from stats import lognormplusexp
import pathogenbase as pth
import map_transfer_matrix as mtm
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
           'LTAC': 'LTAC'}

constantFileNameMap = {'NURSINGHOME': 'nursinghome',
                       'HOSPITAL': 'hospital',
                       'COMMUNITY': 'community',
                       'LTAC': 'ltac'}

careTiers = ['HOME', 'NURSING', 'LTAC', 'HOSP', 'ICU']


defaultNotesPath = '/home/welling/workspace/pyRHEA/src/sim/notes.pkl'
defaultDescrDirs = ['/home/welling/workspace/pyRHEA/models/OrangeCounty/facilityfactsCurrent',
                    '/home/welling/workspace/pyRHEA/models/OrangeCounty/communitiesCurrent']
defaultConstDir = '/home/welling/workspace/pyRHEA/src/sim/facilityImplementations'


def fullCRVFromMeanLOS(fitParms):
    mean = fitParms[0]
    sigma = fitParms[1]
    mu = math.log(mean) - (0.5 * sigma * sigma)
    return lognorm(sigma, scale=math.exp(mu), loc=0.0)


def fullCRVFromLOSModel(losModel):
    """
    Returns a function of the signature: pscore = losFun([x0, x1, x2, ...]) for use
    in plotting analytic PDFs.  The curve is shifted right by 0.5 because the bar
    chart it must overlay centers the bar for integer N at x=N, but that bar really
    represents the integral of the PDF from (N-1) to N and so should be centered at
    x = (N - 0.5).
    """
    if losModel['pdf'] == 'lognorm(mu=$0,sigma=$1)':
        mu, sigma = losModel['parms']
        return lognorm(sigma, scale=math.exp(mu), loc=0.0)
    elif losModel['pdf'] == '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)':
        k, mu, sigma, lmda = losModel['parms']
        return lognormplusexp(s=sigma, mu=mu, k=k, lmda=lmda)
    elif losModel['pdf'] == 'expon(lambda=$0)':
        lmda = losModel['parms'][0]
        return expon(scale=1.0/lmda)
    else:
        raise RuntimeError('Unknown LOS model %s' % losModel['pdf'])


class LOSPlotter(object):
    def __init__(self, descr, constants):
        self.fullCRVs = {}
        if descr['category'] == 'HOSPITAL':
            if 'meanLOSICU' in descr:
                self.fullCRVs['ICU'] = fullCRVFromMeanLOS([descr['meanLOSICU'],
                                                           constants['icuLOSLogNormSigma']])
            if 'losModel' in descr:
                self.fullCRVs['HOSP'] = fullCRVFromLOSModel(descr['losModel'])
        elif descr['category'] == 'LTAC':
            if 'losModel' in descr:
                self.fullCRVs['LTAC'] = fullCRVFromLOSModel(descr['losModel'])
        elif descr['category'] == 'NURSINGHOME':
            if 'losModel' in descr:
                self.fullCRVs['NURSING'] = fullCRVFromLOSModel(descr['losModel'])
            else:
                self.fullCRVs['NURSING'] = fullCRVFromLOSModel(constants['nhLOSModel'])
        elif descr['category'] == 'COMMUNITY':
            self.fullCRVs['HOME'] = fullCRVFromLOSModel(constants['communityLOSModel'])
        else:
            raise RuntimeError("facility %s has unknown category %s - cannot plot LOS"
                               % (descr['abbrev'], descr['category']))

    def plot(self, tier, axes, nBins, rMin, rMax, scale, pattern='r-'):
        """The curve is shifted right by 0.5 because the bar chart it must overlay centers
        the bar for integer N at x=N, but that bar really represents the integral of the PDF
        from (N-1) to N and so should be centered at x = (N - 0.5)."""
        if tier in self.fullCRVs:
            curveX = np.linspace(rMin, rMax, nBins)
            boundedScale = scale / (self.fullCRVs[tier].cdf(rMax)
                                    - self.fullCRVs[tier].cdf(rMin))
            curveY = self.fullCRVs[tier].pdf(curveX - 0.5) * boundedScale
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


def collectBarSamples(histoVal):
    bins = []
    counts = []
    pairList = histoVal.histogram().items()
    pairList.sort()
    quantum = histoVal.d['quantum']
    barWidth = 1.6 * quantum
    for b, c in pairList:
        bins.append(b - 0.5*(quantum + barWidth))
        counts.append(c)
    return bins, counts, barWidth


def overallLOSFig(catNames, allOfCategoryDict):
    figs1, axes1 = plt.subplots(nrows=len(allOfCategoryDict), ncols=1)
    for offset, cat in enumerate(catNames):
        constants = loadFacilityTypeConstants(constantFileNameMap[cat])
        losPlotter = LOSPlotter({'abbrev': 'all', 'category': cat}, constants)
        bigHisto = HistoVal([])
        for tier in careTiers:
            try:
                bigHisto += allOfCategoryDict[cat][tier + '_LOS']
            except:
                pass
        bins, counts, barWidth = collectBarSamples(bigHisto)
        rects = axes1[offset].bar(bins, counts, width=barWidth, color='b')  # @UnusedVariable
        axes1[offset].set_ylabel('Counts')
        axes1[offset].set_xlabel('Days')
        axes1[offset].set_title('LOS for ' + cat)
        if bins:
            losPlotter.plot(nameMap[cat], axes1[offset], 300, 0, max(bins), sum(counts))
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
                tK = tier + '_LOS'
                if tK in notesDict[k]:
                    bins, counts, barWidth = collectBarSamples(notesDict[k][tK])
                    rects = axes1b[idx].bar(bins, counts, width=barWidth,  # @UnusedVariable
                                            color='b')
                    losPlotter.plot(tier, axes1b[idx], 300, 0.0, max(bins),
                                    sum(counts))
                else:
                    rects = axes1b[idx].bar([], [], color='b')  # @UnusedVariable
                axes1b[idx].set_ylabel('Counts')
                axes1b[idx].set_xlabel('Days')
                axes1b[idx].set_title('LOS for ' + tier)
            break
    figs1b.tight_layout()
    figs1b.canvas.set_window_title("%s LOS Histograms By Category" % abbrev)


def bedBounceFig(allOfCategoryDict):
    catNames = allOfCategoryDict.keys()[:]
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
    catNames = allOfCategoryDict.keys()[:]
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
    clrMap = {'death': 'black',
              'HOME': 'green',
              'NURSING': 'red',
              'HOSP': 'blue',
              'ICU': 'cyan',
              'LTAC': 'yellow'}
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
        clrs = []
        for ct, lbl in zip(counts, labels):
            if ct != 0:
                ct2.append(ct)
                lbl2.append(lbl)
                clrs.append(clrMap[lbl])
        axes4[offset].pie(ct2, labels=lbl2, autopct='%1.1f%%', startangle=90, colors=clrs)
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
                assert len(l) == len(dayList), (('field %s is the wrong length in special data %s'
                                                 '(%d vs. %d)')
                                                % (k, patchName, len(l), len(dayList)))
                axes5[offset].plot(dayList, l, label=k)
            axes5[offset].set_xlabel('Days')
            axes5[offset].set_ylabel('Occupancy')
            axes5[offset].legend()
        except Exception, e:
            print e
        axes5[offset].set_title(patchName)
    figs5.tight_layout()
    figs5.canvas.set_window_title("Time History of Occupancy")


# def pathogenTimeFig(specialDict):
#     figs6, axes6 = plt.subplots(nrows=1, ncols=len(specialDict))
#     if len(specialDict) == 1:
#         axes6 = [axes6]
#     for offset, (patchName, data) in enumerate(specialDict.items()):
#         try:
#             pthDataList = data['pathogen']
#             assert isinstance(pthDataList, types.ListType), \
#                 'Special data %s is not a list' % patchName
#             fields = {}
#             for d in pthDataList:
#                 for k, v in d.items():
#                     if k not in fields:
#                         fields[k] = []
#                     fields[k].append(v)
#             assert 'day' in fields, 'Date field is missing for special data %s' % patchName
#             dayList = fields['day']
#             del fields['day']
#             keys = fields.keys()
#             keys.sort()
#             for k in keys:
#                 facClassNm, enumStr = k.split('_')
#                 l = fields[k]
#                 lbl = '%s %s' % (facClassNm, pth.Status.names[int(enumStr)])
#                 assert len(l) == len(dayList), (('field %s is the wrong length in special data %s'
#                                                  '(%d vs. %d)')
#                                                 % (k, patchName, len(l), len(dayList)))
#                 if not all([ll == 0 for ll in l]) and facClassNm == 'Hospital':
#                     axes6[offset].plot(dayList, l, label=lbl)
#             axes6[offset].set_xlabel('Days')
#             axes6[offset].set_ylabel('Pathogen Prevalence')
#             axes6[offset].legend()
#         except Exception, e:
#             print e
#         axes6[offset].set_title(patchName)
#     figs6.tight_layout()
#     figs6.canvas.set_window_title("Time History of Occupancy")


def pathogenTimeFig(specialDict):
    patchList = specialDict.keys()[:]
    patchList.sort()
    catList = []
    for patchName, data in specialDict.items():
        pthDataList = data['pathogen']
        assert isinstance(pthDataList, types.ListType), 'Special data %s is not a list' % patchName
        for d in pthDataList:
            for k in d.keys():
                if k != 'day':
                    cat = k.split('_')[0]
                    if cat not in catList:
                        catList.append(cat)
    catList.sort()
    figs6, axes6 = plt.subplots(nrows=len(catList), ncols=len(patchList))
    axes6.reshape((len(catList), len(patchList)))
    if len(catList) == 1:
        axes6 = axes6[np.newaxis, :]
    if len(patchList) == 1:
        axes6 = axes6[:, np.newaxis]
    for colOff, patchName in enumerate(patchList):
        try:
            pthDataList = specialDict[patchName]['pathogen']
            assert isinstance(pthDataList, types.ListType), \
                'Special data %s is not a list' % patchName
            fields = {}
            for d in pthDataList:
                for k, v in d.items():
                    if k not in fields:
                        fields[k] = []
                    fields[k].append(v)
            assert 'day' in fields, 'Date field is missing for special data %s' % patchName
            dayList = fields['day']
            del fields['day']

            for rowOff, cat in enumerate(catList):
                for pthLvl in xrange(len(pth.Status.names)):
                    key = '%s_%d' % (cat, pthLvl)
                    if key in fields:
                        lbl = '%s' % pth.Status.names[pthLvl]
                        l = fields[key]
                        assert len(l) == len(dayList), (('field %s is the wrong length in special'
                                                         ' data %s (%d vs. %d)')
                                                        % (k, patchName, len(l), len(dayList)))
                        if not all([ll == 0 for ll in l]):
                            axes6[rowOff, colOff].plot(dayList, l, label=lbl)
                axes6[rowOff, colOff].set_xlabel('Days')
                axes6[rowOff, colOff].set_ylabel('Pathogen Prevalence')
                axes6[rowOff, colOff].legend()
                axes6[rowOff, colOff].set_title(cat)
        except Exception, e:
            print e
    figs6.tight_layout()
    figs6.canvas.set_window_title("Time History of Infection Status")


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


def buildTransferMap(catNames, categoryDict):
    transDict = {}
    for catNm in catNames:
        for fac, d in categoryDict[catNm].items():
            if fac not in transDict:
                transDict[fac] = {}
            for k, v in d.items():
                if k.endswith('_transfer'):
                    dest = k[:-9].lower()
                    if dest in transDict[fac]:
                        transDict[fac][dest] += v
                    else:
                        transDict[fac][dest] = v
    return transDict


def writeTransferMapAsCSV(transferMap, fname):
    recs = []
    for src, d in transferMap.items():
        rec = {'': src.upper()}
        for dst, count in d.items():
            rec['To_' + dst.upper()] = count
        recs.append(rec)
    allKeys = set(recs[0].keys())
    for rec in recs[1:]:
        allKeys.update(rec.keys())
    kL = list(allKeys)
    kL.sort()
    for rec in recs:
        for k in kL:
            if k not in rec:
                rec[k] = 0
    with open('sim_transfer_matrix.csv', 'w') as f:
        csv_tools.writeCSV(f, kL, recs)


def writeTransferMapAsDot(transferDict, fname):
    facDict = mtm.parseFacilityData('/home/welling/workspace/pyRHEA/models/OrangeCounty/'
                                    'facilityfactsCurrent')

    mtm.initializeMapCoordinates(facDict.values())

    inclusionSet = [('NURSINGHOME', 'HOSPITAL'),
                    ('NURSINGHOME', 'LTAC'),
                    ('LTAC', 'NURSINGHOME'),
                    ('HOSPITAL', 'NURSINGHOME'),
                    ('HOSPITAL', 'HOSPITAL'),
                    ('NURSINGHOME', 'NURSINGHOME'),
                    ('LTAC', 'LTAC'),
                    ('LTAC', 'HOSPITAL'),
                    ('HOSPITAL', 'LTAC')
                    ]

    mtm.writeDotGraph('sim_graph.dot', 'Simulated patient transfers',
                      facDict, transferDict, inclusionSet)


def main():
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

    writeTransferMapAsDot(buildTransferMap(catNames, categoryDict),
                          'sim_transfer_matrix.csv')

    countBirthsDeaths(catNames, allOfCategoryDict)

    overallLOSFig(catNames, allOfCategoryDict)

    singleLOSFig('SJUD', notesDict)
    singleLOSFig('WAEC', notesDict)
    singleLOSFig('CM69', notesDict)
    singleLOSFig('COLL', notesDict)
    bedBounceFig(allOfCategoryDict)
    patientFlowFig(allOfCategoryDict)
    patientFateFig(catNames, allOfCategoryDict)
    occupancyTimeFig(specialDict)
    pathogenTimeFig(specialDict)

    plt.show()

if __name__ == "__main__":
    main()
