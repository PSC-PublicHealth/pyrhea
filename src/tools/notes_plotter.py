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

import sys
import os.path
from optparse import OptionParser

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))

import re
import yaml
import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.notes.noteholder as noteholder
from pyrheautils import importConstants
from facilitybase import CareTier as CareTierEnum
import schemautils
from phacsl.utils.notes.statval import HistoVal
from stats import lognormplusexp
import pathogenbase as pth
import map_transfer_matrix as mtm
import math
import pickle
import types

import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import lognorm, expon

SCHEMA_DIR = '../schemata'
INPUT_SCHEMA = 'rhea_input_schema.yaml'

DEFAULT_NOTES_FNAME = 'notes.pkl'

CARE_TIERS = CareTierEnum.names.values()[:]

def checkInputFileSchema(fname, schemaFname):
    try:
        with open(fname, 'rU') as f:
            inputJSON = yaml.safe_load(f)
        validator = schemautils.getValidator(os.path.join(os.path.dirname(__file__), schemaFname))
        nErrors = sum([1 for e in validator.iter_errors(inputJSON)])  # @UnusedVariable
        if nErrors:
            print 'Input file violates schema:'
            for e in validator.iter_errors(inputJSON):
                print ('Schema violation: %s: %s' %
                       (' '.join([str(word) for word in e.path]), e.message))
            sys.exit('Input file violates schema')
        else:
            return inputJSON
    except Exception, e:
        sys.exit('Error checking input against its schema: %s' % e)


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
    def __init__(self, descr, constants, facToImplDict):
        self.fullCRVs = {}
        if facToImplDict[descr['category']] == 'HOSPITAL':
            if 'meanLOSICU' in descr:
                self.fullCRVs['ICU'] = fullCRVFromMeanLOS([descr['meanLOSICU'],
                                                           constants['icuLOSLogNormSigma']])
            if 'losModel' in descr:
                self.fullCRVs['HOSP'] = fullCRVFromLOSModel(descr['losModel'])
        elif facToImplDict[descr['category']] == 'LTAC':
            if 'losModel' in descr:
                self.fullCRVs['LTAC'] = fullCRVFromLOSModel(descr['losModel'])
        elif facToImplDict[descr['category']] == 'NURSINGHOME':
            if 'losModel' in descr:
                self.fullCRVs['NURSING'] = fullCRVFromLOSModel(descr['losModel'])
            else:
                self.fullCRVs['NURSING'] = fullCRVFromLOSModel(constants['nhLOSModel'])
        elif facToImplDict[descr['category']] == 'COMMUNITY':
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


def loadFacilityDescription(abbrev, facilityDirs):
    for descDir in facilityDirs:
        fname = os.path.join(descDir, abbrev + '.yaml')
        if os.path.exists(fname):
            with open(fname, 'rU') as f:
                return yaml_tools._simplify(yaml.safe_load(f))
    raise RuntimeError('No description file found for %s' % abbrev)


def loadFacilityTypeConstants(category, implementationDir):
    jsonDict = importConstants(os.path.join(implementationDir,
                                            category.lower() + '_constants.yaml'),
                               os.path.join(SCHEMA_DIR,
                                            category.lower() + '_constants_schema.yaml'))
    return yaml_tools._simplify(jsonDict)


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


def overallLOSFig(catNames, allOfCategoryDict, facToImplDict, implDir):
    nPlots = 0
    for cat in catNames:
        for tier in CARE_TIERS:
            if cat in allOfCategoryDict and (tier + '_LOS') in allOfCategoryDict[cat]:
                nPlots += 1
    figs1, axes1 = plt.subplots(nrows=nPlots, ncols=1)
    offset = 0
    for cat in catNames:
        if cat not in allOfCategoryDict:
            continue
        constants = loadFacilityTypeConstants(facToImplDict[cat].lower(), implDir)
        losPlotter = LOSPlotter({'abbrev': 'all', 'category': cat}, constants,
                                facToImplDict)
        for tier in CARE_TIERS:
            if (tier + '_LOS') in allOfCategoryDict[cat]:
                bins, counts, barWidth = collectBarSamples(allOfCategoryDict[cat][tier + '_LOS'])
                axes1[offset].bar(bins, counts, width=barWidth, color='b')
                axes1[offset].set_ylabel('Counts')
                axes1[offset].set_xlabel('Days')
                axes1[offset].set_title('LOS for %s %s' % (cat, tier))
                if bins:
                    losPlotter.plot(tier, axes1[offset], 300, 0, max(bins), sum(counts))
                offset += 1
    figs1.tight_layout()
    figs1.canvas.set_window_title("LOS Histograms By Category")


def singleLOSFig(abbrev, notesDict, facilityDirs, facToImplDict, implDir):
    figs1b, axes1b = plt.subplots(nrows=len(CARE_TIERS), ncols=1)
    losDescr = loadFacilityDescription(abbrev, facilityDirs)
    constants = loadFacilityTypeConstants(facToImplDict[losDescr['category']].lower(), implDir)
    losPlotter = LOSPlotter(losDescr, constants, facToImplDict)
    for k in notesDict.keys():
        if k.endswith(abbrev):
            for idx, tier in enumerate(CARE_TIERS):
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
    figs2, axes2 = plt.subplots(nrows=len(CARE_TIERS), ncols=1)
    for offset, tier in enumerate(CARE_TIERS):
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
        axes2[offset].bar(bins, counts, color='r')  # @UnusedVariable
        axes2[offset].set_ylabel('Counts')
        axes2[offset].set_xlabel('Bounces')
        axes2[offset].set_title('Bounces for Bed Requests for ' + tier)
    # figs2.tight_layout()
    figs2.canvas.set_window_title("Bed Request Bounce Histograms")


def patientFlowFig(allOfCategoryDict):
    catNames = allOfCategoryDict.keys()[:]
    figs3, axes3 = plt.subplots(nrows=1, ncols=len(CARE_TIERS))
    for offset, tier in enumerate(CARE_TIERS):
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
        keys = ['death'] + ['%s_found' % tier for tier in CARE_TIERS]
        labels = ['death'] + CARE_TIERS
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
                if not k.endswith('_all'):
                    # The '_all' curves count every arrival at the facility and
                    # are no longer of interest
                    axes5[offset].plot(dayList, l, label=k)
            axes5[offset].set_xlabel('Days')
            axes5[offset].set_ylabel('Occupancy')
            axes5[offset].legend()
        except Exception, e:
            print e
        axes5[offset].set_title(patchName)
    figs5.tight_layout()
    figs5.canvas.set_window_title("Time History of Occupancy")


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


def writeTransferMapAsDot(transferDict, fname, facilityDirs, catToImplDict):
    facDict = mtm.parseFacilityData(facilityDirs)

    mtm.initializeMapCoordinates(facDict.values())

    inclusionSet = [
        ('NURSINGHOME', 'HOSPITAL'),
        ('NURSINGHOME', 'LTAC'),
        ('LTAC', 'NURSINGHOME'),
        ('HOSPITAL', 'NURSINGHOME'),
        ('HOSPITAL', 'HOSPITAL'),
        ('NURSINGHOME', 'NURSINGHOME'),
        ('LTAC', 'LTAC'),
        ('LTAC', 'HOSPITAL'),
        ('HOSPITAL', 'LTAC')
        ]
    
    catL = catToImplDict.keys()[:]
    translatedInclusionSet = []
    for src in catL:
        for dst in catL:
            if (catToImplDict[src], catToImplDict[dst]) in inclusionSet:
                translatedInclusionSet.append((src, dst))

    mtm.writeDotGraph('sim_graph.dot', 'Simulated patient transfers',
                      facDict, transferDict, translatedInclusionSet)


def findFacImplCategory(facImplDict,
                        facImplRules,
                        category):
    """
    The 'category' parameter comes from a facility description 'category' attribute.
    Map it to its matching facility implementation 'category' attribute using the
    supplied rules.
    """
    for facRegex, implStr in facImplRules:
        if facRegex.match(category):
            return implStr
    return None

def main():
    """
    main
    """
    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl] run_descr.yaml
    """)
    parser.add_option('-n', '--notes', action='store', type='string',
                      help="Notes filename (overrides any name in the run description)")
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    parser.destroy()

    inputDict = checkInputFileSchema(args[0],
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    implDir = inputDict['facilityImplementationDir']
    
    if opts.notes:
        notesFName = opts.notes
    elif 'notesFileName' in inputDict:
        notesFName = inputDict['notesFileName']
    else:
        notesFName = DEFAULT_NOTES_FNAME

    notesDict = importNotes(notesFName)
    categoryDict = {}
    specialDict = {}
    for nm, dct in notesDict.items():
        try:
            if '_' in nm and not nm.startswith('Patch'):
                category, abbrev = tuple(nm.split('_', 1))
                abbrev = abbrev.lower()
                if category not in categoryDict:
                    categoryDict[category] = {}
                categoryDict[category][abbrev] = dct
            else:
                specialDict[nm] = dct
        except:
            specialDict[nm] = dct

    noteHolderGroup = noteholder.NoteHolderGroup()
    allOfCategoryDict = {}
    for category in categoryDict.keys():
        allOfCategoryDict[category] = noteHolderGroup.createNoteHolder()
        for abbrev, nhDict in categoryDict[category].items():
            allOfCategoryDict[category].addNote({k: v for k, v in nhDict.items() if k != 'name'})

    catNames = allOfCategoryDict.keys()[:]

    if 'facilitySelectors' in inputDict:
        facImplRules = [(re.compile(rule['category']), rule['implementation'])
                       for rule in inputDict['facilitySelectors']]
    else:
        facImplRules = [(re.compile(cat), cat)
                        for cat in catNames]  # an identity map

    catToImplDict = {cat: findFacImplCategory(implDir, facImplRules, cat)
                     for cat in catNames}

    writeTransferMapAsDot(buildTransferMap(catNames, categoryDict),
                          'sim_transfer_matrix.csv',
                          inputDict['facilityDirs'], catToImplDict)

    countBirthsDeaths(catNames, allOfCategoryDict)

    overallLOSFig(catNames, allOfCategoryDict, catToImplDict, implDir)

#     singleLOSFig('SJUD', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
#     singleLOSFig('WAEC', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
#     singleLOSFig('CM69', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
#     singleLOSFig('COLL', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
    singleLOSFig('MANO_512_S', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)

    bedBounceFig(allOfCategoryDict)
    patientFlowFig(allOfCategoryDict)
    patientFateFig(catNames, allOfCategoryDict)
    occupancyTimeFig(specialDict)
    pathogenTimeFig(specialDict)

    plt.show()

if __name__ == "__main__":
    main()
