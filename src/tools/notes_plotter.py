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

import sys
import os.path
import logging
import logging.config
from optparse import OptionParser

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))

import re
import yaml
import math
import pickle
import ujson
import types
import glob
from imp import load_source
from collections import defaultdict

if 'line_profiler' not in dir():
    def profile(func):
        """
        Dummy substitute for __builtin__.profile line profiler
        """
        return func

import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.notes.noteholder as noteholder
import pyrheautils
from facilitybase import CareTier as CareTierEnum
from facilitybase import PatientOverallHealth as OverallHealthEnum
import schemautils
from phacsl.utils.notes.statval import HistoVal
from stats import lognormplusexp, doubleweibull, doubleexpon
from stats import fullCRVFromPDFModel, fullLogNormCRVFromMean
import pathogenbase as pth
import map_transfer_matrix as mtm
import tools_util as tu
from pyrhea import getLoggerConfig

import numpy as np
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

from scipy.stats import lognorm, expon

logger = None

SCHEMA_DIR = '../schemata'
INPUT_SCHEMA = 'rhea_input_schema.yaml'

DEFAULT_NOTES_FNAME = 'notes.pkl'

CARE_TIERS = CareTierEnum.names.values()[:]

_CACHED_PARSE_DICT = {}
_CACHED_PARSE_ID = None   # the id of the info from which the cached dict was generated

FAC_TYPE_TO_CATEGORY_MAP = {'NursingHome': ['SNF', 'NURSINGHOME'],
                            'LTAC': ['LTACH', 'LTAC'],
                            'Community': ['COMMUNITY'],
                            'VentSNF': ['VSNF'],
                            'Hospital': ['HOSPITAL']}

def checkInputFileSchema(fname, schemaFname):
    try:
        with open(fname, 'rU') as f:
            inputJSON = yaml.safe_load(f)
        if os.name != "nt":
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
        else:
            return inputJSON
    except Exception, e:
        sys.exit('Error checking input against its schema: %s' % e)


def fullCRVFromMeanLOS(fitParms):
    return fullLogNormCRVFromMean(fitParms[0], fitParms[1])


def fullCRVFromLOSModel(losModel):
    return fullCRVFromPDFModel(losModel)


def buildEnumNameStr(name, ind):
    return '%s_%d' % (name, ind)

def splitEnumNameStr(enumName):
    offset = enumName.rfind('_')
    patchName = enumName[:offset]
    ind = int(enumName[offset+1:])
    return patchName, ind

def buildEnumEnumNameStr(name, ind1, ind2):
    return '%s_%d_%d' % (name, ind1, ind2)

def splitEnumEnumNameStr(enumName):
    frontStr, ind2 = splitEnumNameStr(enumName)
    patchName, ind1 = splitEnumNameStr(frontStr)
    return patchName, ind1, ind2


@profile
def mergeNotesFiles(notesPathList, globFlag=False):
    """
    Given a list of path strings pointing to notes files, produce a dict containing all the
    time series info and all 'unusual' info in those notes.
    """
    specialDict = {}
    if globFlag:
        newNotesPathList = []
        for notesFName in notesPathList:
            print '%s yields %s' % (notesFName, glob.glob(notesFName))
            newNotesPathList += glob.glob(notesFName)
        newNotesPathList.sort()
        notesPathList = newNotesPathList
    if notesPathList:
        for ind, notesFName in enumerate(notesPathList):
            notesDict = importNotes(notesFName)
            for nm, dct in notesDict.items():
                if '_' not in nm or nm.startswith('Patch'):
                    specialDict[buildEnumNameStr(nm, ind)] = dct
    return specialDict


@profile
def getTimeSeriesList(locKey, specialDict, specialDictKey):
    """
    key specifies the target entity, for example a loc abbrev or a facility category. These are
        typically in uppercase, but that is not enforced.
    specialDict is of the form output by mergeNotesFiles().
    specialDictKey is one of the second-level keys of specialDict, for example 'localpathogen'.
    Returns a list of tuples, of three possible forms.

    The first form appears if the time series contains
    entries for multiple different status levels for multiple tiers, for example
    populations sorted by CareTier at different PthLevel values.
    That form is:

          [(dayArray, valArrayD), (dayArray, valArrayD), ...]

    where dayArray is a numpy array of dates and valArrayD has the form:

          {(intLvl, intLvl): valArray, (intLvl, intLvl): valArray, ...}

    and the ints in the tuple represent the two indices, for example (tier, pthStatus).

    The second form appears if the time series contains
    entries for multiple different status levels, for example populations at different PthLevel values.
    That form is:

          [(dayArray, valArrayD), (dayArray, valArrayD), ...]

    where dayArray is a numpy array of dates and valArrayD has the form:

          {intLvl: valArray, intLvl: valArray, ...}

    with pathLvl being an integer status index (like a PthLevel) and fracArray being a numpy array counts at
    that level.

    The third form appears if the time series data is not classified by level, for example a simple population
    count.  That form is:

          [(dayArray, valArray), (dayArray, valArray), ...]

    """

    global _CACHED_PARSE_DICT
    global _CACHED_PARSE_ID

    if _CACHED_PARSE_ID != id(specialDict) or specialDictKey not in _CACHED_PARSE_DICT:
        _CACHED_PARSE_ID = id(specialDict)
        thisCD = _CACHED_PARSE_DICT[specialDictKey] = defaultdict(list)
        for tpl in timeSeriesListGenerator(specialDict, specialDictKey):
            notesInd, patchName, (loc, ind1, ind2), dayV, dataV = tpl
            thisCD[loc].append((notesInd, ind1, ind2, dayV, dataV))

    tplL = _CACHED_PARSE_DICT[specialDictKey][locKey]
    resultD = {}
    for tpl in tplL:
        notesInd, ind1, ind2, dayV, dataV = tpl
        if ind1 is None:
            assert notesInd not in resultD, 'Generator returned duplicate notes index?'
            resultD[notesInd] = (dayV, dataV)
        else:
            innerKey = ind1 if ind2 is None else (ind1, ind2)
            if notesInd not in resultD:
                resultD[notesInd] = (dayV, {innerKey: dataV})
            else:
                oldDayV, innerD = resultD[notesInd]
                assert dayV is oldDayV, 'Generator provided inconsistent day vectors?'
                innerD[innerKey] = dataV
    return resultD.values()


def timeSeriesListGenerator(specialDict, specialDictKey):
    """
    This generator returns all the time series data in specialDict having the associated
    top-level key (for example 'occupancy' or 'localtierpathogen'.  For each time series,
    the data returned is:

    (seqId, patchName, (locKey, ind1, ind2), dayVec, valVec)

    where seqId is the sequence number in order of the series of notes.pkl files being
    parsed, patchName is the name of the current patch within that notes.pkl file,
    dayVec and valVec are numpy arrays containing date and value vectors,
    and the tuple has one of the following forms:

    category, None, None         for aggregate info like 'occupancy'

    category, lvl, None          for aggregate info like 'pathogen'

    category, tier, None         for aggregate info like 'localtierarrivals'

    abbrev, lvl, None            for aggregate info like 'localpathogen'

    abbrev, tier, lvl            for aggregate info like 'localtierpathogen'

    where category is for example 'NursingHome', lvl is an integer PthLevel, and tier is
    an integer CareTier.
    """
    patchList = []
    indList = []
    for enumPatchName in specialDict:
        patchName, ind = splitEnumNameStr(enumPatchName)
        if patchName not in patchList:
            patchList.append(patchName)
        if ind not in indList:
            indList.append(ind)
    patchList.sort()
    indList.sort()
    rsltL = []
    for ind in indList:
        for patchName in patchList:
            enumPatchName = buildEnumNameStr(patchName, ind)
            if enumPatchName not in specialDict or specialDictKey not in specialDict[enumPatchName]:
                continue  # No guarantee all patches have all inds, or the field of interest
            dataList = specialDict[enumPatchName][specialDictKey]
            assert isinstance(dataList, types.ListType), \
                'Special data %s is not a list' % enumPatchName
            fields = defaultdict(list)
            for d in dataList:
                for k, v in d.items():
                    fields[k].append(v)
            assert 'day' in fields, 'Date field is missing for special data %s' % patchName
            dayV = np.array(fields['day'])
            del fields['day']

            for key, fld in fields.items():
                fld = np.array(fld)
                try:
                    loc, subI1, subI2 = splitEnumEnumNameStr(key)
                except ValueError:
                    try:
                        loc, subI1 = splitEnumNameStr(key)
                        subI2 = None
                    except ValueError:
                        loc = key
                        subI1 = subI2 = None
                yield (ind, patchName, (loc, subI1, subI2), dayV, fld)


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
        elif facToImplDict[descr['category']] == 'VSNF':
            if 'losModel' in descr:
                self.fullCRVs['VSNF'] = fullCRVFromLOSModel(descr['losModel'])
            else:
                self.fullCRVs['VSNF'] = fullCRVFromLOSModel(constants['nhLOSModel'])
        elif facToImplDict[descr['category']] == 'COMMUNITY':
            if 'communityLOSModel' in constants:
                self.fullCRVs['HOME'] = fullCRVFromLOSModel(constants['communityLOSModel'])
            elif 'losModelMap' in constants and 'initialUnhealthyFrac' in constants:
                # This is a hack based on a specific 'classifier' in the community implementation
                assert ('HEALTHY_base' in constants['losModelMap']
                        and constants['losModelMap']['HEALTHY_base']['pdf'] == 'expon(lambda=$0)'
                        and 'UNHEALTHY_base' in constants['losModelMap']
                        and constants['losModelMap']['HEALTHY_base']['pdf'] == 'expon(lambda=$0)'
                        ), 'Unexpected losModelMap format for community'
                lmda1 = constants['losModelMap']['HEALTHY_base']['parms'][0]
                lmda2 = constants['losModelMap']['UNHEALTHY_base']['parms'][0]
                k = constants['initialUnhealthyFrac']
                print 'creating doubleexpon(%s, %s, %s)' % (k, lmda1, lmda2)
                self.fullCRVs['HOME'] = doubleexpon(k, lmda1, lmda2)
            else:
                self.fullCRVs['HOME'] = None
        else:
            raise RuntimeError("facility %s has unknown category %s - cannot plot LOS"
                               % (descr['abbrev'], descr['category']))

    def plot(self, tier, axes, nBins, rMin, rMax, scale, pattern='r-'):
        """The curve is shifted right by 0.5 because the bar chart it must overlay centers
        the bar for integer N at x=N, but that bar really represents the integral of the PDF
        from (N-1) to N and so should be centered at x = (N - 0.5)."""
        if tier in self.fullCRVs:
            if self.fullCRVs[tier] is not None:  # meaning we have no LOS model
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
    sys.path.append(implementationDir)

    fullPath = pyrheautils.pathTranslate(os.path.join(implementationDir,
                                                      category.lower() + '.py'))
    newMod = load_source(category.lower(), fullPath)
    assert hasattr(newMod, 'generateFull'), ('%s does not contain a facility implementation' %
                                             fullPath)
    assert hasattr(newMod, 'category') and newMod.category.lower() == category.lower(), \
        ('%s does not contain an implementation of %s' % (fullPath, category))
    assert hasattr(newMod, '_constants_values') and hasattr(newMod, '_constants_schema'), \
        ('%s does not point to its associated constants' % fullPath)

    constPath = pyrheautils.pathTranslate(newMod._constants_values)
    schemaPath = pyrheautils.pathTranslate(newMod._constants_schema)
    jsonDict = pyrheautils.importConstants(constPath, schemaPath)
    sys.path.pop()  # drop implementaionDir
    return yaml_tools._simplify(jsonDict)


def importNotes(fname):
    try:
        print "loading pkl"
        with open(fname, 'r') as f:
            stuff = pickle.load(f)
    except (KeyError, pickle.UnpicklingError):
        print "loading json"
        with open(fname, 'r') as f:
            stuff = ujson.load(f)
    return stuff


def collectBarSamples(histoVal):
    bins = []
    counts = []
    pairList = histoVal.histogram().items()
    pairList.sort()
    quantum = histoVal.d['quantum']
    barWidth = 1.6 * quantum
    for b, c in pairList:
        bins.append(b - 0.5*quantum)
        counts.append(c)
    return bins, counts, barWidth


def overallLOSFig(catNames, allOfCategoryDict, facToImplDict, constDir):
    nPlots = 0
    for cat in catNames:
        for tier in CARE_TIERS:
            if cat in allOfCategoryDict and (tier + '_LOS') in allOfCategoryDict[cat]:
                nPlots += 1
    figs1, axes1 = plt.subplots(nrows=nPlots, ncols=1)
    if nPlots == 1:
        axes1 = [axes1]
    offset = 0
    for cat in catNames:
        if cat not in allOfCategoryDict:
            continue
        constants = loadFacilityTypeConstants(facToImplDict[cat].lower(), constDir)
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


def singleLOSFig(abbrev, notesDict, facilityDirs, facToImplDict, constDir):
    figs1b, axes1b = plt.subplots(nrows=len(CARE_TIERS), ncols=1)
    losDescr = loadFacilityDescription(abbrev, facilityDirs)
    constants = loadFacilityTypeConstants(facToImplDict[losDescr['category']].lower(), constDir)
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


def patientFateFigClassic(catNames, allOfCategoryDict, allFacInfo, catToImplDict):
#     fig, axes = plt.subplot()
    fig = plt.figure(figsize=(len(catNames), 2))
    axes = fig.gca()
    axes.set_xlim((-0.5, len(catNames)-0.5))
    axes.set_ylim((-0.5, 1.5))
    axes.set_aspect('equal', 'datalim')
    axes.set_xticks(range(len(catNames)))
    axes.set_yticks([0.0, 1.0])
    axes.set_xticklabels(catNames)
    axes.set_yticklabels(['real', 'sim'])
    clrMap = {'death': 'black',
              'HOME': 'green',
              'other': 'green',
              'NURSING': 'red',
              'VSNF': 'pink',
              'SKILNRS': 'gray',
              'VENT': 'purple',
              'NURSING+SKILNRS+VENT': 'gray',
              'HOSP': 'blue',
              'HOSP+ICU': 'blue',
              'ICU': 'cyan',
              'LTAC': 'yellow'}
    implTierMap = {'HOSPITAL': ['HOSP', 'ICU'],
                   'LTAC': ['LTAC'],
                   'NURSINGHOME': ['NURSING'],
                   'COMMUNITY': ['HOME'],
                   'VSNF': ['NURSING', 'SKILNRS', 'VENT']}
    wrapOrderMap = {'ICU': 0,
                    'HOSP': 1,
                    'HOSP+ICU': 1,
                    'LTAC': 2,
                    'NURSING': 3,
                    'VSNF': 4,
                    'SKILNRS': 5,
                    'NURSING+SKILNRS+VENT': 4,
                    'VENT': 6,
                    'HOME': 7,
                    'other': 7,
                    'death': 8}
    for offset, cat in enumerate(catNames):
        keys = ['death'] + ['%s_found' % tier for tier in CARE_TIERS]
        labels = ['death'] + CARE_TIERS
        counts = []
        for k in keys:
            if k in allOfCategoryDict[cat]:
                counts.append(allOfCategoryDict[cat][k])
            else:
                counts.append(0)
        sortMe = [(wrapOrderMap[lbl], (ct, lbl, clrMap[lbl])) for ct, lbl in zip(counts, labels)
                  if ct != 0]
        sortMe.sort(reverse=True)
        ct2, lbl2, clrs = [list(tpl) for tpl in zip(*[tpl for idx, tpl in sortMe])]
        row = 1
        axes.pie(ct2, labels=lbl2, autopct='%1.1f%%', startangle=90, colors=clrs,
                radius=0.25, center=(offset, row), frame=True)
        row = 0
        labels = []
        values = []
        clrs = []
        pairDict = {}
        if cat in allFacInfo:
            for lbl, val in allFacInfo[cat].items():
                if lbl in catToImplDict:
                    toImpl = catToImplDict[lbl]
                    toLbl = '+'.join(implTierMap[toImpl])
                else:
                    toLbl = lbl
#                 print '%s %s' % (toLbl, pairDict)
                if toLbl in pairDict:
                    pairDict[toLbl] += val
                else:
                    pairDict[toLbl] = val
                clrs.append(clrMap[toLbl])
            if pairDict:
                sortMe = [(wrapOrderMap[lbl], (ct, lbl, clrMap[lbl])) for lbl, ct in pairDict.items()
                          if ct != 0]
                sortMe.sort(reverse=True)
                values, labels, clrs = [list(tpl) for tpl in zip(*[tpl for idx, tpl in sortMe])]
            else:
                labels = []
                values = []
                clrs = []
        axes.pie(values, labels=labels, autopct='%1.1f%%', startangle=90, colors=clrs,
                radius=0.25, center=(offset, row), frame=True)
 
    fig.tight_layout()
    fig.canvas.set_window_title("Patient Fates (Classic mode)")
#     fig.tight_layout()
#     fig.canvas.set_window_title("Patient Fates (Classic mode)")


FLOW_CLR_MAP = {'death': 'black',
              'birth': 'black',
              'COMMUNITY': 'green',
              'HOME': 'green',
              'other': 'green',
              'NURSINGHOME': 'red',
              'NURSING': 'red',
              'VSNF': 'pink',
              'SKILNRS': 'gray',
              'VENT': 'purple',
              'NURSING+SKILNRS+VENT': 'gray',
              'HOSPITAL': 'blue',
              'HOSP': 'blue',
              'HOSP+ICU': 'blue',
              'ICU': 'cyan',
              'LTACH': 'yellow',
              'LTAC': 'yellow'}


def pieHelper(pairDict):
    implTierMap = {'HOSPITAL': ['HOSP', 'ICU'],
                   'LTAC': ['LTAC'],
                   'NURSINGHOME': ['NURSING'],
                   'COMMUNITY': ['HOME'],
                   'VSNF': ['NURSING', 'SKILNRS', 'VENT']}
    wrapOrderMap = {'ICU': 0,
                    'HOSPITAL': 1,
                    'HOSP': 1,
                    'HOSP+ICU': 1,
                    'LTACH': 2,
                    'LTAC': 2,
                    'NURSINGHOME': 3,
                    'NURSING': 3,
                    'VSNF': 4,
                    'SKILNRS': 5,
                    'NURSING+SKILNRS+VENT': 4,
                    'VENT': 6,
                    'COMMUNITY': 7,
                    'HOME': 7,
                    'other': 7,
                    'death': 8,
                    'birth': 8}
    if pairDict:
        sortMe = [(wrapOrderMap[lbl], (ct, lbl, FLOW_CLR_MAP[lbl])) for lbl, ct in pairDict.items()
                  if ct != 0]
        sortMe.sort(reverse=True)
        #print 'sortMe: ',sortMe
        #print 'zip: ', [tpl for idx, tpl in sortMe]
        values, labels, clrs = [list(tpl) for tpl in zip(*[tpl for idx, tpl in sortMe])]
    else:
        labels = []
        values = []
        clrs = []
    return labels, values, clrs


def drawFlowBarPlot(axes, countMtx, simCountMtx, categoryL, nRealCols):
    colWid = 0.25
    xV = np.arange(nRealCols)
    baseV0 = np.zeros(nRealCols)
    baseV1 = np.zeros(nRealCols)
    artistL = [Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)]
    labelL = ['measured | sim']
    for idx, ctg in enumerate(categoryL):
        yV0 = countMtx[:nRealCols, idx]
        artistL.append(axes.bar(xV - 0.5 * colWid, yV0, colWid, bottom=baseV0,
                                color=FLOW_CLR_MAP[ctg], edgecolor='black'))
        labelL.append(ctg)
        baseV0 += yV0
        yV1 = simCountMtx[:nRealCols, idx]
        axes.bar(xV + 0.5 * colWid, yV1, colWid, bottom=baseV1, 
                 color=FLOW_CLR_MAP[ctg], edgecolor='black')
        baseV1 += yV1
    axes.set_xticks(xV)
    axes.set_xticklabels(categoryL[:nRealCols])
    axes.legend(artistL, labelL)

    


def patientFateFig(categoryDict, facDict, numDays, pieMode=True):
    categoryS = set([rec['category'] for rec in facDict.values()] + ['COMMUNITY'])
    categoryL = list(categoryS)
    nRealCols = len(categoryL)
    categoryL.sort()
    categoryL.append('death')
    # indices for countMtx are [destinationCat, source count for patients in destinationCat]
    catIdxD = {catNm:idx for idx, catNm in enumerate(categoryL)}
    # indices for countMtx are [srcCat, destination count for patients in srcCat]
    countMtx = np.zeros((len(categoryL), len(categoryL)))
    simCountMtx = np.zeros((len(categoryL), len(categoryL)))

    for rec in facDict.values():
        srcCat = rec['category']
        if srcCat == 'COMMUNITY':
            pass
        else:
            dischPerYear = rec['totalDischarges']['value']
            if 'deathRate' in rec:
                deathsPerYear = rec['deathRate']['value'] * dischPerYear
                dischPerYear -= deathsPerYear
            else:
                deathsPerYear = 0.0
            tTO = 0.0
            for recRow in rec['totalTransfersOut']:
                tTO += recRow['count']['value']
                countMtx[catIdxD[srcCat], catIdxD[recRow['category']]] += recRow['count']['value']
            countMtx[catIdxD[srcCat], catIdxD['COMMUNITY']] += dischPerYear - tTO
            # estimate of number of arrivals from community per year
            fmComCt = rec['totalDischarges']['value'] - rec['totalTransfersIn']['value']
            countMtx[catIdxD['COMMUNITY'], catIdxD[srcCat]] += fmComCt
            countMtx[catIdxD[srcCat], catIdxD['death']] += deathsPerYear

    for category in categoryL:
        if category not in categoryDict:  #special case, for ex. 'death'
            continue
        for subKey in categoryDict[category]:
            srcLoc = subKey.upper()
            assert srcLoc in facDict, 'Unknown srcLoc %s' % srcLoc
            srcCat = facDict[srcLoc]['category']
            assert srcCat == category, 'Table structure is not as expected'
            srcIdx = catIdxD[srcCat]
            for subSubKey in categoryDict[category][subKey]:
                if subSubKey.startswith('fromTo:'):
                    dstLoc, tier = subSubKey.split(':')[1:]
                    assert dstLoc in facDict, 'Unknown dstLoc %s' % dstLoc
                    dstCat = facDict[dstLoc]['category']
                    dstIdx = catIdxD[dstCat]
                    simCountMtx[srcIdx, dstIdx] += categoryDict[category][subKey][subSubKey]
                elif subSubKey == 'death':
                    ct = categoryDict[category][subKey][subSubKey]
                    simCountMtx[srcIdx, catIdxD['death']] += ct

    countScale = 365.0 / numDays

    print 'Patient fate matrices follow'
    print 'approx reality:'
    print countMtx[0:nRealCols, :]

    print 'simulation:'
    print (countScale * simCountMtx)[0:nRealCols, :]

    if pieMode:
        fig = plt.figure(figsize= (nRealCols, 2))
        axes = fig.gca()
        axes.set_xlim((-0.5, nRealCols - 0.5))
        axes.set_ylim((-0.5, 1.5))
        axes.set_aspect('equal', 'datalim')
        axes.set_xticks(range(nRealCols))
        axes.set_yticks([0.0, 1.0])
        axes.set_xticklabels(categoryL)
        axes.set_yticklabels(['approx real', 'sim'])
        for srcIdx, srcCat in enumerate(categoryL[:nRealCols]):

            row = 0
            pairDict = {dstCat : countMtx[srcIdx, dstIdx]
                        for dstIdx, dstCat in enumerate(categoryL)}
            labels, values, clrs = pieHelper(pairDict)
            axes.pie(values, labels=labels, colors=clrs, autopct='%1.1f%%', startangle=90,
                     radius=0.25, frame=True, center=(srcIdx, row))
            row = 1
            pairDict = {dstCat : simCountMtx[srcIdx, dstIdx]
                        for dstIdx, dstCat in enumerate(categoryL)}
            labels, values, clrs = pieHelper(pairDict)
            scaleFac = countScale * (np.sum(simCountMtx[srcIdx, :])/np.sum(countMtx[srcIdx, :]))
            axes.pie(values, labels=labels, colors=clrs, autopct='%1.1f%%', startangle=90,
                     radius=(scaleFac * 0.25), frame=True, center=(srcIdx, row))
    else:
        fig = plt.figure()
        axes = fig.gca()

        drawFlowBarPlot(axes, countMtx, countScale * simCountMtx, categoryL, nRealCols)

    fig.tight_layout()
    fig.canvas.set_window_title("Patient Fates")


def patientSourceFig(categoryDict, facDict, numDays, pieMode=True):
    categoryS = set([rec['category'] for rec in facDict.values()] + ['COMMUNITY'])
    categoryL = list(categoryS)
    nRealCols = len(categoryL)
    categoryL.sort()
    categoryL.append('birth')
    catIdxD = {catNm:idx for idx, catNm in enumerate(categoryL)}
    # indices for countMtx are [destinationCat, source count for patients in destinationCat]
    countMtx = np.zeros((len(categoryL), len(categoryL)))
    simCountMtx = np.zeros((len(categoryL), len(categoryL)))
    for rec in facDict.values():
        srcCat = rec['category']
        if srcCat == 'COMMUNITY':
            pass
        else:
            dischPerYear = rec['totalDischarges']['value']
            if 'deathRate' in rec:
                deathsPerYear = rec['deathRate']['value'] * dischPerYear
                dischPerYear -= deathsPerYear
            else:
                deathsPerYear = 0.0
            tTO = 0.0
            for recRow in rec['totalTransfersOut']:
                tTO += recRow['count']['value']
                countMtx[catIdxD[recRow['category']], catIdxD[srcCat]] += recRow['count']['value']
            countMtx[catIdxD['COMMUNITY'], catIdxD[srcCat]] += dischPerYear - tTO
            # estimate of number of arrivals from community per year
            fmComCt = rec['totalDischarges']['value'] - rec['totalTransfersIn']['value']
            countMtx[catIdxD[srcCat], catIdxD['COMMUNITY']] += fmComCt
            countMtx[catIdxD['COMMUNITY'], catIdxD['birth']] += deathsPerYear

    for category in categoryL:
        if category not in categoryDict:  # special case, for ex 'birth'
            continue
        for subKey in categoryDict[category]:
            srcLoc = subKey.upper()
            assert srcLoc in facDict, 'Unknown srcLoc %s' % srcLoc
            srcCat = facDict[srcLoc]['category']
            assert srcCat == category, 'Table structure is not as expected'
            srcIdx = catIdxD[srcCat]
            for subSubKey in categoryDict[category][subKey]:
                if subSubKey.startswith('fromTo:'):
                    dstLoc, tier = subSubKey.split(':')[1:]
                    assert dstLoc in facDict, 'Unknown dstLoc %s' % dstLoc
                    dstCat = facDict[dstLoc]['category']
#                     if category != 'COMMUNITY':
#                         print '%s %s (%s %s)-> %s (%s) %s' %(srcLoc, tier, srcCat, category, dstLoc, dstCat, 
#                                                              categoryDict[category][subKey][subSubKey])
                    dstIdx = catIdxD[dstCat]
                    simCountMtx[dstIdx, srcIdx] += categoryDict[category][subKey][subSubKey]
                elif subSubKey == 'birth':
                    ct = categoryDict[category][subKey][subSubKey]
                    assert srcCat == 'COMMUNITY', ('%s has %d births but is not a community?'
                                                   % (subKey, ct))
                    simCountMtx[srcIdx, catIdxD['birth']] += ct

    countScale = 365.0 / numDays

    print 'Patient source matrices follow'
    print 'approx reality:'
    print countMtx[:, 0:nRealCols]

    print 'simulation:'
    print (countScale * simCountMtx)[:, 0:nRealCols]

    if pieMode:
        fig = plt.figure(figsize=(nRealCols, 2))
        axes = fig.gca()
        axes.set_xlim((-0.5, nRealCols-0.5))
        axes.set_ylim((-0.5, 1.5))
        axes.set_aspect('equal', 'datalim')
        axes.set_xticks(range(nRealCols))
        axes.set_yticks([0.0, 1.0])
        axes.set_xticklabels(categoryL[:nRealCols])
        axes.set_yticklabels(['approx real', 'sim'])
    
        for dstIdx, dstCat in enumerate(categoryL[:nRealCols]):
            row = 0
            pairDict = {srcCat : countMtx[dstIdx, srcIdx]
                        for srcIdx, srcCat in enumerate(categoryL)}        
            labels, values, clrs = pieHelper(pairDict)
            axes.pie(values, labels=labels, colors=clrs, autopct='%1.1f%%', startangle=90,
                     radius=0.25, frame=True, center=(dstIdx, row))
            row = 1
            pairDict = {srcCat : simCountMtx[dstIdx, srcIdx]
                        for srcIdx, srcCat in enumerate(categoryL)}
            labels, values, clrs = pieHelper(pairDict)
            scaleFac = countScale * (np.sum(simCountMtx[dstIdx, :])/np.sum(countMtx[dstIdx, :]))
            axes.pie(values, labels=labels, colors=clrs, autopct='%1.1f%%', startangle=90,
                     radius=(scaleFac * 0.25), frame=True, center=(dstIdx, row))
    else:
        fig = plt.figure()
        axes = fig.gca()

        drawFlowBarPlot(axes, countMtx, countScale * simCountMtx, categoryL, nRealCols)        

    fig.tight_layout()
    fig.canvas.set_window_title("Patient Sources")


def getDayRange(inputDict):
    minDay = inputDict['burnInDays']
    maxDay = minDay + inputDict['runDurationDays'] - 1
    return minDay, maxDay

def occupancyTimeFig(specialDict, meanPopByCat=None):
    figs5, axes5 = plt.subplots(nrows=1, ncols=len(specialDict))
    if len(specialDict) == 1:
        axes5 = [axes5]
    for offset, (patchName, data) in enumerate(specialDict.items()):
        try:
            occDataList = data['occupancy']
            assert isinstance(occDataList, types.ListType), \
                'Special data %s is not a list' % patchName
            fields = defaultdict(list)
            for d in occDataList:
                for k, v in d.items():
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
                    baseLine, = axes5[offset].plot(dayList, l, label=k)
                    if meanPopByCat is not None and k in FAC_TYPE_TO_CATEGORY_MAP:
                        for cat in FAC_TYPE_TO_CATEGORY_MAP[k]:
                            if cat in meanPopByCat:
                                meanPop = meanPopByCat[cat]
                                axes5[offset].plot(dayList, [meanPop] * len(dayList),
                                                   color=baseLine.get_color(),
                                                   linestyle='--')

            axes5[offset].set_xlabel('Days')
            axes5[offset].set_ylabel('Occupancy')
            axes5[offset].legend()
        except Exception, e:
            print e
        axes5[offset].set_title(patchName)
    figs5.tight_layout()
    figs5.canvas.set_window_title("Time History of Occupancy")


def pathogenTimeFig(specialDict):
    try:
        patchList = specialDict.keys()[:]
        patchList.sort()
        catList = []
        for patchName, data in specialDict.items():
            pthDataList = data['pathogen']
            assert isinstance(pthDataList, types.ListType), ('Special data %s is not a list'
                                                             % patchName)
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
                dayVec = np.array(dayList)
                del fields['day']

                for rowOff, cat in enumerate(catList):
                    curvesThisCat = {}
                    for pthLvl in xrange(len(pth.PthStatus.names)):
                        key = '%s_%d' % (cat, pthLvl)
                        if key in fields:
                            l = fields[key]
                            assert len(l) == len(dayList), (('field %s is the wrong length in special'
                                                             ' data %s (%d vs. %d)')
                                                            % (key, patchName, len(l), len(dayList)))
                            curvesThisCat[pthLvl] = np.array(l)
                    totsThisCat = sum(curvesThisCat.values())
                    for pthLvl, lVec in curvesThisCat.items():
                        lbl = '%s' % pth.PthStatus.names[pthLvl]
                        if np.count_nonzero(lVec):
                            with np.errstate(divide='ignore', invalid='ignore'):
                                scaleV = np.true_divide(np.asfarray(lVec), np.asfarray(totsThisCat))
                                scaleV[scaleV == np.inf] = 0.0
                                scaleV = np.nan_to_num(scaleV)
                                axes6[rowOff, colOff].plot(dayVec, scaleV, label=lbl)
                    axes6[rowOff, colOff].set_xlabel('Days')
                    axes6[rowOff, colOff].set_ylabel('Pathogen Prevalence')
                    axes6[rowOff, colOff].legend()
                    axes6[rowOff, colOff].set_title(cat)
            except Exception, e:
                print e
        figs6.tight_layout()
        figs6.canvas.set_window_title("Time History of Infection Status")
    except KeyError as e:
        if e.message == 'pathogen':
            print('pathogen data not available; skipping time history of infection status')


def colonizationTimeFig(specialDict, facDict):
    try:
        patchSet = set()
        catSet = set()
        tierSet = set()
        for tSL in timeSeriesListGenerator(specialDict, 'localtiernewcolonized'):
            ind, patchName, tpl, dayV, valV = tSL
            fac, tier, shouldBeNone = tpl
            assert ind == 0, 'There should be only one notes file?'
            assert shouldBeNone is None, 'localtiernewcolonized data has unexpected format'
            if valV.any():
                patchSet.add(patchName)
                catSet.add(facDict[fac]['category'])
                tierSet.add(tier)
        patchL = list(patchSet)
        patchL.sort()
        catL = list(catSet)
        catL.sort()
        tierL = list(tierSet)
        tierL.sort()
        fig, axes = plt.subplots(nrows=len(tierL), ncols=len(patchL))
        if len(tierL) == 1 and len(patchL) == 1:
            axes = np.asarray([[axes]])
        else:
            axes.reshape((len(tierL), len(patchL)))
            if len(tierL) == 1:
                axes = axes[np.newaxis, :]
            if len(patchL) == 1:
                axes = axes[:, np.newaxis]
        dayVD = {}
        valVD = {}
        for tSL in timeSeriesListGenerator(specialDict, 'localtiernewcolonized'):
            ind, patchName, tpl, dayV, valV = tSL
            fac, tier, shouldBeNone = tpl
            if valV.any():
                vKey = (patchName, tier)
#                 print 'before: %s' % valV[-10:]
#                 valV[1:] -= valV[0:-1]
#                 print 'after: %s' % valV[-10:]
                if vKey in dayVD:
                    assert np.all(dayV == dayVD[vKey]), ('Day vectors do not match for patch %s %s'
                                                         % (patchName, CareTierEnum.names[tier]))
                    valVD[vKey] += valV
                else:
                    dayVD[vKey] = dayV.copy()
                    valVD[vKey] = valV.copy()
        for colOff, patchName in enumerate(patchL):
            for rowOff, tier in enumerate(tierL):
                vKey = (patchName, tier)
                axes[rowOff, colOff].plot(dayVD[vKey], valVD[vKey])
                axes[rowOff, colOff].set_xlabel('Days')
                axes[rowOff, colOff].set_ylabel('New Colonizations')
                axes[rowOff, colOff].legend()
                axes[rowOff, colOff].set_title(CareTierEnum.names[tier])
        fig.tight_layout()
        fig.canvas.set_window_title("Time History of New Colonizations")
    except Exception as e:
        print('new colonizations fig failed with exception %s' % e)
        raise


def overallHealthTimeFig(specialDict):
    try:
        patchList = specialDict.keys()[:]
        patchList.sort()
        catList = []
        for patchName, data in specialDict.items():
            ohDataList = data['occupancyByOH']
            assert isinstance(ohDataList, types.ListType), ('Special data %s is not a list'
                                                            % patchName)
            for d in ohDataList:
                for k in d.keys():
                    if k != 'day':
                        cat = k.split('_')[0]
                        if cat not in catList:
                            catList.append(cat)
        catList.sort()
        figsOH, axesOH = plt.subplots(nrows=len(catList), ncols=len(patchList))
        axesOH.reshape((len(catList), len(patchList)))
        if len(catList) == 1:
            axesOH = axesOH[np.newaxis, :]
        if len(patchList) == 1:
            axesOH = axesOH[:, np.newaxis]
        for colOff, patchName in enumerate(patchList):
            try:
                ohDataList = specialDict[patchName]['occupancyByOH']
                assert isinstance(ohDataList, types.ListType), \
                    'Special data %s is not a list' % patchName
                fields = {}
                for d in ohDataList:
                    for k, v in d.items():
                        if k not in fields:
                            fields[k] = []
                        fields[k].append(v)
                assert 'day' in fields, 'Date field is missing for special data %s' % patchName
                dayList = fields['day']
                dayVec = np.array(dayList)
                del fields['day']

                for rowOff, cat in enumerate(catList):
                    curvesThisCat = {}
                    for ohLvl in xrange(len(OverallHealthEnum.names)):
                        key = '%s_%d' % (cat, ohLvl)
                        if key in fields:
                            l = fields[key]
                            assert len(l) == len(dayList), (('field %s is the wrong length in special'
                                                             ' data %s (%d vs. %d)')
                                                            % (key, patchName, len(l), len(dayList)))
                            curvesThisCat[ohLvl] = np.array(l)
                    totsThisCat = sum(curvesThisCat.values())
                    for ohLvl, lVec in curvesThisCat.items():
                        lbl = '%s' % OverallHealthEnum.names[ohLvl]
                        if np.count_nonzero(lVec):
                            with np.errstate(divide='ignore', invalid='ignore'):
                                scaleV = np.true_divide(np.asfarray(lVec), np.asfarray(totsThisCat))
                                scaleV[scaleV == np.inf] = 0.0
                                scaleV = np.nan_to_num(scaleV)
                                axesOH[rowOff, colOff].plot(dayVec, lVec, label=lbl,
                                                            color='C%d' % ohLvl)
                    axesOH[rowOff, colOff].set_xlabel('Days')
                    axesOH[rowOff, colOff].set_ylabel('Overall Health')
                    axesOH[rowOff, colOff].legend()
                    axesOH[rowOff, colOff].set_title(cat)
            except Exception, e:
                print e
        figsOH.tight_layout()
        figsOH.canvas.set_window_title("Time History of Patient Population Overall Health")
    except KeyError as e:
        if e.message == 'occupancyByOH':
            print('occupancyByOH data not available; skipping time history of'
                  ' overall health')


def miscTimeFig(specialDict):
    patchList = specialDict.keys()[:]
    patchList.sort()
    try:
        catList = []
        for patchName, data in specialDict.items():
            for topic in ['localtiernthawed', 'bedHoldStats']:
                if topic in data:
                    topicDataList = data[topic]
                    assert isinstance(topicDataList, types.ListType), ('Special data %s is not a list'
                                                                     % patchName)
                    for d in topicDataList:
                        for k in d.keys():
                            if k != 'day':
                                cat = k.split('_')[0]
                                if cat not in catList:
                                    catList.append(cat)
        catList.sort()
        for colOff, patchName in enumerate(patchList):
            allFields = {}
            for topic in ['localtiernthawed', 'bedHoldStats']:
                if topic in data:
                    try:
                        topicDataList = specialDict[patchName][topic]
                        assert isinstance(topicDataList, types.ListType), \
                            'Special data %s %s is not a list' % (topic, patchName)
                        fields = {}
                        for d in topicDataList:
                            for key, val in d.items():
                                if key in catList:
                                    fieldNm = topic
                                    catNm = key
                                elif key == 'day':
                                    fieldNm = key
                                    catNm = None
                                else:
                                    words = key.split('_', 1)
                                    if words[0] in catList:
                                        fieldNm = words[1]
                                        catNm = words[0]
                                    else:
                                        raise RuntimeError('Unparsable key %s' % key)
                                if catNm not in fields:
                                    fields[catNm] = {}
                                if fieldNm not in fields[catNm]:
                                    fields[catNm][fieldNm] = []
                                fields[catNm][fieldNm].append(val)
                            assert 'day' in fields[None], ('Date field is missing for special data %s %s'
                                                           % (patchName))
                        dayVec = np.asarray(fields[None]['day'])
                        del fields[None]['day']
                        del fields[None]
                        for catNm in fields.keys():
                            for key, lst in fields[catNm].items():
                                valV = np.asarray(lst)
                                if np.count_nonzero(valV):
                                    fields[catNm][key] = valV
                                else:
                                    del fields[catNm][key]
                            fields[catNm]['day'] = dayVec
                        allFields.update(fields)
                    except Exception, e:
                        print e
                        raise

        # Trim out categories for which we have only a 'day' vector
        catList = [catNm for catNm in catList if len(allFields[catNm]) > 1]
        figsMisc, axesMisc = plt.subplots(nrows=len(catList), ncols=len(patchList))
        if hasattr(axesMisc, 'reshape'):
            axesMisc.reshape((len(catList), len(patchList)))
            if len(catList) == 1:
                axesMisc = axesMisc[np.newaxis, :]
            if len(patchList) == 1:
                axesMisc = axesMisc[:, np.newaxis]
        else:
            axesMisc = np.asarray([[axesMisc]])
        for colOff, patchName in enumerate(patchList):
            for rowOff, catNm in enumerate(catList):
                if catNm in allFields:
                    fldD = allFields[catNm]
                    dayVec = fldD['day']
                    toPlotL = fldD.keys()
                    toPlotL.remove('day')
                    if 'frail' in toPlotL and 'frail_held' in toPlotL:
                        fldD['commited_frail'] = fldD['frail'] + fldD['frail_held']
                        del fldD['frail']
                        del fldD['frail_held']
                    if 'non_frail' in toPlotL and 'non_frail_held' in toPlotL:
                        fldD['commited_non_frail'] = fldD['non_frail'] + fldD['non_frail_held']
                        del fldD['non_frail']
                        del fldD['non_frail_held']
                    labelD = {'commited_frail' : 'commited frail',
                             'commited_non_frail' : 'committed non-frail',
                             'localtiernthawed' : 'num thawed'
                             }
                    for fldNm, yVec in fldD.items():
                        if fldNm == 'day':
                            continue
                        assert dayVec.shape == yVec.shape, ('day vec mismatch for %s' % fldNm)
                        label = labelD[fldNm] if fldNm in labelD else fldNm
                        if np.count_nonzero(yVec):
                            axesMisc[rowOff, colOff].plot(dayVec, yVec, '-', label=label)
                axesMisc[rowOff, colOff].set_xlabel('Days')
                axesMisc[rowOff, colOff].set_title(catNm)
                axesMisc[rowOff, colOff].legend()
            figsMisc.tight_layout()
            figsMisc.canvas.set_window_title("Misc Time Histories")
    except KeyError as e:
        if e.message == 'localtiernthawed':
            print('localtiernthawed data not available; skipping time history of thaws')
        else:
            raise


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

def readFacFiles(facilityDirs):
    return mtm.parseFacilityData(facilityDirs)

def scanAllFacilities(facilityDirs, facDict=None):
    transOutByCat = defaultdict(dict)
    meanPopByCat = defaultdict(lambda: 0.0)
    if not facDict:
        facDict = readFacFiles(facilityDirs)
    for fac in facDict.values():
        cat = fac['category']
        if 'meanPop' in fac:
            meanPopByCat[cat] += fac['meanPop']['value']
        else:
            print 'Excluding %s from summaries; no meanPop' % fac['abbrev']
            continue
        if 'totalDischarges' in fac:
            totDisch = fac['totalDischarges']['value']
        else:
            #if fac['category'] != 'COMMUNITY':
            #    print fac['category']
            #    print fac
            assert fac['category'] == 'COMMUNITY', '%s should have totDisch and does not' % fac['abbrev']
            totDisch = None
        if 'totalTransfersOut' in fac:
            knownDisch = 0
            for dct in fac['totalTransfersOut']:
                toCat = dct['category']
                toN = dct['count']['value']
                if toCat not in transOutByCat[cat]:
                    transOutByCat[cat][toCat] = 0
                transOutByCat[cat][toCat] += toN
                knownDisch += toN
            if totDisch is not None:
                delta = totDisch - knownDisch
                if 'other' in transOutByCat[cat]:
                    transOutByCat[cat]['other'] += delta
                else:
                    transOutByCat[cat]['other'] = delta

    return transOutByCat, meanPopByCat

def main():
    """
    main
    """

    global logger
    logging.config.dictConfig(getLoggerConfig())
    logger = logging.getLogger(__name__)

    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl] run_descr.yaml
    """)
    parser.add_option('-n', '--notes', action='store', type='string',
                      help="Notes filename (overrides any name in the run description)")
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    parser.destroy()

    runDesc = args[0]

    inputDict = tu.readModelInputs(args[0])
    pyrheautils.prepPathTranslations(inputDict)
    facDict = tu.getFacDict(inputDict)

    implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])

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
            if '_' not in nm:
                specialDict[nm] = dct
                continue
            if nm.startswith('Patch') and not nm.endswith('Registry'):
                specialDict[nm] = dct
                continue
            if nm.startswith('Patch') and nm.endswith('Registry'):
                continue
            category, abbrev = tuple(nm.split('_', 1))
            abbrev = abbrev.lower()
            if category not in categoryDict:
                categoryDict[category] = {}
            categoryDict[category][abbrev] = dct
        except:
            specialDict[nm] = dct

    noteHolderGroup = noteholder.NoteHolderGroup()
    allOfCategoryDict = {}
    for category in categoryDict.keys():
        allOfCategoryDict[category] = noteHolderGroup.createNoteHolder()
        for abbrev, nhDict in categoryDict[category].items():
            allOfCategoryDict[category].addNote({k: v for k, v in nhDict.items() if k != 'name'})

    catNames = allOfCategoryDict.keys()[:]

    facDirList = [pyrheautils.pathTranslate(pth) for pth in inputDict['facilityDirs']]
    allOfCategoryFacilityInfo, meanPopByCategory = scanAllFacilities(facDirList)

    if 'facilitySelectors' in inputDict:
        facImplRules = [(re.compile(rule['category']), rule['implementation'])
                       for rule in inputDict['facilitySelectors']]
    else:
        facImplRules = [(re.compile(cat), cat)
                        for cat in catNames]  # an identity map

    catToImplDict = {cat: findFacImplCategory(implDir, facImplRules, cat)
                     for cat in catNames}

#     writeTransferMapAsDot(buildTransferMap(catNames, categoryDict),
#                           'sim_transfer_matrix.csv',
#                           facDirList, catToImplDict)

    minDay, maxDay = getDayRange(inputDict)
    print 'day range: ', minDay, maxDay

    countBirthsDeaths(catNames, allOfCategoryDict)

    try:
        overallLOSFig(catNames, allOfCategoryDict, catToImplDict, implDir)
    except Exception, e:
        logger.error('Exception in overallLOSFig: %s' % e)
 
    if "2013" in runDesc:
        #singleLOSFig('SJUD', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
        #singleLOSFig('WAEC', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
        #singleLOSFig('CM69', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
        #singleLOSFig('COLL', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
        pass
 
    if 'trackedFacilities' in inputDict:
        for abbrev in inputDict['trackedFacilities']:
            try:
                singleLOSFig(abbrev, notesDict, facDirList, catToImplDict, implDir)
            except Exception, e:
                logger.error('Exception in singleLOSFig for %s: %s' % (abbrev, e))
 
    try:
        bedBounceFig(allOfCategoryDict)
    except Exception, e:
        logger.error('Exception in bedBounceFig: %s' % e)
    try:
        patientFlowFig(allOfCategoryDict)
    except Exception, e:
        logger.error('Exception in patientFlowFig: %s' % e)
    try:
        patientSourceFig(categoryDict, facDict, maxDay + 1 - minDay, pieMode=False)
    except Exception, e:
        logger.error('Exception in patientSourceFig: %s' % e)
    try:
        patientFateFig(categoryDict, facDict, maxDay + 1 - minDay, pieMode=False)
        # patientFateFigClassic(catNames, allOfCategoryDict, allOfCategoryFacilityInfo, catToImplDict)
    except Exception, e:
        logger.error('Exception in patientFateFig: %s' % e)
    try:
        occupancyTimeFig(specialDict, meanPopByCat=meanPopByCategory)
    except Exception, e:
        logger.error('Exception in occupancyTimeFig: %s' % e)
    try:
        pathogenTimeFig(specialDict)
    except Exception, e:
        logger.error('Exception in pathogenTimeFig: %s' % e)
    try:
        colonizationTimeFig(specialDict, facDict)
    except Exception, e:
        logger.error('Exception in colonizationTimeFig: %s' % e)
    try:
        overallHealthTimeFig(specialDict)
    except Exception, e:
        logger.error('Exception in overallHealthTimeFig: %s' % e)
    try:
        miscTimeFig(specialDict)
    except Exception, e:
        logger.error('Exception in miscTimeFig: %s' % e)

    plt.show()

if __name__ == "__main__":
    main()
