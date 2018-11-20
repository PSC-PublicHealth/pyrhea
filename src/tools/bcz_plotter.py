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
import numpy as np
if __name__ == "__main__":
    import matplotlib.pyplot as plt

from scipy.stats import lognorm, expon
import pandas as pd
import bcolz

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
from bcz_monitor import generateBCZDirName
from notes_plotter import findFacImplCategory, checkInputFileSchema, loadFacilityDescription
from notes_plotter import loadFacilityTypeConstants

logger = None

SCHEMA_DIR = '../schemata'
INPUT_SCHEMA = 'rhea_input_schema.yaml'

DEFAULT_NOTES_FNAME = 'notes.pkl'

CARE_TIERS = CareTierEnum.names.values()[:]


def importBCZ(basePath):
    patchId = 0
    pth = generateBCZDirName(basePath, patchId)
    df = bcolz.ctable(rootdir=pth, mode='r').todataframe()
    if df is None:
        sys.exit('Failed to read %s' % pth)
    df['patch'] = patchId
    while True:
        patchId += 1
        pth = generateBCZDirName(basePath, patchId)
        if os.path.isdir(pth):
            subDF = bcolz.ctable(rootdir=pth, mode='r').todataframe()
            subDF['patch'] = patchId
            if tu.pandasVersionIsAtLeast23():
                df = pd.concat((df, subDF), sort=True)
            else:
                df = pd.concat((df, subDF))
        else:
            break
    return df

def addCategory(fac, facDict):
    return facDict[fac]['category']


def occupancyTimeFig(allDF, facDict, meanPopByCat=None):
    allDF['category'] = allDF['fac'].apply(addCategory, args=(facDict,))
    df = allDF.groupby(['day', 'category', 'patch']).sum().reset_index()
    nPatch = df['patch'].nunique()
    nCat = df['category'].nunique()

    fig, axes = plt.subplots(nrows=1, ncols=nPatch)
    if nPatch == 1:
        axes = [axes]
    patchOff = 0
    for patch, patchDF in df.groupby(['patch']):
        for cat, patchCatDF in patchDF.groupby(['category']):
            if np.any(patchCatDF['TOTAL'].values):
                baseline, = axes[patchOff].plot(patchCatDF['day'].values,
                                                patchCatDF['TOTAL'].values,
                                                '-', label=cat)
                meanPop = meanPopByCat[cat]
                yline = np.full(patchCatDF['day'].values.shape, meanPop)
                axes[patchOff].plot(patchCatDF['day'].values, yline,
                                    color=baseline.get_color(),
                                    linestyle='--')
            axes[patchOff].legend()
            axes[patchOff].set_xlabel('Days')
            axes[patchOff].set_ylabel('count')
            axes[patchOff].set_title('patch %s' % patchOff)
        patchOff += 1
        fig.suptitle('Population by Category')
    fig.tight_layout()
    fig.canvas.set_window_title("Time History of Occupancy")


def pathogenTimeFig(allDF):
    df = allDF.groupby(['day', 'tier', 'patch']).sum().reset_index()
    nPatch = df['patch'].nunique()
    nTier = df['tier'].nunique()

    fig, axes = plt.subplots(nrows=nTier, ncols=nPatch)
    axes.reshape((nTier, nPatch))
    if nTier == 1:
        axes = axes[np.newaxis, :]
    if nPatch == 1:
        axes = axes[:, np.newaxis]
    patchOff = 0
    for patch, patchDF in df.groupby(['patch']):
        tierOff = 0
        for tier, patchTierDF in patchDF.groupby(['tier']):
            totV = patchTierDF['TOTAL'].values.astype(np.float)
            for pthStatus in pth.PthStatus.names.values():
                if np.any(patchTierDF[pthStatus].values):
                    yV = patchTierDF[pthStatus].values.astype(np.float)/totV
                    axes[tierOff, patchOff].plot(patchTierDF['day'].values,
                                                 yV, '-', label=pthStatus)
            axes[tierOff, patchOff].legend()
            axes[tierOff, patchOff].set_xlabel('Days')
            axes[tierOff, patchOff].set_ylabel('count')
            axes[tierOff, patchOff].set_title(tier)
            tierOff += 1
        patchOff += 1
        fig.suptitle('Pathogen Status by Tier')


def colonizationTimeFig(allDF):
    df = allDF.groupby(['day', 'tier', 'patch']).sum().reset_index()
    nPatch = df['patch'].nunique()
    nTier = df['tier'].nunique()

    fig, axes = plt.subplots(nrows=nTier, ncols=nPatch)
    axes.reshape((nTier, nPatch))
    if nTier == 1:
        axes = axes[np.newaxis, :]
    if nPatch == 1:
        axes = axes[:, np.newaxis]
    patchOff = 0
    for patch, patchDF in df.groupby(['patch']):
        tierOff = 0
        for tier, patchTierDF in patchDF.groupby(['tier']):
            axes[tierOff, patchOff].plot(patchTierDF['day'].values,
                                         patchTierDF['localtiernewcolonized'].values,
                                         '-')
            axes[tierOff, patchOff].set_xlabel('Days')
            axes[tierOff, patchOff].set_ylabel('count')
            axes[tierOff, patchOff].set_title(tier)
            tierOff += 1
        patchOff += 1
        fig.suptitle('New Colonizations by Tier')


def cpTimeFig(allDF):
    df = allDF.groupby(['day', 'tier', 'patch']).sum().reset_index()
    nPatch = df['patch'].nunique()
    nTier = df['tier'].nunique()

    fig, axes = plt.subplots(nrows=nTier, ncols=nPatch)
    axes.reshape((nTier, nPatch))
    if nTier == 1:
        axes = axes[np.newaxis, :]
    if nPatch == 1:
        axes = axes[:, np.newaxis]
    patchOff = 0
    for patch, patchDF in df.groupby(['patch']):
        tierOff = 0
        for tier, patchTierDF in patchDF.groupby(['tier']):
            axes[tierOff, patchOff].plot(patchTierDF['day'].values,
                                         patchTierDF['localtierCP'].values,
                                         '-')
            axes[tierOff, patchOff].set_xlabel('Days')
            axes[tierOff, patchOff].set_ylabel('count')
            axes[tierOff, patchOff].set_title(tier)
            tierOff += 1
        patchOff += 1
        fig.suptitle('Contact Precautions by Tier')


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
    %prog --bczname basename run_descr.yaml
    """)
    parser.add_option('-b', '--bczname', action='store', type='string',
                      help="base name for bcolz directories", default=None)
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')
    if opts.bczname is None:
        parser.error('--bczname is a required argument')
    bczBaseName = opts.bczname

    parser.destroy()

    runDesc = args[0]

    inputDict = tu.readModelInputs(runDesc)
    pyrheautils.prepPathTranslations(inputDict)
    facDict = tu.getFacDict(inputDict)

    implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])

    allDF = importBCZ(bczBaseName)

    catNames = list(set([rec['category'] for rec in facDict.values()]))

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

    try:
        occupancyTimeFig(allDF, facDict, meanPopByCat=meanPopByCategory)
    except Exception, e:
        logger.error('Exception in occupancyTimeFig: %s' % e)
        raise
    try:
        pathogenTimeFig(allDF)
    except Exception, e:
        logger.error('Exception in pathogenTimeFig: %s' % e)
    if 'trackedValues' in inputDict:
        if 'newColonizations' in inputDict['trackedValues']:
            try:
                colonizationTimeFig(allDF)
            except Exception, e:
                logger.error('Exception in colonizationTimeFig: %s' % e)
        if 'creCounters' in inputDict['trackedValues']:
            try:
                cpTimeFig(allDF)
            except Exception, e:
                logger.error('Exception in cpTimeFig: %s' % e)

    plt.show()

if __name__ == "__main__":
    main()
