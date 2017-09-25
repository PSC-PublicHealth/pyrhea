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
import logging
import logging.config
from optparse import OptionParser

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))

import re
import yaml
import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.notes.noteholder as noteholder
import pyrheautils
from facilitybase import CareTier as CareTierEnum
import schemautils
from phacsl.utils.notes.statval import HistoVal
from stats import lognormplusexp
import pathogenbase as pth
import map_transfer_matrix as mtm
import math
import pickle
import ujson
import types
from imp import load_source
from collections import defaultdict
from pyrhea import getLoggerConfig

import numpy as np

from scipy.stats import lognorm, expon

logger = None

SCHEMA_DIR = '../schemata'
INPUT_SCHEMA = 'rhea_input_schema.yaml'

DEFAULT_NOTES_FNAME = 'notes.pkl'

CARE_TIERS = CareTierEnum.names.values()[:]

FAC_TYPE_TO_CATEGORY_MAP = {'NursingHome': 'SNF',
                            'LTAC': 'LTACH',
                            'Community': 'COMMUNITY',
                            'VentSNF': 'VSNF',
                            'Hospital': 'HOSPITAL'}

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
        bins.append(b - 0.5*(quantum + barWidth))
        counts.append(c)
    return bins, counts, barWidth

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
        assert 'meanPop' in fac, '%s has no meanPop' % fac['abbrev']
        meanPopByCat[cat] += fac['meanPop']['value']
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
        
    schemautils.setSchemaBasePath(SCHEMA_DIR)
    inputDict = checkInputFileSchema(args[0],
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    modelDir = inputDict['modelDir']
    pyrheautils.PATH_STRING_MAP['MODELDIR'] = modelDir
    implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])
    pyrheautils.PATH_STRING_MAP['IMPLDIR'] = implDir
    constDir = os.path.join(modelDir, 'constants')
    
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
    
    facDirList = [pyrheautils.pathTranslate(pth) for pth in inputDict['facilityDirs']]
    if "ChicagoLand" in runDesc:
        allOfCategoryFacilityInfo, meanPopByCategory = scanAllFacilities(facDirList)

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
                          facDirList, catToImplDict)
    
    countBirthsDeaths(catNames, allOfCategoryDict)


    if "2013" in runDesc:
        #singleLOSFig('SJUD', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
        #singleLOSFig('WAEC', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
        #singleLOSFig('CM69', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
        #singleLOSFig('COLL', notesDict, inputDict['facilityDirs'], catToImplDict, implDir)
        pass
    


if __name__ == "__main__":
    main()
