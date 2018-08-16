#! /usr/bin/env python

import sys
import os
import pandas as pd
import numpy as np
import yaml

import matplotlib
matplotlib.use('Agg')  # hack for no display with pyplot
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.legend_handler import HandlerTuple
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from pathogenbase import PthStatus
from plotting_utilities import pltCurvesWithBounds

def loadMpzFiles(mpzDir, inFacS):
    selSumL = []
    idx = 0
    for fname in os.listdir(mpzDir):
        if not fname.endswith('.mpz'):
            continue
        print 'loading %s' % fname
        df = pd.read_msgpack(os.path.join(mpzDir, fname))
        df['TOTAL'] = df[PthStatus.names.values()].sum(axis=1)
        selDF = df[df['abbrev'].isin(inFacS)]
        selGrps = selDF.groupby(['tier', 'day'])
        selSums = selGrps.sum().add_suffix('_sum').reset_index()
        selSums['run'] = idx
        selSums['prevalence'] = (selSums['COLONIZED_sum'].astype(float)
                                 /selSums['TOTAL_sum'].astype(float))
        selSumL.append(selSums)
        idx += 1
    selDF = pd.concat(selSumL)
    return selDF

def readSelFacs(yamlPath):
    with open(yamlPath, 'rU') as f:
        json = yaml.load(f)
    return json['locationsImplementingScenario']['locAbbrevList']

baseDir = os.path.join(os.path.dirname(__file__),
                       'cre_bundle_201808')
scenarioL = ['baseline_CRE_Prevalence_13_Mile',
             'CRE_Prevalence_13_Mile',
             'Min_SNA_CRE_13_Mile']
inFacS = frozenset(readSelFacs(os.path.join(baseDir,
                                            scenarioL[0],
                                            'xdro_facs.yaml')))
fig, axes = plt.subplots(1, 1)

allArtists = []
allLabels = []
tierNm = 'LTAC'
for scen in scenarioL:
    dirPth = os.path.join(baseDir, scen)
    df = loadMpzFiles(dirPth, inFacS)

    artistL, labelL = pltCurvesWithBounds(df, axes,
                                          'prevalence', 'day', 'tier',
                                          [tierNm],
                                          ['%s %s' % (tierNm, scen)])
    allArtists += artistL
    allLabels += labelL

axes.legend(allArtists, allLabels, handler_map={tuple: HandlerTuple()})
axes.set_title('Prevalence for ' + tierNm)
plt.savefig('realitycheck.png')

