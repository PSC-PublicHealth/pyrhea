#! /usr/bin/env python

import pylab
import numpy as np
import pandas as pd
import scipy.stats
import scipy.sparse
from ggplot import *
import concurrent.futures
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os.path

import schemautils
import pyrheautils
from pyrhea import checkInputFileSchema
from map_transfer_matrix import parseFacilityData

SCHEMA_DIR = os.path.join(os.path.dirname(__file__), '../schemata')

# Constants for bootstrap sampling
BOOTSTRAP_SIZE = 1000
TOTAL_SAMPLES = 16000
NUM_WORKERS = 16

def softmax(X, theta = 1.0, axis = None):
    """
    Compute the softmax of each element along an axis of X.

    Parameters
    ----------
    X: ND-Array. Probably should be floats. 
    theta (optional): float parameter, used as a multiplier
        prior to exponentiation. Default = 1.0
    axis (optional): axis to compute values along. Default is the 
        first non-singleton axis.

    Returns an array the same size as X. The result will sum to 1
    along the specified axis.
    
    Implementation from: https://nolanbconaway.github.io/blog/2017/softmax-numpy
    """

    # make X at least 2d
    y = np.atleast_2d(X)

    # find axis
    if axis is None:
        axis = next(j[0] for j in enumerate(y.shape) if j[1] > 1)

    # multiply y against the theta parameter, 
    y = y * float(theta)

    # subtract the max for numerical stability
    y = y - np.expand_dims(np.max(y, axis = axis), axis)
    
    # exponentiate y
    y = np.exp(y)

    # take the sum along the specified axis
    ax_sum = np.expand_dims(np.sum(y, axis = axis), axis)

    # finally: divide elementwise
    p = y / ax_sum

    # flatten if X was 1D
    if len(X.shape) == 1: p = p.flatten()

    return p

def bootstrap_test(ma, mb, bootstrap_size=10000, num_samples=1000):
    def bootstrap_samples(ma, mb, bootstrap_size, num_samples):
        aa = ma.flatten()
        ba = mb.flatten()
        
        assert(len(aa)==len(ba))
        
        aprobs = aa.astype(float) / float(sum(aa))
        bprobs = ba.astype(float) / float(sum(ba))

        for sample in range(num_samples):
            aidxs = np.random.choice(len(aa), size=bootstrap_size, replace=True, p=aprobs)
            bidxs = np.random.choice(len(ba), size=bootstrap_size, replace=True, p=bprobs)
            permute_idxs = np.random.randint(0,2,bootstrap_size,dtype=bool)
            
            null_a_s = np.zeros(len(aa), dtype=float)
            null_b_s = np.zeros(len(aa), dtype=float)
            a_s = np.zeros(len(aa), dtype=float)
            b_s = np.zeros(len(aa), dtype=float)

            np.add.at(a_s, aidxs, 1)
            np.add.at(b_s, bidxs, 1)
            np.add.at(null_a_s, np.where(permute_idxs, aidxs, bidxs), 1)
            np.add.at(null_b_s, np.where(permute_idxs, bidxs, aidxs), 1)
            
            yield tuple(x.reshape(ma.shape) for x in [a_s, b_s, null_a_s, null_b_s])
            
    result = []
    for a,b,null_a,null_b in bootstrap_samples(ma, mb, bootstrap_size, num_samples):
        result.append(
            (scipy.stats.entropy(softmax(a).flatten(),softmax(b).flatten()),
             scipy.stats.entropy(softmax(null_a).flatten(),softmax(null_b).flatten())))
    return pd.DataFrame(result, columns=['test','null'])

def sampleAndPlot(m1, m2, title, bootstrap_size=BOOTSTRAP_SIZE,
                  total_samples=TOTAL_SAMPLES, num_workers=NUM_WORKERS):
    futures = []
    with concurrent.futures.ProcessPoolExecutor(num_workers) as executor:
        for s in range(0,TOTAL_SAMPLES,(total_samples // num_workers)):
            futures.append(
                executor.submit(
                    bootstrap_test,
                    m1, m2, bootstrap_size, (total_samples // num_workers)))
    df_direct = pd.concat([future.result() for future in futures])
    return ggplot(df_direct, aes()) \
        + xlab('Summary Statistic Distribution') \
        + ggtitle('%s\n(Null Distribution shown in gray)' % title) \
        + geom_histogram(aes(x='null'), fill='black', color='black', alpha=0.5, binwidth=0.001) \
        + geom_histogram(aes(x='test'), fill='red', color='red', alpha=0.5, binwidth=0.001)

def normalize(m):
    return m.astype(np.float32)/np.sum(m)

def findCategoryEdges(facL, facDict):
    """This assumes that the facility abbrev list facL is sorted by category"""
    edgeD = {}
    oCat = facDict[facL[0]]['category']
    edgeD[oCat.lower() + 'L'] = 0
    for idxM1, abbrev in enumerate(facL[1:]):
        idx = idxM1 + 1 # we want the un-shifted index
        if abbrev in facDict:
            nCat = facDict[abbrev]['category']
        else:
            nCat = abbrev
        if nCat != oCat:
            edgeD[oCat.lower() + 'H'] = idx - 1
            edgeD[nCat.lower() + 'L'] = idx
            oCat = nCat
    edgeD[nCat.lower() + 'H'] = idx

    return edgeD

def partition(m, facL, facDict):
    """matrix to partition, ordered list of facilities, facility dictionary"""
    edgeD = findCategoryEdges(facL, facDict)
    hospL = edgeD['hospitalL']
    hospH = edgeD['hospitalH']
    if 'ltacL' in edgeD:
        ltacL = edgeD['ltacL']
        ltacH = edgeD['ltacH']
    else:
        ltacL = edgeD['ltachL']
        ltacH = edgeD['ltachH']
    if 'snfL' in edgeD:
        snfL = edgeD['snfL']
        snfH = edgeD['snfH']
    else:
        snfL = edgeD['nursinghomeL']
        snfH = edgeD['nursinghomeH']
    vsnfL = edgeD['vsnfL']
    vsnfH = edgeD['vsnfH']
    mm = np.zeros([m.shape[0], 4])
    mm[:,0] = m[:,hospL:hospH].sum(1)
    mm[:,1] = m[:,ltacL:ltacH].sum(1)
    mm[:,2] = m[:,snfL:snfH].sum(1)
    mm[:,3] = m[:,vsnfL:vsnfH].sum(1)
    return mm

def myPie(mx, rowIdx, ax, idxFacTbl, facDict):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    normRw = normalize(mx[rowIdx,:])
    catClrMap = {ctg: colors[i] 
                 for i,ctg in enumerate(['HOSPITAL', 'LTACH', 'VSNF', 'SNF', 'COMMUNITY']) }
    segTbl = {key: [] for key in catClrMap.keys()}
    labels = []
    sizes = []
    cL = []
    for idx, val in enumerate(normRw.flatten()):
        abbrev = idxFacTbl[idx]
        if abbrev == 'COMMUNITY':
            category = 'COMMUNITY'
        else:
            category = facDict[abbrev]['category']
        segTbl[category].append((val, abbrev))
    kL = segTbl.keys()
    kL.sort()
    for category in kL:
        segL = segTbl[category]
        segL.sort(reverse=True)
        totV = sum([val for val, abbrev in segL])
        szLim = 0.05 * totV
        otherSz = 0.0
        otherCt = 0
        curClr = catClrMap[category]
        for val, abbrev in segL:
            if val > szLim:
                labels.append(abbrev)
                sizes.append(val)
                cL.append(curClr)
            elif val > 0.0:
                otherSz += val
                otherCt += 1
        if otherCt > 0:
            labels.append('%d Other %ss' % (otherCt, category[:4]))
            sizes.append(otherSz)
            cL.append(curClr)
    tpl = ax.pie(sizes, labels=labels, colors=cL, startangle=90)
    patches, texts = tpl
    for txt in texts:
        txt.set_fontsize(7)

def drawAllFour(abbrev, measDirectMx, measIndMx, simDirectMx, simIndMx,
                idxFacTbl, facDict):
    idx = {val:key for key, val in idxFacTbl.items()}[abbrev]
    
    outerFig = plt.figure()
    outerAx = outerFig.gca()
    outerAx.set_title(abbrev)
    outerAx.set_aspect('equal', 'datalim')
    outerAx.set_ylim(0.0,4.0)
    outerAx.set_yticks([1.0, 3.0])
    outerAx.set_yticklabels(['Measured', 'Simulated'])
    outerAx.set_xlim(0.0,4.0)
    outerAx.set_xticks([1.0, 3.0])
    outerAx.set_xticklabels(['Direct', 'Indirect'])

    ax1 = inset_axes(outerAx, width='40%', height='40%', loc=2, borderpad=1.0)
    ax1.set_aspect('equal', 'datalim')
    myPie(simDirectMx, idx, ax1, idxFacTbl, facDict)

    ax2 = inset_axes(outerAx, width='40%', height='40%', loc=3, borderpad=1.0)
    ax2.set_aspect('equal', 'datalim')
    myPie(measDirectMx, idx, ax2, idxFacTbl, facDict)

    ax3 = inset_axes(outerAx, width='40%', height='40%', loc=1, borderpad=1.0)
    ax3.set_aspect('equal', 'datalim')
    myPie(simIndMx, idx, ax3, idxFacTbl, facDict)

    ax4 = inset_axes(outerAx, width='40%', height='40%', loc=4, borderpad=1.0)
    ax4.set_aspect('equal', 'datalim')
    myPie(measIndMx, idx, ax4, idxFacTbl, facDict)

def myBars_old(mxList, rowIdx, ax, colLabelL, idxFacTbl, facDict):
    assert len(colLabelL) == len(mxList), 'number of labels does not match number of columns'
    facIdxTbl = {val: key for key, val in idxFacTbl.items()}
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    # cannot use a completion here because of semantics of generator
    normRwL = []
    for mx in mxList:
        normRwL.append(normalize(mx[rowIdx,:]))
    catClrMap = {ctg: colors[i] for i,ctg in enumerate(['HOSPITAL', 'LTACH', 'VSNF', 'SNF', 'COMMUNITY']) }
    segTbl = {key: [] for key in catClrMap.keys()}
    for idx, val in enumerate(normRwL[0].flatten()):
        abbrev = idxFacTbl[idx]
        if abbrev == 'COMMUNITY':
            category = 'COMMUNITY'
        else:
            category = facDict[abbrev]['category']
        segTbl[category].append((val, abbrev))
    kL = segTbl.keys()
    kL.sort()
    labels = []
    sizeLL = []  # will be a list of lists
    cL = []
    for category in kL:
        segL = segTbl[category]
        segL.sort(reverse=True)
        totV = sum([val for val, abbrev in segL])
        szLim = 0.05 * totV
        otherSz = 0.0
        otherCt = 0
        curClr = catClrMap[category]
        otherSzL = [0.0 for _ in normRwL]
        for val, abbrev in segL:
            idx = facIdxTbl[abbrev] if abbrev in facIdxTbl else facIdxTbl['COMMUNITY']
            szL = [row[idx] for row in normRwL]
            if val > szLim:
                labels.append(abbrev)
                sizeLL.append(szL)
                cL.append(curClr)
            elif val > 0.0:
                otherSzL = [a + b for a, b in zip(szL, otherSzL)]
                otherCt += 1
        if otherCt > 0:
            labels.append('%d Other %ss' % (otherCt, category[:4]))
            sizeLL.append(otherSzL)
            cL.append(curClr)
    xL = np.arange(len(normRwL)) + 2
    baseL = [0.0 for _ in normRwL]
    for offset in xrange(len(normRwL) - 1):
        ax.plot([2.8 + offset, 3.0 + offset],[0.0, 0.0], color='k')
    for szL,clr,lbl in zip(sizeLL, cL, labels):
        artists = ax.bar(xL, szL, color=clr, bottom=baseL)
        ax.text(1.8, baseL[0] + 0.5*szL[0], lbl,
                horizontalalignment='right', verticalalignment='center', fontsize=6)
        #print lbl, xL, yL, clr, baseL
        baseL = [a + b for a, b in zip(szL, baseL)]
        for offset in xrange(len(normRwL) - 1):
            ax.plot([2.8 + offset, 3.0 + offset],baseL[offset:offset+2], color='k')
    ax.set_xlim((0.0, len(mxList) + 2.0))
    ax.set_xticks(np.arange(len(mxList)) + 2.4)
    ax.set_xticklabels(colLabelL)
    ax.set_yticks([])
    
def myBars(mxList, rowIdx, ax, colLabelL, idxFacTbl, facDict,
           hideThisIdx=None):
    assert len(colLabelL) == len(mxList), 'number of labels does not match number of columns'
    facIdxTbl = {val: key for key, val in idxFacTbl.items()}
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    # cannot use a completion here because of semantics of generator
    normRwL = []
    totRwL = []
    for mx in mxList:
        rw = mx[rowIdx,:]
        if hideThisIdx is not None:
            rw[hideThisIdx] = 0
        normRwL.append(normalize(rw))
        totRw = np.sum(mx[rowIdx,:])
        if hideThisIdx is not None:
            totRw -= mx[rowIdx, hideThisIdx]
        totRwL.append(totRw)
    catClrMap = {ctg: colors[i] for i,ctg
                 in enumerate(['HOSPITAL', 'LTACH', 'VSNF', 'SNF', 'COMMUNITY']) }
    segTbl = {key: [] for key in catClrMap.keys()}
    for idx, val in enumerate(normRwL[0].flatten()):
        abbrev = idxFacTbl[idx]
        if abbrev == 'COMMUNITY':
            category = 'COMMUNITY'
        else:
            category = facDict[abbrev]['category']
        if idx != hideThisIdx:
            segTbl[category].append((val, abbrev))
    kL = segTbl.keys()
    kL.sort()
    labels = []
    sizeLL = []  # will be a list of lists
    cL = []
    for category in kL:
        segL = segTbl[category]
        segL.sort(reverse=True)
        totV = sum([val for val, abbrev in segL])
        szLim = 0.05 * totV
        otherCt = 0
        curClr = catClrMap[category]
        otherSzL = [0.0 for _ in normRwL]
        for val, abbrev in segL:
            idx = facIdxTbl[abbrev] if abbrev in facIdxTbl else facIdxTbl['COMMUNITY']
            szL = [row[idx] for row in normRwL]
            if szL[0] > szLim:
                labels.append(abbrev)
                sizeLL.append(szL)
                cL.append(curClr)
            else:
                otherSzL = [a + b for a, b in zip(szL, otherSzL)]
                otherCt += 1
        if otherCt > 0:
            labels.append('%d Other %ss' % (otherCt, category[:4]))
            sizeLL.append(otherSzL)
            cL.append(curClr)
    ax.set_xlim(0.0, 3.0 + len(normRwL))
    xL = np.arange(len(normRwL)) + 2
    zOff = 0.1
    baseL = [zOff for _ in normRwL]
    for offset in xrange(len(normRwL) - 1):
        ax.plot([2.8 + offset, 3.0 + offset],[zOff, zOff], color='k')
    for offset in xrange(len(normRwL)):
        colLbl = colLabelL[offset]
        if len(colLbl) > 12:
            colLbl = '...' + colLbl[-10:]
        ax.text(2.0 + offset, 0.0, colLbl,
                horizontalalignment='center', verticalalignment='bottom',
                fontsize=6)
    for szL,clr,lbl in zip(sizeLL, cL, labels):
        rects = ax.bar(xL, szL, color=clr, edgecolor='k', bottom=baseL)
        for sz, tot, rect in zip(szL, totRwL, rects):
            yOff = (0.5 * rect.get_height()) + rect.get_y()
            if sz > 0.0:
                ax.text(rect.get_x() + rect.get_width()/2., yOff,
                        '%d' % int(round(sz*tot)),
                        ha='center', va='center')
        ax.text(1.5, baseL[0] + 0.5*szL[0], lbl,
                horizontalalignment='right', verticalalignment='center', fontsize=6)
        ax.text(1.5+len(normRwL), baseL[-1] + 0.5*szL[-1], lbl,
                horizontalalignment='left', verticalalignment='center', fontsize=6)
        #print lbl, xL, yL, clr, baseL
        baseL = [a + b for a, b in zip(szL, baseL)]
        for offset in xrange(len(normRwL) - 1):
            ax.plot([2.4 + offset, 2.6 + offset],baseL[offset:offset+2], color='k')
    # Draw totals above each column, using the last set of rectangles to get coordinates
    for tot, base, rect in zip(totRwL, baseL, rects):
        ax.text(rect.get_x() + rect.get_width()/2.0,
                base,
                '%d' % tot, ha='center', va='bottom')
    ax.set_xticks([])
    #ax.set_xticks(np.arange(len(normRwL)+2))
    ax.set_yticks([])

def drawBarSets(abbrev, directMxList, indirectMxList, colLabels, idxFacTbl, facDict,
                hideDirectTransfersToSelf=False):
    facIdxTbl = {val: key for key, val in idxFacTbl.items()}
    idx = facIdxTbl[abbrev]

    figs, axisL = plt.subplots(nrows=2, ncols=1)
    axisL[0].set_title(abbrev)
    if hideDirectTransfersToSelf:
        hideThisIdx = idx
    else:
        hideThisIdx = None
    myBars(directMxList, idx, axisL[0], colLabels, idxFacTbl, facDict,
           hideThisIdx=hideThisIdx)
    axisL[0].axis('off')
    myBars(indirectMxList, idx, axisL[1], colLabels, idxFacTbl, facDict,
           hideThisIdx=None)
    axisL[1].axis('off')
    figs.tight_layout(h_pad=2)
    figs.canvas.set_window_title(abbrev)

def plotLogImg(mtx, minVal, maxVal):
    fltIm = np.log10(mtx.astype(np.float32) + 1.0)
    myMin = np.log10(float(minVal) + 1.0)
    myMax = np.log10(float(maxVal) + 1.0)
#     fltIm = mtx.astype(np.float32)
#     myMin, myMax = minVal, maxVal
    return plt.imshow(fltIm,
#                       cmap=cm.RdYlGn,
                      cmap=cm.viridis,
                      vmin=myMin, vmax=myMax,
#                       interpolation='bilinear',
#                       origin='lower',
#                       extent=[-3, 3, -3, 3]
               )

def plotImg(mtx, minVal, maxVal):
    fltIm = mtx.astype(np.float32)
    myMin, myMax = minVal, maxVal
    return plt.imshow(fltIm,
                      cmap=cm.bwr,
                      vmin=myMin, vmax=myMax,
                      )

def pltMtxImageFig(directMtx, measDirectMtx, indirectMtx, measIndirectMtx):
    maxVal = max(np.max(directMtx), np.max(indirectMtx))

    ax11 = plt.subplot(2, 2, 1)
    ax11.set_title('log direct')
    plotLogImg(directMtx, 0.0, float(maxVal))

    ax12 = plt.subplot(2, 2, 2)
    ax12.set_title('log indirect')
    plotLogImg(indirectMtx, 0.0, float(maxVal))

    plt.colorbar(ax=[ax11, ax12])

    sclMtx = (float(np.sum(directMtx))/float(np.sum(measDirectMtx))) * measDirectMtx
    #deltaMtx = (2.0*(directMtx - sclMtx)/(directMtx + sclMtx))
    directDeltaMtx = directMtx - sclMtx
    #directDeltaMtx[0:100, 0:100] = 0.0

    sclMtx = (float(np.sum(indirectMtx))/float(np.sum(measIndirectMtx))) * measIndirectMtx
    #deltaMtx = (2.0*(indirectMtx - sclMtx)/(directMtx + sclMtx))
    indirectDeltaMtx = indirectMtx - sclMtx
    #indirectDeltaMtx[0:100, 0:100] = 0.0
    lim = max(np.max(np.fabs(directDeltaMtx)), np.max(np.fabs(indirectDeltaMtx)))

    ax21 = plt.subplot(2, 2, 3)
    ax21.set_title('normalized direct delta')
    plotImg(directDeltaMtx, -lim, lim)

    ax22 = plt.subplot(2, 2, 4)
    ax22.set_title('normalized indirect delta')
    plotImg(indirectDeltaMtx, -lim, lim)

    plt.colorbar(ax=[ax21, ax22])


def main():
    input_yaml = '../sim/twoyear_run_ChicagoLand.yaml'
    allD = {}
    for dataCase in ['indirect_0_830_case0', 'indirect_0_830_case1',
                     'indirect_100_465_case0', 'indirect_100_465_case1',
                     'indirect_466_830_case0', 'indirect_466_830_case1',
                     'indirect_0_465', 'directonly_0_465',
                     'directonly_0_830_case0', 'directonly_0_830_case1',
                     'directonly_100_465_case0', 'directonly_100_465_case1',
                     'directonly_466_830_case0', 'directonly_466_830_case1',
                     'indirect_466_830_42abc723_case0', 'indirect_466_830_42abc723_case1',
                     'indirect_466_830_42abc723_case2', 'indirect_466_830_42abc723_case3',
                     'indirect_466_830_42abc723_case4', 'indirect_0_465_2a2a0ef',
                     'indirect_0_830_2a2a0ef', 'indirect_466_830_2a2a0ef'
                    ]:
        d = np.load('arrays_%s.npz' % dataCase)
        allD[dataCase] = d

    print os.path.abspath(input_yaml)
    schemautils.setSchemaBasePath(SCHEMA_DIR)
    inputDict = checkInputFileSchema(input_yaml,
                                     #os.path.join(SCHEMA_DIR,'rhea_input_schema.yaml'),
                                     'rhea_input_schema.yaml',
                                     comm=None)
    if 'modelDir' in inputDict:
        pyrheautils.PATH_STRING_MAP['MODELDIR'] = pyrheautils.pathTranslate(inputDict['modelDir'])
    if 'pathTranslations' in inputDict:
        for elt in inputDict['pathTranslations']:
            pyrheautils.PATH_STRING_MAP[elt['key']] = elt['value']
    facDirs = [pyrheautils.pathTranslate(dct) for dct in inputDict['facilityDirs']]
    pathPrefix = os.path.dirname(os.path.abspath(input_yaml))
    facDirs = [os.path.join(pathPrefix, fD) for fD in facDirs]
    facDict = parseFacilityData(facDirs)
    print 'IMPLEMENTING SPECIAL PATCH FOR WAUK_2615_H'
    facDict['WAUK_2615_H'] = {'category':'HOSPITAL'}
    print 'FIND OUT THE REAL ANSWER AND DELETE THIS!'
    
    # Sort by category, then alphabetically by name
    facL = [(facDict[fac]['category'], fac) for fac in facDict.keys()
            if facDict[fac]['category'] != 'COMMUNITY']
    facL.sort()
    facL = [fac for facCat, fac in facL]
    facIdxTbl = {key:idx for idx, key in enumerate(facL)}
    idxFacTbl = {idx:key for key, idx in facIdxTbl.items()}
    idxFacTbl[len(facL)] = 'COMMUNITY'

    # lclRates for FRAN_20201_H: {'death': 0.0, 'nursinghome': 0.10150579568844112,
    #                             'hospital': 0.048423789405264865, 'ltac': 0.003033257501895786,
    #                             'vsnf': 0.021828620951142887, 'home': 0.8252085364532553}
    #lclRates for FRAN_1423_H: {'death': 0.0, 'nursinghome': 0.13148148148148148,
    #                           'hospital': 0.19074074074074074, 'ltac': 0.0, 'vsnf': 0.0,
    #                           'home': 0.6777777777777778}

    mMx = allD['directonly_466_830_case0']['direct_measured']
    miMx = allD['directonly_466_830_case0']['indirect_measured']
    sMx0 = allD['directonly_466_830_case0']['direct_simulated']
    siMx0 = allD['directonly_466_830_case0']['indirect_simulated']
    sMx1 = allD['directonly_466_830_case1']['direct_simulated']
    siMx1 = allD['directonly_466_830_case1']['indirect_simulated']

    # The following is done as a test of partition
    print 'partition test matrix:'
    print partition(mMx, facL, facDict)

    #sampleAndPlot(mMx, sMx0, 'test of sampleAndPlot').show()

    drawAllFour('FRAN_20201_H', mMx, miMx, sMx0+sMx1, siMx0+siMx1, idxFacTbl, facDict)
    drawAllFour('FRAN_1423_H', mMx, miMx, sMx0+sMx1, siMx0+siMx1, idxFacTbl, facDict)

    caseStr = 'indirect_466_830_42abc723'
    drawBarSets('PRES_333_H',
               [mMx, allD[caseStr+'_case0']['direct_simulated'], allD[caseStr+'_case1']['direct_simulated'],
                allD[caseStr+'_case2']['direct_simulated'], allD[caseStr+'_case3']['direct_simulated']
               ],
               [miMx, allD[caseStr+'_case0']['indirect_simulated'], allD[caseStr+'_case1']['indirect_simulated'],
                allD[caseStr+'_case2']['indirect_simulated'], allD[caseStr+'_case3']['indirect_simulated']
               ],
               ['measured', 'run0', 'run1', 'run2', 'run3'],
               idxFacTbl, facDict)
    
    pltMtxImageFig(sMx0, mMx, siMx0, miMx)

    plt.show()

if __name__ == "__main__":
    main()
