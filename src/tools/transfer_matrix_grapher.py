#! /usr/bin/env python

"""
Generate diagnostic graphs for measured vs. simulated transfer matrices
"""

from optparse import OptionParser
import glob
import os.path

import matplotlib.pyplot as plt
import numpy as np

import schemautils
import pyrheautils
import tools_util as tu
from plotting_utilities import drawBarSets

EXPECTED_MTX_KEYS = frozenset(['direct_measured', 'indirect_measured',
                               'direct_simulated', 'indirect_simulated'])

def loadMtxFile(fname):
    """
    Make sure the file is a .npz file of the sort produced by arrival_time_parser,
    and if so return its contents as a dict.
    """
    assert fname.endswith('.npz'), '%s is not a numpy matrix file?' % fname
    mtxD = np.load(fname)
    assert set(mtxD.keys()) == EXPECTED_MTX_KEYS, '%s has the wrong matrix names' % fname
    return mtxD

def expandGlobbedList(pathList):
    """
    Given a list of path strings, expand with globbing
    """
    newPathList = []
    for fName in pathList:
        print '%s yields %s' % (fName, glob.glob(fName))
        newPathList += glob.glob(fName)
    return newPathList

def main():
    """Main"""
    parser = OptionParser(usage="""
    %prog [--matrix notes_file.pkl] [--glob] [-f abbrev] run_descr.yaml
    """)
    parser.add_option('-m', '--matrix', action='append', type='string',
                      help="matrix filename - may be repeated.  These are numpy .npz files")
    parser.add_option('--glob', action='store_true',
                      help="Apply filename globbing to the given matrix files")

    parser.add_option('-f', '--facility', action='append', type='string',
                      help=('Facility abbrev to plot - may be repeated.  The default is to'
                            'use the trackedFacilities list from the input file'))
    
    parser.add_option('--transpose', action='store_true',
                      help=('Graph the transpose of the transfer matrix, thus showing inputs rather'
                            ' than outputs'))

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    if not opts.matrix:
        parser.error('At least one --matrix option is required')

    parser.destroy()

    runDesc = args[0]

    print 'Reading run info from %s' % os.path.abspath(runDesc)
    inputDict = tu.readModelInputs(args[0])
    facDict = tu.getFacDict(inputDict)

    print 'IMPLEMENTING SPECIAL PATCH FOR WAUK_2615_H'
    facDict['WAUK_2615_H'] = {'category':'HOSPITAL'}
    print 'FIND OUT THE REAL ANSWER AND DELETE THIS!'

    # Sort by category, then alphabetically by name
    facL = [(facDict[fac]['category'], fac) for fac in facDict
            if facDict[fac]['category'] != 'COMMUNITY']
    facL.sort()
    facL = [fac for facCat, fac in facL]  # @UnusedVariable
    facIdxTbl = {key:idx for idx, key in enumerate(facL)}
    idxFacTbl = {idx:key for key, idx in facIdxTbl.items()}
    idxFacTbl[len(facL)] = 'COMMUNITY'

    mtxFileL = opts.matrix
    if opts.glob:
        mtxFileL = expandGlobbedList(mtxFileL)
    allD = {os.path.splitext(os.path.basename(fname))[0]: loadMtxFile(fname)
            for fname in mtxFileL}
    mtxNameL = [os.path.splitext(os.path.basename(fname))[0]
                for fname in mtxFileL]

    sampM = None
    for key, mD in allD.items():
        testM = mD['direct_measured']
        assert testM.shape == (len(facL) + 1, len(facL) + 1), ('%s does not match known facilities'
                                                               % key)
        if sampM is None:
            sampM = testM
        else:
            assert np.allclose(testM,sampM), ('Not all matrices represent the same'
                                              ' direct transfer data table')

    if opts.facility:
        plotThese = opts.facility
    else:
        plotThese = inputDict['trackedFacilities']

    assert plotThese, 'No facilities to plot'

    directL = [allD[mtxNameL[0]]['direct_measured']]
    indirectL = [allD[mtxNameL[0]]['indirect_measured']]
    labelL = ['measured']
    for mtxName in mtxNameL:
        directL.append(allD[mtxName]['direct_simulated'])
        indirectL.append(allD[mtxName]['indirect_simulated'])
        labelL.append(mtxName)

    if opts.transpose:
        directL = [mtx.transpose() for mtx in directL]
        indirectL = [mtx.transpose() for mtx in indirectL]

    for abbrev in plotThese:
        drawBarSets(abbrev, directL, indirectL, labelL, idxFacTbl, facDict,
                    hideDirectTransfersToSelf=True,
                    titleSuffix=('inputs' if opts.transpose else 'outputs'))
    plt.show()

if __name__ == "__main__":
    main()
