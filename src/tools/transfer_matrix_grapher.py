#! /usr/bin/env python

from optparse import OptionParser
import glob

import matplotlib.pyplot as plt
import numpy as np
import os.path

import schemautils
import pyrheautils
from plotting_utilities import drawBarSets
from pyrhea import checkInputFileSchema
from map_transfer_matrix import parseFacilityData

SCHEMA_DIR = os.path.join(os.path.dirname(__file__), '../schemata')

EXPECTED_MTX_KEYS = frozenset(['direct_measured', 'indirect_measured',
                               'direct_simulated', 'indirect_simulated'])

def loadMtxFile(fname):
    """
    Make sure the file is a .npz file of the sort produced by arrival_time_parser,
    and if so return its contents as a dict.
    """
    assert fname.endswith('.npz'), '%s is not a numpy matrix file?' % fname
    d = np.load(fname)
    assert set(d.keys()) == EXPECTED_MTX_KEYS, '%s has the wrong matrix names' % fname
    return d

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

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    if not opts.matrix:
        parser.error('At least one --matrix option is required')

    parser.destroy()

    runDesc = args[0]

    print 'Reading run info from %s' % os.path.abspath(runDesc)
    #
    # Begin boilerplate to import run information
    #
    schemautils.setSchemaBasePath(SCHEMA_DIR)
    inputDict = checkInputFileSchema(runDesc,
                                     'rhea_input_schema.yaml',
                                     comm=None)
    if 'modelDir' in inputDict:
        pyrheautils.PATH_STRING_MAP['MODELDIR'] = pyrheautils.pathTranslate(inputDict['modelDir'])
    if 'pathTranslations' in inputDict:
        for elt in inputDict['pathTranslations']:
            pyrheautils.PATH_STRING_MAP[elt['key']] = elt['value']
    facDirs = [pyrheautils.pathTranslate(dct) for dct in inputDict['facilityDirs']]
    pathPrefix = os.path.dirname(os.path.abspath(runDesc))
    facDirs = [os.path.join(pathPrefix, fD) for fD in facDirs]
    facDict = parseFacilityData(facDirs)
    #
    # End boilerplate to import run information
    #
    
    print 'IMPLEMENTING SPECIAL PATCH FOR WAUK_2615_H'
    facDict['WAUK_2615_H'] = {'category':'HOSPITAL'}
    print 'FIND OUT THE REAL ANSWER AND DELETE THIS!'

    # Sort by category, then alphabetically by name
    facL = [(facDict[fac]['category'], fac) for fac in facDict.keys()
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

    mtxNameL = allD.keys()[:]
    mtxNameL.sort()
    directL = [allD[mtxNameL[0]]['direct_measured']]
    indirectL = [allD[mtxNameL[0]]['indirect_measured']]
    labelL = ['measured']
    for mtxName in mtxNameL:
        directL.append(allD[mtxName]['direct_simulated'])
        indirectL.append(allD[mtxName]['indirect_simulated'])
        labelL.append(mtxName)

    for abbrev in plotThese:
        drawBarSets(abbrev, directL, indirectL, labelL, idxFacTbl, facDict)
    plt.show()

if __name__ == "__main__":
    main()
