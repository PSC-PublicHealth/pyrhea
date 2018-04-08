#! /usr/bin/env python

###################################################################################
# Copyright   2017, Pittsburgh Supercomputing Center (PSC).  All Rights Reserved. #
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

"""
This program parses a CSV file extracted from an xlsx representation of a transfer 
matrix and outputs the transfer data in yaml form.  This version is designed for
direct transfers.
"""

import sys
import os.path
import signal
import optparse
import yaml
from collections import defaultdict
import numpy as np
from scipy.stats import norm
import logging
import logging.config

import schemautils
import pyrheautils
from pyrhea import getLoggerConfig, checkInputFileSchema, loadPathogenImplementations
from facilitybase import CareTier, PthStatus
from map_transfer_matrix import parseFacilityData

import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools

SCHEMA_DIR = os.path.join(os.path.dirname(__file__), '../schemata')
INPUT_SCHEMA = 'rhea_input_schema.yaml'

LOGGER = None

inputCSV = '../../models/ChicagoLand/Matrices_LOS_09292016_cleaned_Transfer3day.csv'

def main():
    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)

    global LOGGER
    #logging.config.dictConfig(getLoggerConfig())
    LOGGER = logging.getLogger(__name__)

    parser = optparse.OptionParser(usage="""
    %prog run_descr.yaml
    """)

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('An input yaml file matching %s is required' % INPUT_SCHEMA)

    parser.destroy()

    inputDict = checkInputFileSchema(args[0],
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    schemautils.setSchemaBasePath(SCHEMA_DIR)
    if 'modelDir' in inputDict:
        pyrheautils.PATH_STRING_MAP['MODELDIR'] = pyrheautils.pathTranslate(inputDict['modelDir'])
    if 'pathTranslations' in inputDict:
        for elt in inputDict['pathTranslations']:
            pyrheautils.PATH_STRING_MAP[elt['key']] = elt['value']

    facDirs = [pyrheautils.pathTranslate(dct) for dct in inputDict['facilityDirs']]
    facDict = parseFacilityData(facDirs)

    with open(os.path.join(os.path.dirname(__file__), inputCSV), 'rU') as f:
        keys, recs = csv_tools.parseCSV(f)
    allTbl = {}
    hospCatList = ['HOSPITAL', 'LTAC', 'LTACH']
    toSnfCt = 0
    toHospCt = 0
    allCt = 0
    lostCt = 0
    selfCt = 0
    for rec in recs:
        src = str(rec['UNIQUE_ID'])
        assert src not in allTbl, 'Duplicate entry for source %s' % src
        allTbl[src] = {}
        for dst in [k for k in keys if k != 'UNIQUE_ID']:
            n = rec[dst]
            if n != 0:
                assert str(dst) not in allTbl[src], 'redundant entry for %s %s' % (src, dst)
                if dst == src:
                    selfCt += n
                else:
                    allTbl[src][str(dst)] = n
                    if dst not in facDict:
                        lostTrans = n
                        if n != 0:
                            print 'Unknown destination %s: %d transfers lost' % (dst, lostTrans)
                        lostCt += lostTrans
                    else:
                        if src not in facDict:
                            lostTrans = n
                            if n != 0:
                                print 'Unknown source %s: %d transfers lost' % (src, lostTrans)
                            lostCt += lostTrans
                        else:
                            if facDict[src]['category'] in hospCatList:
                                print '%s -> %s' % (facDict[src]['category'],
                                                    facDict[str(dst)]['category'])
                                if facDict[str(dst)]['category'] in hospCatList:
                                    toHospCt += n
                                else:
                                    toSnfCt += n
                        allCt += n
    print 'understood %s transfers' % allCt
    print 'excluded %s self transfers' % selfCt
    print 'total transfers lost to unknown endpoints: %d' % lostCt
    print 'hosp -> hosp %s' % toHospCt
    print 'hosp -> snf %s' % toSnfCt
    print 'parsed'
    with open('direct_transfer_counts.yaml', 'w') as f:
        yaml.dump(allTbl, f, indent=4, encoding='utf-8',width=130,explicit_start=True)
    print 'done'


if __name__ == "__main__":
    main()
