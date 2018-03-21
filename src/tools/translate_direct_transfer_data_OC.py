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
This program parses a set of CSV files extracted from an xlsx representation of a transfer 
matrix and outputs the transfer data in yaml form.  This version is designed for
direct transfers in Orange County.
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
import tools_util as tu
from facilitybase import CareTier, PthStatus

import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools

SCHEMA_DIR = os.path.join(os.path.dirname(__file__), '../schemata')
INPUT_SCHEMA = 'rhea_input_schema.yaml'

LOGGER = None

inputDir = '../../models/OrangeCounty2013'
inputCSVL = ['OC_Direct_Transfer_Matrices_for_RHEA_2.0_-_Adult_Only_-_09-15-2017_FINAL_HOSP-HOSP.csv',
             'OC_Direct_Transfer_Matrices_for_RHEA_2.0_-_Adult_Only_-_09-15-2017_FINAL_NH-NH.csv',
             'OC_Direct_Transfer_Matrices_for_RHEA_2.0_-_Adult_Only_-_09-15-2017_FINAL_HOSP-NH.csv',
             'OC_Direct_Transfer_Matrices_for_RHEA_2.0_-_Adult_Only_-_09-15-2017_FINAL_NH-HOSP.csv']


srcColKey = '"Transfers 2 Years of Data Averaged (2013-2014)"'
ignoredKeys = [srcColKey, 'TOTAL']

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

    runDesc = args[0]
    inputDict = tu.readModelInputs(runDesc)
    #print inputDict.keys()
    facDict = tu.getFacDict(inputDict)

    allKeyS = set()
    allRecs = []
    for inputCSV in inputCSVL:
        with open(os.path.join(os.path.dirname(__file__), inputDir, inputCSV), 'rU') as f:
            keys, recs = csv_tools.parseCSV(f)
            allKeyS = allKeyS.union(keys)
            allRecs.extend(recs)
    keys = list(allKeyS)
    recs = allRecs
    allTbl = {}
    hospCatList = ['HOSPITAL', 'LTAC', 'LTACH']
    allCt = 0
    lostCt = 0
    for rec in recs:
        src = str(rec[srcColKey])
        if src in ignoredKeys:
            continue
        elif src not in facDict:
            lostTrans = sum([rec[k] for k in rec if k not in ignoredKeys])
            print 'unknown source %s: %d transfers lost' % (src, lostTrans)
            lostCt += lostTrans
            continue
        else:
            tbl = allTbl
        if src not in tbl:
            tbl[src] = {}
        for dst in [k for k in keys if k not in ignoredKeys]:
            assert dst.startswith('To_'), 'unexpected dst format for %s' % dst
            dst = dst[3:]
            n = rec['To_' + dst] if 'To_' + dst in rec else 0
            if n != 0:
                assert str(dst) not in tbl[src], 'redundant entry for %s %s' % (src, dst)
                #print '%s -> %s' % (src, dst)
                tbl[src][str(dst)] = n
                if dst not in facDict:
                    lostTrans = n
                    if n != 0:
                        print 'Unknown destination %s: %d transfers lost' % (dst, lostTrans)
                    lostCt += lostTrans
                    continue
                allCt += n
    print 'understood %s transfers' % allCt
    print 'total transfers lost to unknown endpoints: %d' % lostCt
    print 'parsed'
    with open('direct_transfer_counts.yaml', 'w') as f:
        yaml.dump(allTbl, f, indent=4, encoding='utf-8',width=130,explicit_start=True)
    print 'done'


if __name__ == "__main__":
    main()
