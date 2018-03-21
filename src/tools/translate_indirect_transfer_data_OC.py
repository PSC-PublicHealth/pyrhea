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
indirect transfers in Orange County.
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

inputCSV = '../../models/OrangeCounty2013/OC_Indirect_Transfer_Matrix_for_RHEA_2.0_-_Adult_Only_-_UPDATED_-_10-20-2017_HOSP_INDIRECT.csv'

srcColKey = '"Hosp-Hosp Indirect Transfers (2013 followed 365 days)"'
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

    with open(os.path.join(os.path.dirname(__file__), inputCSV), 'rU') as f:
        keys, recs = csv_tools.parseCSV(f)
    hospTbl = {}
    nhTbl = {}
    hospCatList = ['HOSPITAL', 'LTAC', 'LTACH']
    toSnfCt = 0
    toHospCt = 0
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
        elif facDict[src]['category'] in hospCatList:
            tbl = hospTbl
        else:
            tbl = nhTbl
        assert src not in tbl, 'Duplicate entry for source %s' % src
        tbl[src] = {}
        for dst in [k for k in keys if k not in ignoredKeys]:
            n = rec[dst]
            if n != 0:
                assert str(dst) not in tbl[src], 'redundant entry for %s %s' % (src, dst)
                #print '%s -> %s' % (src, dst)
                tbl[src][str(dst)] = n
                if dst not in facDict:
                    lostTrans = rec[dst]
                    print 'Unknown destination %s: %d transfers lost' % (dst, lostTrans)
                    lostCt += lostTrans
                    continue
                if facDict[src]['category'] in hospCatList:
                    #print '%s -> %s' % (facDict[src]['category'], facDict[str(dst)]['category'])
                    if facDict[str(dst)]['category'] in hospCatList:
                        toHospCt += n
                    else:
                        toSnfCt += n
    print 'hosp -> hosp %s' % toHospCt
    print 'hosp -> snf %s' % toSnfCt
    print 'total transfers lost to unknown endpoints: %d' % lostCt
    print 'parsed'
    with open('hosp_indirect_transfer_counts.yaml', 'w') as f:
        yaml.dump(hospTbl, f, indent=4, encoding='utf-8',width=130,explicit_start=True)
    with open('nh_readmit_transfer_counts.yaml', 'w') as f:
        yaml.dump(nhTbl, f, indent=4, encoding='utf-8',width=130,explicit_start=True)
    print 'done'


if __name__ == "__main__":
    main()
