#! /usr/bin/env python

"""
This script applies patches to specific facility files.
"""

import os.path
import types
import yaml
from StringIO import StringIO
import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools

def loadCSVByAbbrev(modelDir, fname, key=None):
    if key is None:
        key = 'abbrev'
    fullName = os.path.join(modelDir, fname)
    with open(fullName) as fl:
        keys, recs = csv_tools.parseCSV(fl)
    assert key in keys, ('%s has no "%s" field' % (fullName, key))
    return {rec[key]: rec for rec in recs}

def typeCheck(val):
    return (type(val) in [types.FloatType, types.IntType]) 

modelDir = '/home/welling/git/pyrhea/models/ChicagoLand'

allKeySet, recs = yaml_tools.parse_all(os.path.join(modelDir, 'facilityfacts'))
facDict = {r['abbrev']:r for r in recs}

newRecs = []

snfCategoryLOSModel = {'parms': [0.999, 3.19055827, 1.11488425, 0.03982806],
                       'pdf': '$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)',
                       'prov': 'los_model_fits.py version fbb96209 category average using this LOS model'
                       }
# patch1Facs = ['PETE_520_S', 'LITT_2325_S', 'LITT_80_S', 'LUTH_16220_S']
patch1Facs = []
for abbrev in patch1Facs:
    assert abbrev not in [r['abbrev'] for r in newRecs], '%s has conflicting patches!' % abbrev
    nR = facDict[abbrev].copy()
    nR['losModel'] = snfCategoryLOSModel
    newRecs.append(nR)

patch2Facs = ['ELMH_200_H', 'NORT_660_H', 'BETH_5025_H', 'MERC_901_H',
              'THE_600_H', 'MUNS_901_H', 'FRAN_1201_H', 'FRAN_5454_H',
              'FRAN_24_H', 'FRAN_701_H', 'AURO_10400_H', 'INDI_1007_H',
              'UNIT_6308_H', 'FRAN_301_H', 'JASP_1104_H', 'STC_4321_H',
              'PORT_85_H', 'STM_1500_H', 'PINN_9301_H']

fapdFitSlope = 0.90
fapdFitIntercept = 0.0
fapdProvStr = ('fit line with intercept=0.0 to value vs. (nBedsICU/nBeds) for known HOSPITALS gives slope=%s' %
               fapdFitSlope)
      
medianMeanLOSICU = 3.3
meanLOSICUProvStr = 'median of HOSPITALs with known values'

# This table actually uses the exact estimated slope, which is 0.895575586139
patch2Tbl = """abbrev, estFapdVal, nBeds, bedRatio,
AURO_10400_H, 0.193638, 74.0, 0.216216216216,
BETH_5025_H, 0.620014, 13.0, 0.692307692308,
ELMH_200_H, 0.824872, 38.0, 0.921052631579,
FRAN_1201_H, 0.083486, 236.0, 0.0932203389831,
FRAN_24_H, 0.077395, 162.0, 0.0864197530864,
FRAN_301_H, 0.068514, 183.0, 0.0765027322404,
FRAN_5454_H, 0.065610, 273.0, 0.0732600732601,
FRAN_701_H, 0.146549, 55.0, 0.163636363636,
INDI_1007_H, 0.089558, 200.0, 0.1,
JASP_1104_H, 0.154410, 29.0, 0.172413793103,
MERC_901_H, 0.079021, 34.0, 0.0882352941176,
MUNS_901_H, 0.031286, 458.0, 0.0349344978166,
NORT_660_H, 0.160402, 67.0, 0.179104477612,
PINN_9301_H, 0.099508, 18.0, 0.111111111111,
PORT_85_H, 0.092235, 301.0, 0.102990033223,
STC_4321_H, 0.049411, 290.0, 0.0551724137931,
STM_1500_H, 0.090008, 199.0, 0.100502512563,
THE_600_H, 0.041334, 260.0, 0.0461538461538,
UNIT_6308_H, 0.110170, 252.0, 0.123015873016,
"""

sio = StringIO(patch2Tbl)
sio.name = 'NotARealFile'
keys, recs = csv_tools.parseCSV(sio)
p2Dict = {r['abbrev']: r for r in recs}

for abbrev in patch2Facs:
    assert abbrev not in [r['abbrev'] for r in newRecs], '%s has conflicting patches!' % abbrev
    nR = facDict[abbrev].copy()
    nR['meanLOSICU'] = {'value': medianMeanLOSICU,
                        'prov': meanLOSICUProvStr}
    nR['fracAdultPatientDaysICU'] = {'value': p2Dict[abbrev]['estFapdVal'],
                                     'prov': fapdProvStr}
    newRecs.append(nR)
                
print '%d of %d records modified' % (len(newRecs), len(facDict))
yaml_tools.save_all(os.path.join(modelDir, 'facilityfactsUpdated'), newRecs)

    
    

