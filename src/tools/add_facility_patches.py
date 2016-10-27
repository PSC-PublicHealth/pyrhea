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

def applyPatch(facList, patchTblStr, key, provStr, tblKey=None):
    """
    Returns a list of patched records
    """
    nRList = []
    if tblKey is None:
        tblKey = key
    sio = StringIO(patchTblStr)
    sio.name = 'NotARealFile'
    keys, recs = csv_tools.parseCSV(sio)
    assert tblKey in keys, ("Cannot patch %s with %s because it is not a table column" %
                            (key, tblKey))
    pDict = {r['abbrev']: r for r in recs}
    
    for abbrev in facList:
        if not typeCheck(pDict[abbrev][tblKey]):
            print 'Skipped %s <%s>' % (abbrev, pDict[abbrev][tblKey])
        nR = facDict[abbrev].copy()
        nR[key] = {'value': pDict[abbrev][tblKey], 'prov': provStr}
        nRList.append(nR)
        
    return nRList
    
def noOverlap(rL1, rL2):
    """
    Raises AssertionError if there is an overlap, returns True otherwise
    """
    for abbrev in [rec['abbrev'] for rec in rL1]:
        assert abbrev not in [rec['abbrev'] for rec in rL2], '%s has conflicting patches!' % abbrev
    return True

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

# patch2Facs = ['THE_600_H', 'MUNS_901_H', 'FRAN_1201_H', 'FRAN_5454_H',
#               'FRAN_24_H', 'FRAN_701_H', 'AURO_10400_H', 'INDI_1007_H',
#               'UNIT_6308_H', 'FRAN_301_H', 'JASP_1104_H', 'STC_4321_H',
#               'PORT_85_H', 'STM_1500_H', 'PINN_9301_H']
patch2Facs = []
fapdFitSlope = 0.8955900058
fapdFitIntercept = 0.0
fapdProvStr = ('fit line with intercept=0.0 to value vs. (nBedsICU/nBeds) for known HOSPITALS gives slope=%s' %
               fapdFitSlope)
      
# This table actually uses the exact estimated slope, which is 0.8955900058
patch2Tbl = """abbrev, estFapdVal, nBeds, bedRatio,
AURO_10400_H, 0.193641, 74.0, 0.216216216216,
FRAN_1201_H, 0.083487, 236.0, 0.0932203389831,
FRAN_24_H, 0.077397, 162.0, 0.0864197530864,
FRAN_301_H, 0.068515, 183.0, 0.0765027322404,
FRAN_5454_H, 0.065611, 273.0, 0.0732600732601,
FRAN_701_H, 0.146551, 55.0, 0.163636363636,
INDI_1007_H, 0.089559, 200.0, 0.1,
JASP_1104_H, 0.154412, 29.0, 0.172413793103,
MUNS_901_H, 0.031287, 458.0, 0.0349344978166,
PINN_9301_H, 0.099510, 18.0, 0.111111111111,
PORT_85_H, 0.092237, 301.0, 0.102990033223,
STC_4321_H, 0.049412, 290.0, 0.0551724137931,
STM_1500_H, 0.090009, 199.0, 0.100502512563,
THE_600_H, 0.041335, 260.0, 0.0461538461538,
UNIT_6308_H, 0.110172, 252.0, 0.123015873016,
"""

nRList = applyPatch(patch2Facs, patch2Tbl, 'fracAdultPatientDaysICU', fapdProvStr,
                    tblKey='estFapdVal')
if noOverlap(nRList, newRecs):
    newRecs += nRList

# patch2bFacs = ['AURO_10400_H', 'FRAN_1201_H', 'FRAN_24_H', 'FRAN_301_H',
#                'FRAN_5454_H', 'FRAN_701_H', 'INDI_1007_H', 'JASP_1104_H',
#                'MUNS_901_H', 'PINN_9301_H', 'PORT_85_H', 'STC_4321_H',
#                'STM_1500_H', 'THE_600_H', 'UNIT_6308_H']
patch2bFacs = []
medianMeanLOSICU = 3.27
meanLOSICUProvStr = 'median of HOSPITALs with known values from Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS:$W 14ba183e'

for abbrev in patch2bFacs:
    assert abbrev not in [r['abbrev'] for r in newRecs], '%s has conflicting patches!' % abbrev
    nR = facDict[abbrev].copy()
    nR['meanLOSICU'] = {'value': medianMeanLOSICU, 'prov': meanLOSICUProvStr}
    newRecs.append(nR)



# patch3Facs = ['PINN_9301_H', 'FRAN_5454_H', 'MUNS_901_H', 'FRAN_1201_H', 'FRAN_24_H',
#               'THE_600_H', 'FRAN_301_H', 'JASP_1104_H', 'STM_1500_H', 'STC_4321_H',
#               'AURO_10400_H', 'UNIT_6308_H', 'PORT_85_H', 'INDI_1007_H', 'FRAN_701_H']
patch3Facs = []

patch3Tbl = """abbrev, estMeanPop,
AURO_10400_H, 42.682371,
FRAN_1201_H, 165.238324,
FRAN_24_H, 109.255975,
FRAN_301_H, 125.142858,
FRAN_5454_H, 193.229499,
FRAN_701_H, 28.308525,
INDI_1007_H, 138.003668,
JASP_1104_H, 8.639050,
MUNS_901_H, 333.185372,
PINN_9301_H, 0.317350,
PORT_85_H, 214.412010,
STC_4321_H, 206.090309,
STM_1500_H, 137.247150,
THE_600_H, 183.394762,
UNIT_6308_H, 177.342616,
"""

meanPopProvStr = ("Linear regression over hospitals with known meanPop and nBeds "
                  "from Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS 14ba183e "
                  "yields slope: 0.75651823171, intercept: -13.2999782334, "
                  "std_err: 0.0254409192632, r_value: 0.961113716723, "
                  "p_value: 1.54689587629e-42")
nRList = applyPatch(patch3Facs, patch3Tbl, 'meanPop', meanPopProvStr,
                    tblKey='estMeanPop')
if noOverlap(nRList, newRecs):
    newRecs += nRList

# patch4Facs = ['ADAM_119_S', 'BEVE_1703_S', 'BROO_3506_S', 'CHIC_6685_S', 'COMM_503_S',
#               'CONS_1000_S', 'DECA_353_S', 'EXTE_8633_S', 'FAL_9630_S', 'GEOR_3623_S',
#               'HANC_101_S', 'HANC_203_S', 'HEAR_3100_S', 'HEND_1700_S', 'HEND_1900_S',
#               'HEND_251_S', 'HEND_3175_S', 'HEND_8800_S', 'JOHN_2901_S', 'JOHN_5909_S',
#               'KIND_2300_S', 'KIND_3415_S', 'KIND_8400_S', 'LIFE_1000_S', 'MAJO_10352_S',
#               'MAJO_2350_S', 'MAJO_3301_S', 'MAJO_4410_S', 'MAJO_5025_S', 'MAJO_601_S',
#               'MAJO_7935_S', 'MAJO_8380_S', 'MCAL_18300_S', 'MICH_802_S', 'MIDW_1519_S',
#               'MIDW_8540_S', 'OAK_221_S', 'RENS_1309_S', 'RESU_2380_S', 'SSC_10330_S',
#               'ST_9244_S', 'THE_110_S', 'THE_606_S', 'THE_710_S', 'TRIL_1101_S', 'VALP_3405_S',
#               'WITT_1200_S', 'WOOD_7444_S']
patch4Facs = []
patch4Tbl = """abbrev, meanPop, note
ADAM_119_S, 47.822603,
BEVE_1703_S, 85.482558,
BROO_3506_S, 137.839080,
CHIC_6685_S, 131.409332,
COMM_503_S, 100.179126,
CONS_1000_S, 60.682100,
DECA_353_S, 62.519171,
EXTE_8633_S, 88.238164,
FAL_9630_S, 32.207500,
GEOR_3623_S, 61.600636,
HANC_101_S, 45.985532,
HANC_203_S, 168.150751,
HEAR_3100_S, 107.527409,
HEND_1700_S, 77.215739,
HEND_1900_S, 158.965396,
HEND_251_S, 75.378668,
HEND_3175_S, 168.150751,
HEND_8800_S, 147.942971,
JOHN_2901_S, 98.342055,
JOHN_5909_S, 57.926494,
KIND_2300_S, 147.942971,
KIND_3415_S, 76.297203,
KIND_8400_S, 68.948919,
LIFE_1000_S, 65.274777,
MAJO_10352_S, 82.726952,
MAJO_2350_S, 162.639538,
MAJO_3301_S, 135.083474,
MAJO_4410_S, 124.061048,
MAJO_5025_S, 104.771803,
MAJO_601_S, 145.187364,
MAJO_7935_S, 203.973635,
MAJO_8380_S, 89.156700,
MCAL_18300_S, 87.319629,
MICH_802_S, 107.527409,
MIDW_1519_S, 77.215739,
MIDW_8540_S, 304.094002,
OAK_221_S, 51.496745,
RENS_1309_S, 141.513222,
RESU_2380_S,,no nBeds
SSC_10330_S, 105.690339,
ST_9244_S, 79.971345,
THE_110_S, 89.156700,
THE_606_S, 147.942971,
THE_710_S, 76.297203,
TRIL_1101_S, 162.639538,
VALP_3405_S, 78.134274,
WITT_1200_S, 139.676151,
WOOD_7444_S,,no nBeds
"""
patch4ProvStr=("Linear regression over SNFs with known meanPop and nBeds"
               " if nBeds >= meanPop "
               " from Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS 14ba183e"
               " yields slope: 0.91853547892, intercept: -2.69684797399, "
                  "std_err: 0.0191476061526, r_value: 0.940750150331, "
                  "p_value: 1.81252301046e-142")

nRList = applyPatch(patch4Facs, patch4Tbl, 'meanPop', patch4ProvStr)
if noOverlap(nRList, newRecs):
    newRecs += nRList

patch5Facs = ['VIBR_9509_L']
patch5Tbl = """abbrev, meanPop, note
VIBR_9509_L, 18.244795,
"""
patch5ProvStr=("Linear regression over LTACHs with known meanPop and nBeds"
               " from Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016.xlsx:BedLOS 14ba183e"
               " yields slope: 0.963881278536, intercept: -20.3104566208, "
                  "std_err: 0.122932980501, r_value: 0.954505298723, "
                  "p_value: 0.00022744937517")

nRList = applyPatch(patch5Facs, patch5Tbl, 'meanPop', patch5ProvStr)
if noOverlap(nRList, newRecs):
    newRecs += nRList

print '%d of %d records modified' % (len(newRecs), len(facDict))
yaml_tools.save_all(os.path.join(modelDir, 'facilityfactsUpdated'), newRecs)

    
    

