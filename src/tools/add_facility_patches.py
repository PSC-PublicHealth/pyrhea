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

# patch2Facs = ['ELMH_200_H', 'NORT_660_H', 'BETH_5025_H', 'MERC_901_H',
#               'THE_600_H', 'MUNS_901_H', 'FRAN_1201_H', 'FRAN_5454_H',
#               'FRAN_24_H', 'FRAN_701_H', 'AURO_10400_H', 'INDI_1007_H',
#               'UNIT_6308_H', 'FRAN_301_H', 'JASP_1104_H', 'STC_4321_H',
#               'PORT_85_H', 'STM_1500_H', 'PINN_9301_H']
patch2Facs = []

fapdFitSlope = 0.90
fapdFitIntercept = 0.0
fapdProvStr = ('fit line with intercept=0.0 to value vs. (nBedsICU/nBeds) for known HOSPITALS gives slope=%s' %
               fapdFitSlope)
      
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

nRList = applyPatch(patch2Facs, patch2Tbl, 'fracAdultPatientDaysICU', fapdProvStr,
                    tblKey='estFapdVal')
if noOverlap(nRList, newRecs):
    newRecs += nRList

# patch2bFacs = ['ELMH_200_H', 'NORT_660_H', 'BETH_5025_H', 'MERC_901_H',
#                'THE_600_H', 'MUNS_901_H', 'FRAN_1201_H', 'FRAN_5454_H',
#                'FRAN_24_H', 'FRAN_701_H', 'AURO_10400_H', 'INDI_1007_H',
#                'UNIT_6308_H', 'FRAN_301_H', 'JASP_1104_H', 'STC_4321_H',
#                'PORT_85_H', 'STM_1500_H', 'PINN_9301_H']
patch2bFacs = []
medianMeanLOSICU = 3.3
meanLOSICUProvStr = 'median of HOSPITALs with known values'

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
AURO_10400_H, 44.876450,
FRAN_1201_H, 165.849467,
FRAN_24_H, 110.590187,
FRAN_301_H, 126.271875,
FRAN_5454_H, 193.479106,
FRAN_701_H, 30.688257,
INDI_1007_H, 138.966574,
JASP_1104_H, 11.272834,
MUNS_901_H, 331.627304,
PINN_9301_H, 3.058617,
PORT_85_H, 214.388023,
STC_4321_H, 206.173806,
STM_1500_H, 138.219827,
THE_600_H, 183.771395,
UNIT_6308_H, 177.797419,
"""

meanPopProvStr = ("Linear regression over hospitals with known meanPop and nBeds"
                  "yields slope: 0.746747016753, intercept: -10.3828292868, "
                  "std_err: 0.0244206131928, r_value: 0.963111004474, "
                  "p_value: 2.33822158602e-43")
nRList = applyPatch(patch3Facs, patch3Tbl, 'meanPop', meanPopProvStr,
                    tblKey='estMeanPop')
if noOverlap(nRList, newRecs):
    newRecs += nRList

# patch4Facs = ['HANC_203_S', 'JOHN_5909_S', 'SSC_10330_S', 'RENS_1309_S', 'MAJO_2350_S',
#               'MAJO_8380_S', 'MICH_802_S', 'THE_110_S', 'THE_710_S', 'KIND_2300_S',
#               'MAJO_7935_S', 'MIDW_1519_S', 'HEND_251_S', 'HEND_1700_S', 'MAJO_10352_S',
#               'CONS_1000_S', 'HANC_101_S', 'GEOR_3623_S', 'WOOD_7444_S', 'HEAR_3100_S',
#               'ADAM_119_S', 'JOHN_2901_S', 'TRIL_1101_S', 'MAJO_5025_S', 'EXTE_8633_S',
#               'THE_606_S', 'BROO_3506_S', 'ST_9244_S', 'HEND_8800_S', 'WITT_1200_S',
#               'CHIC_6685_S', 'KIND_3415_S', 'HEND_3175_S', 'KIND_8400_S', 'OAK_221_S',
#               'DECA_353_S', 'HEND_1900_S', 'MAJO_3301_S', 'FAL_9630_S', 'VALP_3405_S',
#               'LIFE_1000_S', 'MAJO_4410_S', 'MAJO_601_S', 'COMM_503_S', 'BEVE_1703_S',
#               'MIDW_8540_S', 'MCAL_18300_S']
patch4Facs = []
patch4Tbl = """abbrev, meanPop, note
ADAM_119_S, 48.916559,
BEVE_1703_S, 86.152309,
BROO_3506_S, 137.919082,
CHIC_6685_S, 131.561759,
COMM_503_S, 100.683333,
CONS_1000_S, 61.631205,
DECA_353_S, 63.447583,
EXTE_8633_S, 88.876876,
FAL_9630_S, 33.477346,
GEOR_3623_S, 62.539394,
HANC_101_S, 47.100181,
HANC_203_S, 167.889320,
HEAR_3100_S, 107.948845,
HEND_1700_S, 77.978607,
HEND_1900_S, 158.807430,
HEND_251_S, 76.162229,
HEND_3175_S, 167.889320,
HEND_8800_S, 147.909161,
JOHN_2901_S, 98.866955,
JOHN_5909_S, 58.906638,
KIND_2300_S, 147.909161,
KIND_3415_S, 77.070418,
KIND_8400_S, 69.804906,
LIFE_1000_S, 66.172150,
MAJO_10352_S, 83.427742,
MAJO_2350_S, 162.440186,
MAJO_3301_S, 135.194515,
MAJO_4410_S, 124.296247,
MAJO_5025_S, 105.224278,
MAJO_601_S, 145.184594,
MAJO_7935_S, 203.308691,
MAJO_8380_S, 89.785065,
MCAL_18300_S, 87.968687,
MICH_802_S, 107.948845,
MIDW_1519_S, 77.978607,
MIDW_8540_S, 302.301294,
OAK_221_S, 52.549315,
RENS_1309_S, 141.551838,
SSC_10330_S, 106.132467,
ST_9244_S, 80.703174,
THE_110_S, 89.785065,
THE_606_S, 147.909161,
THE_710_S, 77.070418,
TRIL_1101_S, 162.440186,
VALP_3405_S, 78.886796,
WITT_1200_S, 139.735460,
WOOD_7444_S,,no nBeds
"""
patch4ProvStr=("Linear regression over SNFs and VSNFs with known meanPop and nBeds"
               " yields slope: 0.908189012909, intercept: -1.03383666603, "
                  "std_err: 0.0178734280511, r_value: 0.944229978962, "
                  "p_value: 1.57105127932e-153")

nRList = applyPatch(patch4Facs, patch4Tbl, 'meanPop', patch4ProvStr)
if noOverlap(nRList, newRecs):
    newRecs += nRList

patch5Facs = ['VIBR_9509_L']
patch5Tbl = """abbrev, meanPop, note
VIBR_9509_L, 19.448458,
"""
patch5ProvStr=("Linear regression over LTACHs with known meanPop and nBeds"
               " yields slope: 0.915393514265, intercept: -17.1672830436, "
                  "std_err: 0.0896885497054, r_value: 0.972388088035, "
                  "p_value: 5.1545642312e-05")

nRList = applyPatch(patch5Facs, patch5Tbl, 'meanPop', patch5ProvStr)
if noOverlap(nRList, newRecs):
    newRecs += nRList

print '%d of %d records modified' % (len(newRecs), len(facDict))
yaml_tools.save_all(os.path.join(modelDir, 'facilityfactsUpdated'), newRecs)

    
    

