#! /usr/bin/env python

"""
This recreates a set of facilityfacts YAML files with LOS models added.  The
sources of the data are csv files; different sources are used by facility type.
"""

import os.path
import yaml
import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools

def loadCSVByAbbrev(modelDir, fname):
    fullName = os.path.join(modelDir, fname)
    with open(fullName) as fl:
        keys, recs = csv_tools.parseCSV(fl)
    assert 'abbrev' in keys, ('%s has no "abbrev" field' % fullName)
    return {rec['abbrev']: rec for rec in recs}

modelDir = '/home/welling/git/pyrhea/models/ChicagoLand'
bimodalFitsFName = 'los_model_fit_parms_bimodal.csv'
lognormFitsFName = 'los_model_fit_parms_lognorm.csv'
twopopFitsFName = 'los_model_fit_parms_twopop.csv'

allKeySet, recs = yaml_tools.parse_all(os.path.join(modelDir, 'facilityfacts'))
facDict = {r['abbrev']:r for r in recs}

bimodalFits = loadCSVByAbbrev(modelDir, bimodalFitsFName)
lognormFits = loadCSVByAbbrev(modelDir, lognormFitsFName)
twopopFits = loadCSVByAbbrev(modelDir, twopopFitsFName)

losModelProvStr = 'los_model_fits.py version fbb96209 using the given LOS model'

newRecs = []
for abbrev, rec in facDict.items():
    modelToUse, modelNm = {'LTACH': (bimodalFits, 'bimodal'),
                           'HOSPITAL': (lognormFits, 'lognorm'),
                           'SNF': (twopopFits, 'twopop'),
                           'VSNF': (twopopFits, 'twopop')}[rec['category']]
    if abbrev in modelToUse:                       
        losR = modelToUse[abbrev]
        lmD = {'prov': losModelProvStr}
        if modelNm == 'bimodal':
            lmD['pdf'] = "$0*lognorm(mu=$1,sigma=$2)+(1-$0)*lognorm(mu=$3, sigma=$4)"
            lmD['parms'] = [losR[nm] for nm in ['k', 'mu', 'sigma', 'mu2', 'sigma2']]
        elif modelNm == 'lognorm':
            lmD['pdf'] = "lognorm(mu=$0,sigma=$1)"
            lmD['parms'] = [losR[nm] for nm in ['mu', 'sigma']]
        elif modelNm == 'twopop':
            lmD['pdf'] = "$0*lognorm(mu=$1,sigma=$2)+(1-$0)*expon(lambda=$3)"
            lmD['parms'] = [losR[nm] for nm in ['k', 'mu', 'sigma', 'lmda']]
        else:
            raise RuntimeError('Unexpected los model name %s' % modelNm)
        nR = rec.copy()
        nR['losModel'] = lmD
        newRecs.append(nR)
    else:
        print 'No %s data for %s' % (modelNm, abbrev)

print '%d of %d records modified' % (len(newRecs), len(recs))
yaml_tools.save_all(os.path.join(modelDir, 'facilityfactsUpdated'), newRecs)

    
    

