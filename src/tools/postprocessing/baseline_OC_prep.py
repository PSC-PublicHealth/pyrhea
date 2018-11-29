import csv
import os
import sys
import cPickle as pickle
from collections import defaultdict
from shutil import copytree

import pyrheautils
import tauadjuster
from lockserver import Lock

HOME_CLEAR_COLONIZED = 0.75
COLONIZED_TRANSFER_PROB_SCALE = [('HOSP', 'LTAC', 9.6),
                                 ('LTAC', 'VENT', 7.9),
                                 ('LTAC', 'ICU', 2.0),
                                 ('SKILNRS', 'ICU', 1.148110773),
                                 ('VENT', 'ICU', 1.364772223)]
EXPOSURE_CUTOFF = [('HOSPITAL', 'ICU', 20.0),
                   ('HOSPITAL', 'HOSP', 136.0),
                   ('LTACH', 'LTAC', 59.0),
                   ('VSNF', 'VENT', 62.0),
                   ('VSNF', 'SKILNRS', 51.0),
                   ('VSNF', 'NURSING', 76.0),
                   ('SNF', 'NURSING', 126.0),
                   ('COMMUNITY', 'HOME', 10000.0)]

ICU_WARD_SZ = 20
SNF_DEATH_RATE = 0.07
VSNF_DEATH_RATE = 0.3
COMMUNITY_VERY_SICK_RATE = 0.020

print "*******************"
print "*******************"
print "*******************"
print "moving cache to a better disk"
print "setting initial tau values if init_taudict.pkl is present"
print "setting mrsa_constants.yaml homeClearColonizedStatusProb to %s" % HOME_CLEAR_COLONIZED
#print "setting mrsa_constants.yaml colonizedTransferProbScale elements %s" % COLONIZED_TRANSFER_PROB_SCALE
#print 'setting mrsa_constants.yaml exposureCutoff to %s' % EXPOSURE_CUTOFF
#print "setting ICU ward size to %s" % ICU_WARD_SZ
print "setting SNF deathRate to %s" % SNF_DEATH_RATE
#print "setting VSNF deathRate to %s" % VSNF_DEATH_RATE
print "setting communityVerySickRate to %s" % COMMUNITY_VERY_SICK_RATE
print "*******************"
print "*******************"
print "*******************"

DEFAULT_TAU_DICT = {'HOSP': 0.0092, 'ICU': 0.023, 'NURSING': 6.7e-5,
                    'SKILNRS': 1.0e-4, 'VENT': 1.0e-4, 'LTAC': 5.0e-3}

LTAC_TAU_DICT = {
    'THC_225_L': 1.0e-6,
    'THC_4058_L': 0.007518350336520556,
    'RML_5601_L': 1.0e-06,
    'ADVO_3435_L': 0.005974283433319996,
    'VIBR_9509_L': 1.0e-6,
    'THC_2544_L': 0.013116518349767783,
    'THC_365_L': 0.009897885169845785,
    'THC_6130_L': .005943139421954516,
    'PRES_100_L': .0022323397189482853
}

def mvCache():
    cache = pyrheautils.pathTranslate('$(COMMUNITYCACHEDIR)')
    #newCache = os.path.join(os.environ['LOCAL'], os.path.split(cache)[1])
    newCache = os.path.join('/mnt/consus/users/welling/cache',
                            os.path.split(cache)[1])
    #lockName = '%s:%s' % (os.environ['HOSTNAME'], newCache)
    lockName = '%s:%s' % ('consus', newCache)
 #   with Lock(lockName):
 #       if os.path.exists(newCache):
 #           print 'Someone else copied the cache to %s' % newCache
 #       else:
 #           copytree(cache, newCache)
 #           print 'finished copying'
    # And tell pyrhea that we moved this:
    pyrheautils.PATH_STRING_MAP['COMMUNITYCACHEDIR'] = newCache
    #cpCmd = 'cp -R %s %s'%(cache, localDir)
    #print cpCmd
    #os.system(cpCmd)
    # cache = os.path.join(pyrheautils.pathTranslate('$(COMMUNITYCACHEDIR)'),
    #                      '*')
    # #cache = '/pylon5/pscstaff/welling/git/pyrhea/src/sim/cache/ChicagoLand/*'
    # cpCmd = 'cp -R %s %s'%(cache, localDir)
    # print cpCmd
    # os.system(cpCmd)
    # And tell pyrhea that we moved this:
    #pyrheautils.PATH_STRING_MAP['COMMUNITYCACHEDIR'] = localDir
    #print "finished copying"


def catFromTiers(tierL):
    if 'HOSP' in tierL:
        return 'HOSPITAL'
    elif 'LTAC' in tierL:
        return 'LTAC'
    elif 'HOME' in tierL:
        return 'COMMUNITY'
    elif 'NURSING' in tierL:
        return 'NURSINGHOME'
    else:
        raise RuntimeError('Cannot classify the tier list %s' % tierL)

def filterTaus(tau, fac, tier):
    if tau < 1.0e-8:
        return 1.0e-8
    else:
        return float(tau)
    #return float(tau)
    # if tau <= 0.8 and tau >= 1.0e-8:
    #     return float(tau)
    # else:
    #     return DEFAULT_TAU_DICT[tier]
    #if tier == 'LTAC':
    #    return LTAC_TAU_DICT[fac]
    #else:
    #    return float(tau)
    
def maybeSetInitialTaus(constantsD):
    tauPath = os.path.join(os.path.dirname(__file__), 'init_taudict.pkl')
    if os.path.exists(tauPath):
        print '%s found; initializing taus with these values' % tauPath
        with open(tauPath, 'rU') as f:
            tauDict = pickle.load(f)
        tiersByFac = defaultdict(list)
        for fac, tier in tauDict.keys():
            tiersByFac[fac].append(tier)
        facsByCat = defaultdict(list)
        for fac, tierL in tiersByFac.items():
            facsByCat[catFromTiers(tierL)].append(fac)

        deltas = []
        for cat, facL in facsByCat.items():
            facDeltas = []
            for fac in facL:
                facTierDeltas = []
                for tier in tiersByFac[fac]:
                    val = tauDict[(fac, tier)]
                    val = filterTaus(val, fac, tier)
                    facTierDeltas.append({'frac': {'value': val, 'prov': 'continue taumod run'},
                                          'tier': tier})
                facDeltas.append({'abbrev': fac, 'tiers': facTierDeltas})
            deltas.append({'category': cat, 'facilities': facDeltas})
        
        constantsD['$(CONSTANTS)/mrsa_constants.yaml'] = [[['facilityTau'], deltas]]
    
def buildColTransProbDeltas():
    deltas = []
    for tFrom, tTo, val in COLONIZED_TRANSFER_PROB_SCALE:
        deltas.append({'tierFrom': tFrom, 'tierTo': tTo,
                       'scale': {'prov': 'taumod trial', 'value': val}})
    return deltas
    
def buildExposureCutoffDeltas():
    catTierD = defaultdict(dict)
    for cat, tier, val in EXPOSURE_CUTOFF:
        catTierD[cat][tier] = val
    deltas = []
    for cat, tierD in catTierD.items():
        subLst = []
        for tier, val in tierD.items():
            subLst.append({'tier': tier,
                           'scale': {'prov': 'taumod trial', 'value': val}})
        deltas.append({'category': cat,
                       'tiers': subLst})
    return deltas

COL_DISCH_DELAYS = { 'HOSP': 0.0, 'LTAC': 0.0 }

def buildColDischDelayDeltas():
    deltas = []
    for tier, val in COL_DISCH_DELAYS.items():
        deltas.append({'tier': tier, 'value': val, 'prov': 'DEBUGGING ONLY'})
    return deltas

mvCache()

constantsReplacementData = {}
constantsReplacementData = {'$(CONSTANTS)/mrsa_constants.yaml': []}
maybeSetInitialTaus(constantsReplacementData)

transModList = [[['homeClearColonizedStatusProb', 'value'],
                 HOME_CLEAR_COLONIZED],
                [['homeClearColonizedStatusProb', 'prov'], 
                 "stabilize community prevalence"]]

# transModList.append([['colonizedTransferProbScale'], buildColTransProbDeltas()])

# transModList.append([['exposureCutoff'], buildExposureCutoffDeltas()])

#transModList.append([['colonizedDischargeDelayTime'], buildColDischDelayDeltas()])

constantsReplacementData["$(CONSTANTS)/mrsa_constants.yaml"].extend(transModList)

# constantsReplacementData["$(CONSTANTS)/hospital_constants.yaml"] = [
#     [['bedsPerICUWard', 'value'], ICU_WARD_SZ],
#     [['bedsPerICUWard', 'prov'], 'taumod test case']
# ]

constantsReplacementData["$(CONSTANTS)/nursinghome_constants.yaml"] = [
    [['deathRate', 'value'], SNF_DEATH_RATE],
    [['deathRate', 'prov'], 'taumod test case']
]

constantsReplacementData["$(CONSTANTS)/community_constants.yaml"] = [
    [['communityVerySickRate', 'value'], COMMUNITY_VERY_SICK_RATE],
    [['deathRate', 'prov'], 'test case']
]


#print constantsReplacementData

facilitiesReplacementData = {}
