import csv
import os

import sys
sys.path.append("/pylon5/pscstaff/jleonard/pyrhea/pyrhea/src/sim")

import pyrheautils

onn = pyrheautils.outputNotesName
#onn = '/pylon5/pscstaff/jleonard/pyrhea/output/cre_bundle/initial_scenarios/CRE_Prevalence_13_Mile/pyrhea__SNA_CRE_Prev_13_Mile__20.pkl'


scenario = onn.split('__')[1]

print "*******************"
print "*******************"
print "*******************"
print "cre bundle plus xdro"
print ""
print "scenario: %s"%scenario
print "*******************"
print "*******************"
print "*******************"




def getScenarioFacilities(scenario):
    ret = []
    with open('/pylon5/pscstaff/jleonard/pyrhea/output/cre_bundle/scenario_groups.csv') as f:
        reader = csv.DictReader(f)
        for l in reader:
            if l[scenario] == '':
                continue
            ret.append(l[scenario])
            
    return ret

def mvCache():
    cache = pyrheautils.pathTranslate('$(COMMUNITYCACHEDIR)')
    newCache = os.path.join('/pylon5/pscstaff/jleonard/pyrhea/output/201812_Chicago_xdro_bundle/cache', str(os.environ['SLURM_JOBID']))
    print "copying cache from %s to %s"%(cache, newCache)
    cmd = "mkdir %s"%newCache
    print cmd
    os.system(cmd)
    cmd = "cp -dR %s/* %s"%(cache, newCache)
    print cmd
    os.system(cmd)
    print "finished making copy"
    pyrheautils.PATH_STRING_MAP['COMMUNITYCACHEDIR'] = newCache
    return

def bridges_mvCache():
    localDir = os.environ['LOCAL']
    cache = '/pylon5/pscstaff/jleonard/pyrhea/pyrhea/src/sim/cache/ChicagoLand/*'
    cpCmd = 'cp -R %s %s'%(cache, localDir)
    print cpCmd
    os.system(cpCmd)
    print "finished copying"
    # And tell pyrhea that we moved this:
    pyrheautils.PATH_STRING_MAP['COMMUNITYCACHEDIR'] = localDir
    

mvCache()
facList = getScenarioFacilities(scenario)
print "facilities in scenario:"
print facList

transModList = []

facEltL = []
for fac in facList:
    facEltL.append({'abbrev': fac,
                    'times': {'startDate': 1, 'endDate': 50000}})

transModList.append([['locationsImplementingScenario', 'facilities'],
                     facEltL])
transModList.append([['locationsImplementingScenario', 'prov'],
                     '%s scenario'%scenario])

constantsReplacementData = {"$(CONSTANTS)/xdro_plus_cre_bundle_scenario_constants.yaml":
                            transModList,
                            "$(CONSTANTS)/xdro_registry_scenario_constants.yaml": # strictly for gather_counts.py
                            [[['locationsImplementingScenario', 'locAbbrevList'],facList]]}
facilitiesReplacementData = {}
print constantsReplacementData
