from qsub_run_utils import RunEnvironment 
import yaml
import sys


def getFacList():
    with open("fac_lists.yaml", 'rU') as f:
        return yaml.load(f)


scenarioDef = [
    [100, 100, 100],
    [63, 63, 63],
    [25, 63, 63],
    [100, 63, 63],
    [63, 25, 63],
    [63, 100, 63],
    [63, 63, 25],
    [63, 63, 100],
    ]

def makeReplPairs(instNum, expNum, facLists):
    (capture, searchComp, addComp) = scenarioDef[expNum]
    listTag = 'fac_%s'%capture
    facList = facLists[listTag]

    rpDict = {}

    rp = []
    rp.append([['locationsImplementingScenario', 'locAbbrevList'], facList])
    rpDict['xdro_registry_scenario_constants.yaml'] = rp

    rp = []
    if 0:
        rp.append([['registryAddCompliance', 'HOSPITAL', 'frac', 'value'], addComp * 0.01])
        rp.append([['registryAddCompliance', 'LTACH', 'frac', 'value'], addComp * 0.01])
        rp.append([['registryAddCompliance', 'VSNF', 'frac', 'value'], addComp * 0.01])
        rp.append([['registryAddCompliance', 'SNF', 'frac', 'value'], addComp * 0.01])

        rp.append([['registrySearchCompliance', 'HOSPITAL', 'frac', 'value'], searchComp * 0.01])
        rp.append([['registrySearchCompliance', 'LTACH', 'frac', 'value'], searchComp * 0.01])
        rp.append([['registrySearchCompliance', 'VSNF', 'frac', 'value'], searchComp * 0.01])
        rp.append([['registrySearchCompliance', 'SNF', 'frac', 'value'], searchComp * 0.01])
    if 1:
        rp.append([['registryAddCompliance', 0, 'frac', 'value'], addComp * 0.01])
        rp.append([['registryAddCompliance', 1, 'frac', 'value'], addComp * 0.01])
        rp.append([['registryAddCompliance', 2, 'frac', 'value'], addComp * 0.01])
        rp.append([['registryAddCompliance', 3, 'frac', 'value'], addComp * 0.01])

        rp.append([['registrySearchCompliance', 0, 'frac', 'value'], searchComp * 0.01])
        rp.append([['registrySearchCompliance', 1, 'frac', 'value'], searchComp * 0.01])
        rp.append([['registrySearchCompliance', 2, 'frac', 'value'], searchComp * 0.01])
        rp.append([['registrySearchCompliance', 3, 'frac', 'value'], searchComp * 0.01])
    rpDict['generic_diagnosis_constants.yaml'] = rp

    return rpDict

def main():
    print "starting main"
    runNum = int(sys.argv[1])

    facList = getFacList()
    fac100 = facList['fac_100']

    expCount = len(scenarioDef)
    
    expNum = runNum % expCount
    instNum = runNum / expCount

    replPairs = makeReplPairs(instNum, expNum, facList)

    if 0:
        rEnv = RunEnvironment('ChicagoLand', 'exp_xdro100_ckpt', expNum, instNum,
                              'multiyear_allfac_ChicagoLand_xdro.yaml', replPairs,
                              pyrheaOpts="-k 730", pyrheaPrefix="dmtcp_launch --new-coordinator ")
    else:
        rEnv = RunEnvironment('ChicagoLand', 'exp_xdro', expNum, instNum,
                              'multiyear_allfac_ChicagoLand_xdro.yaml', replPairs)
    rEnv.buildEverything()
    rEnv.runSim(reallyRun=True)

    
    
    

if __name__ == "__main__":
    main()
