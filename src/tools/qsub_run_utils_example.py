from qsub_run_utils import RunEnvironment 
import yaml
import sys


def getExpList():
    with open("../../models/ChicagoLand/facility_scenario_lists.yaml", 'rU') as f:
        return yaml.load(f)


def makeReplPairs(expNum, expList):
    expKeys = expList.keys()
    expKeys.sort()

    rp = []
    rp.append([['locationsImplementingScenario', 'locAbbrevList'], expList[expKeys[expNum]]])

    rpDict = {'xdro_registry_scenario_constants.yaml': rp,
              'cre_bundle_scenario_constants.yaml': rp,
    }
    return rpDict

def main():
    runNum = int(sys.argv[1])

    expList = getExpList()
    expCount = len(expList)
    
    expNum = runNum % expCount
    instNum = int(runNum / expCount)
    replPairs = makeReplPairs(expNum, expList)

    rEnv = RunEnvironment('ChicagoLand', 'exp_cre_scen', expNum, instNum, 'year_allfac_ChicagoLand.yaml', replPairs)
    rEnv.buildEverything()
    rEnv.runSim(reallyRun=False)

    
    
    

if __name__ == "__main__":
    main()
