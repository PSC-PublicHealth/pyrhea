#! /usr/bin/env python

import os.path
import sys
import yaml
import tools_util as tu
import pyrheautils as pu
import csv

def getCapturedLocs(topDir, capture):
    fn = os.path.join(topDir, 'facility_capture_percentages.csv')
    key = '{0}% of facilities'.format(capture)
    locL = []
    with open(fn, 'rU') as f:
        reader = csv.DictReader(f)
        for l in reader:
            if l[key] is not None and l[key] != '':
                locL.append(l[key])
    return locL


def getDiagMods(search, add):
    modList = []
    modList.append([['registryAddCompliance', 0, 'frac', 'value'], add / 100.0])
    modList.append([['registryAddCompliance', 1, 'frac', 'value'], add / 100.0])
    modList.append([['registryAddCompliance', 2, 'frac', 'value'], add / 100.0])
    modList.append([['registryAddCompliance', 3, 'frac', 'value'], add / 100.0])
    modList.append([['registrySearchCompliance', 0, 'frac', 'value'], search / 100.0])
    modList.append([['registrySearchCompliance', 1, 'frac', 'value'], search / 100.0])
    modList.append([['registrySearchCompliance', 2, 'frac', 'value'], search / 100.0])
    modList.append([['registrySearchCompliance', 3, 'frac', 'value'], search / 100.0])

    return modList


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    runYaml = argv[0]
    topDir = argv[1]
    bundleDir = argv[2]
    pyrheaDir = argv[3]
    workDir = argv[4]
    scenario = argv[5]

    fn = os.path.join(bundleDir, runYaml)
    inputDict = tu.readModelInputs(fn)
    pu.prepPathTranslations(inputDict)

    capture = int(scenario.split('capture_')[1].split('_')[0])
    search = int(scenario.split('search_')[1].split('_')[0])
    add = int(scenario.split('add_')[1].split('_')[0])
    if 'eff_' in scenario:
        enhancedDetectionFraction = int(scenario.split('eff_')[1].split('_')[0])

    else:
        enhancedDetectionFraction = None
    
    captureLocL = getCapturedLocs(topDir, capture)

    xdroModList = [[['locationsImplementingScenario', 'locAbbrevList'], captureLocL]]
    if enhancedDetectionFraction is not None:
        xdroModList.append([['enhancedDetectionFraction', 'value'], enhancedDetectionFraction/100.0])
    xdroPath = pu.pathTranslate('$(CONSTANTS)/xdro_registry_scenario_constants.yaml')
    diagModList = getDiagMods(search, add)
    diagPath = pu.pathTranslate('$(CONSTANTS)/$(PATHOGEN)/generic_diagnosis_constants.yaml')

    allData = { 'constantsReplacementData': {xdroPath: xdroModList,
                                             diagPath: diagModList},
                'facilitiesReplacementData': {}
                }

    yaml.dump(allData, sys.stdout,
              indent=4, encoding='utf-8',width=130,explicit_start=True)

if __name__ == "__main__":
    main()
