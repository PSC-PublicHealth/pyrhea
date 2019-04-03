#! /usr/bin/env python

import os.path
import sys
import yaml
import tools_util as tu
import pyrheautils as pu
import csv


def getCapturedLocs(fileName, key):
    locL = []
    with open(fileName, 'rU') as f:
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

    sCD = pu.pathTranslate('$(MODELDIR)/scenario_constants')
    if 'capture_' in scenario:
        # this is XDRO-style
        mode = 'xdro'
        capture = int(scenario.split('capture_')[1].split('_')[0])
        search = int(scenario.split('search_')[1].split('_')[0])
        add = int(scenario.split('add_')[1].split('_')[0])
        if 'eff_' in scenario:
            enhancedDetectionFraction = int(scenario.split('eff_')[1].split('_')[0])
        else:
            enhancedDetectionFraction = None
        key = '{0}% of facilities'.format(capture)
        scenariosCSV = os.path.join(sCD,
                                    'facility_capture_percentages.csv')
    else:
        # this is CRE-style
        mode = 'xdropluscrebundle'
        enhancedDetectionFraction = None
        if scenario.startswith('baseline_'):
            key = scenario[len('baseline_'):]
        else:
            key = scenario
        scenariosCSV = os.path.join(sCD, 'scenario_groups.csv')

    captureLocL = getCapturedLocs(scenariosCSV, key)

    allData = {}
    if mode == 'xdro':
        xdroModList = [[['locationsImplementingScenario', 'locAbbrevList'], captureLocL]]
        if enhancedDetectionFraction is not None:
            xdroModList.append([['enhancedDetectionFraction', 'value'],
                                enhancedDetectionFraction/100.0])
        xdroPath = pu.pathTranslate('$(CONSTANTS)/xdro_registry_scenario_constants.yaml')
        diagModList = getDiagMods(search, add)
        diagPath = pu.pathTranslate('$(CONSTANTS)/$(PATHOGEN)/generic_diagnosis_constants.yaml')

        allData['constantsReplacementData'] = {xdroPath: xdroModList,
                                               diagPath: diagModList
                                               }
    elif mode == 'xdropluscrebundle':
        transModList = []
        facEltL = []
        for fac in captureLocL:
            facEltL.append({'abbrev': fac,
                            'times': {'startDate': 1, 'endDate': 50000}})
        transModList.append([['locationsImplementingScenario', 'facilities'], facEltL])
        transModList.append([['locationsImplementingScenario', 'prov'],
                             '%s scenario' % scenario])
        allData['constantsReplacementData'] = {"$(CONSTANTS)/xdro_plus_cre_bundle_scenario_constants.yaml":
                                               transModList,
                                               "$(CONSTANTS)/xdro_registry_scenario_constants.yaml": # strictly for gather_counts.py
                                               [[['locationsImplementingScenario', 'locAbbrevList'], captureLocL]]}
    else:
        raise RuntimeError('Unknown mode %s' % mode)

    allData['facilitiesReplacementData'] = {}

    yaml.dump(allData, sys.stdout,
              indent=4, encoding='utf-8',width=130,explicit_start=True)

if __name__ == "__main__":
    main()
