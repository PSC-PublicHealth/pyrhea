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

    fn = pu.pathTranslate('$(CONSTANTS)/xdro_registry_scenario_constants.yaml')
    if not os.path.isabs(fn):
        fn = os.path.join(os.path.dirname(__file__), fn)
    with open(fn, 'rU') as f:
        xdroJSON = yaml.load(f)

    sCD = pu.pathTranslate('$(MODELDIR)/scenario_constants')
    if 'capture_' in scenario:
        # this is XDRO-style
        capture = int(scenario.split('capture_')[1].split('_')[0])
        key = '{0}% of facilities'.format(capture)
        scenariosCSV = os.path.join(sCD,
                                    'facility_capture_percentages.csv')
    else:
        # this is CRE-style
        if scenario.startswith('baseline_'):
            key = scenario[len('baseline_'):]
        else:
            key = scenario
        scenariosCSV = os.path.join(sCD, 'scenario_groups.csv')

    locL = getCapturedLocs(scenariosCSV, key)

    xdroJSON['locationsImplementingScenario'] = {
        'prov': 'scenario %s' % scenario,
        'locAbbrevList': locL
    }

    yaml.dump(xdroJSON, sys.stdout,
              indent=4, encoding='utf-8',width=130,explicit_start=True)


if __name__ == "__main__":
    main()
