#! /usr/bin/env python

import os.path
import sys
import yaml
import tools_util as tu
import pyrheautils as pu
import csv

masterScenarioDir = '/pylon5/pscstaff/jleonard/pyrhea/output/cre_bundle'

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
    scenario = sys.argv[2]
    fn = sys.argv[1]
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

    if 'capture_' in scenario:
        # this is XDRO-style
        scenariosCSV = os.path.join(topDir,
                                    'facility_capture_percentages.csv')
        capture = int(scenario.split('capture_')[1].split('_')[0])
        key = '{0}% of facilities'.format(capture)
    else:
        # this is CRE-style
        scenariosCSV = os.path.join(masterScenarioDir, 'scenario_groups.csv')
        key = scenario

    locL = getCapturedLocs(scenariosCSV, key)

    xdroJSON['locationsImplementingScenario'] = {
        'prov': 'scenario %s' % scenario,
        'locAbbrevList': locL
    }

    yaml.dump(xdroJSON, sys.stdout,
              indent=4, encoding='utf-8',width=130,explicit_start=True)


if __name__ == "__main__":
    main()
