#! /usr/bin/env python

import os.path
import sys
import yaml
import tools_util as tu
import pyrheautils
import csv

simDir = '/pylon5/pscstaff/jleonard/pyrhea/pyrhea/src/sim'
outBaseDir = '/pylon5/pscstaff/jleonard/pyrhea/output/cre_bundle'

scenario = sys.argv[2]
fn = sys.argv[1]

if not os.path.isabs(fn):
    fn = os.path.join(outBaseDir, fn)
inputDict = tu.readModelInputs(fn)

pyrheautils.prepPathTranslations(inputDict)
fn = pyrheautils.pathTranslate('$(CONSTANTS)/xdro_registry_scenario_constants.yaml')
if not os.path.isabs(fn):
    fn = os.path.join(simDir, fn)
with open(fn, 'rU') as f:
    xdroJSON = yaml.load(f)

scenariosCSV = os.path.join(outBaseDir, 'scenario_groups.csv')
locL = []
with open(scenariosCSV, 'rU') as f:
    reader = csv.DictReader(f)
    for l in reader:
        if l[scenario] != '':
            locL.append(l[scenario])

xdroJSON['locationsImplementingScenario'] = {
    'prov': 'scenario %s' % scenario,
    'locAbbrevList': locL
    }

yaml.dump(xdroJSON, sys.stdout,
          indent=4, encoding='utf-8',width=130,explicit_start=True)
