#! /usr/bin/env python

###################################################################################
# Copyright   2015, Pittsburgh Supercomputing Center (PSC).  All Rights Reserved. #
# =============================================================================== #
#                                                                                 #
# Permission to use, copy, and modify this software and its documentation without #
# fee for personal use within your organization is hereby granted, provided that  #
# the above copyright notice is preserved in all copies and that the copyright    #
# and this permission notice appear in supporting documentation.  All other       #
# restrictions and obligations are defined in the GNU Affero General Public       #
# License v3 (AGPL-3.0) located at http://www.gnu.org/licenses/agpl-3.0.html  A   #
# copy of the license is also provided in the top level of the source directory,  #
# in the file LICENSE.txt.                                                        #
#                                                                                 #
###################################################################################

import sys
import os
import os.path
import logging
from imp import load_source
import yaml
import phacsl.utils.formats.yaml_tools as yaml_tools
from collections import defaultdict

import schemautils

logger = logging.getLogger(__name__)

PATH_STRING_MAP = {}

def prepPathTranslations(inp):
    """
    load up PATH_STRING_MAP based on the (already read in) pyrhea input yaml file/input dict
    """
    global PATH_STRING_MAP

    simDir = os.path.dirname(os.path.abspath(__file__))
    baseDir = os.path.normpath(os.path.join(simDir, os.path.pardir, os.path.pardir))


    primaryKeys = [('modelDir', 'MODELDIR'),
                   ('model', 'MODEL'),
                   ('facilityImplementationDir', 'IMPLDIR'),
                   ('policyImplementationDir', 'POLICYDIR'),
                   ('pathogenImplementationDir', 'PATHOGENDIR'),
                   ('pathogen', 'PATHOGEN'),]

    defaultTranslations = [('BASEDIR', baseDir),
                           ('SIMDIR', simDir),
                           ('MODELDIR', '$(BASEDIR)/models/$(MODEL)'),
                           ('CONSTANTS', '$(MODELDIR)/constants'),
                           ('COMMUNITYCACHEDIR', '$(SIMDIR)/cache/$(MODEL)'),
                           ('AGENTDIR', '$(SIMDIR)/agents/$(MODEL)'),
                           ]

    for k,v in defaultTranslations:
        PATH_STRING_MAP[k] = v

    for inputKey, xlateKey in primaryKeys:
        if inputKey in inp:
            PATH_STRING_MAP[xlateKey] = inp[inputKey]

    if 'pathTranslations' in inp:
        for elt in inp['pathTranslations']:
            PATH_STRING_MAP[elt['key']] = elt['value']


def pathTranslate(rawPath, lookupDict=None):
    """
    Do recursive string substitution on a path based on a dict.
    """
    if lookupDict is None:
        lookupDict = PATH_STRING_MAP
    path = rawPath
    changed = True
    while changed:
        changed = False
        for key, val in lookupDict.items():
            if '$(%s)' % key in path:
                path = path.replace('$(%s)' % key, val)
                changed = True

    #return os.path.abspath(path)
    return path

constantsReplacementData = {}
facilitiesReplacementData = {}
outputNotesName = ""
saveNewConstants = None

def readConstantsReplacementFile(fileName):
    global constantsReplacementData
    global facilitiesReplacementData

    fileName = pathTranslate(fileName)
    file = os.path.basename(fileName)
    name, ext = os.path.splitext(file)

    newMod = load_source(name, fileName)

    if hasattr(newMod, "constantsReplacementData"):
        constantsReplacementData = getattr(newMod, "constantsReplacementData")
    if hasattr(newMod, "facilitiesReplacementData"):
        facilitiesReplacementData = getattr(newMod, "facilitiesReplacementData")

def replaceData(fileName, yData):
    global constantsReplacementData
    global saveNewConstants

    if fileName not in constantsReplacementData:
        return yData

    repl = constantsReplacementData[fileName]

    yaml = replaceYamlData(repl, yData)

    if saveNewConstants is not None:
        fName = os.path.join(saveNewConstants, os.path.basename(fileName))
        yaml_tools.save_one(fName, yaml)

    return yaml

        

def replaceYamlData(replacementData, yData):
    for path, val in replacementData:
        pointer = yData
        idx = 0
        while idx < len(path) - 1:
            step = path[idx]
            if step == '#KeyVal':
                keyKey = path[idx+1]
                keyValue = path[idx+2]
                idx += 3
                for d in pointer:
                    if d[keyKey] == keyValue:
                        pointer = d
                        break
                else:
                    raise KeyError(keyKey)
            else:  #normal case
                try:
                    pointer = pointer[step]
                except:
                    pointer[step] = {}
                    pointer = pointer[step]
                idx += 1
        pointer[path[-1]] = val

    return yData

def importConstants(valuePath, schemaPath, pathLookupDict=None):
    """
    Import a set of constants in YAML format, checking against its schema.
    """
    print valuePath
    with open(pathTranslate(valuePath, pathLookupDict), 'rU') as f:
        cJSON = yaml.safe_load(f)
        cJSON = replaceData(valuePath, cJSON)
        if os.name != 'nt':
            validator = schemautils.getValidator(pathTranslate(schemaPath, pathLookupDict))
            try:
                nErrors = 0
                for e in validator.iter_errors(cJSON):
                    logger.error('Schema violation: %s: %s',
                                 ' '.join([str(word) for word in e.path]), e.message)
                    nErrors += 1
                if nErrors:
                    raise RuntimeError('%s does not satisfy the schema %s' %
                                       (valuePath, schemaPath))
            except AttributeError:
                logger.error('An error occurred validating %s against schema %s',
                             pathTranslate(valuePath, pathLookupDict),
                             pathTranslate(schemaPath, pathLookupDict))
                raise RuntimeError(('An error occurred validating %s against schema %s'
                                    % (pathTranslate(valuePath, pathLookupDict),
                                       pathTranslate(schemaPath, pathLookupDict))))
    return cJSON


def loadModulesFromDir(implementationDir,
                       requiredAttrList=None):
    """
    Scan the given directory for .py files and import them as modules one by one.
    """
    if requiredAttrList is None:
        requiredAttrList = []
    sys.path.append(implementationDir)
    for fname in os.listdir(implementationDir):
        name, ext = os.path.splitext(fname)
        if ext == '.py':
            try:
                newMod = load_source(name, os.path.join(implementationDir, fname))
                for requiredAttr in requiredAttrList:
                    if not hasattr(newMod, requiredAttr):
                        raise RuntimeError('The implementation in '
                                           '%s is missing the required attribute %s' %
                                           (os.path.join(implementationDir, fname),
                                            requiredAttr))
                yield newMod
            except IOError, e:
                logger.warning('Unable to load %s - is a config file missing? %s',
                               fname, str(e))
        elif ext.startswith('.py'):
            pass  # .pyc files for example
        else:
            logger.info('Skipping non-python file %s' % fname)
    sys.path.pop()  # drop implementaionDir
