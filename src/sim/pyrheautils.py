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

import schemautils

logger = logging.getLogger(__name__)

PATH_STRING_MAP = {}

def prepPathTranslations(inp):
    """
    load up PATH_STRING_MAP based on the (already read in) pyrhea input yaml file/input dict
    """
    primaryKeys = [('modelDir', 'MODELDIR'),
                   ('facilityImplementationDir', 'IMPLDIR'),
                   ('policyImplementationDir', 'POLICYDIR'),
                   ('pathogenImplementationDir', 'PATHOGENDIR')]


    for inputKey, xlateKey in primaryKeys:
        if inputKey in inp:
            PATH_STRING_MAP[xlateKey] = inp[inputKey]

    if 'pathTranslations' in inp:
        for elt in inputDict['pathTranslations']:
            pyrheautils.PATH_STRING_MAP[elt['key']] = elt['value']


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
            if ('$(%s)' % key) in path:
                path = path.replace('$(%s)' % key, val)
                changed = True

    return os.path.abspath(path)


def importConstants(valuePath, schemaPath, pathLookupDict=None):
    """
    Import a set of constants in YAML format, checking against its schema.
    """
    with open(pathTranslate(valuePath, pathLookupDict), 'rU') as f:
        cJSON = yaml.safe_load(f)
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
            newMod = load_source(name, os.path.join(implementationDir, fname))
            for requiredAttr in requiredAttrList:
                if not hasattr(newMod, requiredAttr):
                    raise RuntimeError('The implementation in '
                                       '%s is missing the required attribute %s' %
                                       (os.path.join(implementationDir, fname),
                                        requiredAttr))
            yield newMod
        elif ext.startswith('.py'):
            pass  # .pyc files for example
        else:
            logger.info('Skipping non-python file %s' % fname)
    sys.path.pop()  # drop implementaionDir
