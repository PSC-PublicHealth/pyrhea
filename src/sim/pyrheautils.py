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

import yaml
import logging
import sys
import os
from imp import load_source

import schemautils

logger = logging.getLogger(__name__)


def importConstants(valuePath, schemaPath):
    with open(valuePath, 'rU') as f:
        cJSON = yaml.safe_load(f)
    validator = schemautils.getValidator(schemaPath)
    nErrors = sum([1 for e in validator.iter_errors(cJSON)])  # @UnusedVariable
    if nErrors:
        for e in validator.iter_errors(cJSON):
            logger.error('Schema violation: %s: %s' %
                         (' '.join([str(word) for word in e.path]), e.message))
        raise RuntimeError('%s does not satisfy the schema %s' %
                           (valuePath, schemaPath))
    return cJSON


def loadModulesFromDir(implementationDir,
                       requiredAttrList=[]):
    """
    Scan the given directory for .py files and import them as modules one by one.
    """
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
