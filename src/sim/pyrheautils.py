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
from facilitybase import FacilityRegistry

import schemautils

logger = logging.getLogger(__name__)

PATH_STRING_MAP = {}

def pathTranslate(rawPath, lookupDict=None):
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
    return path
    

def syncXDRORegistries(facilities):
    import time
    syncRegistry = FacilityRegistry("sync reg") 
    time1 = time.time()
    for fac in facilities:
        if fac.registry.hasRegistry("xdroRegistry"):
            syncRegistry.transferFromOtherRegistry(fac.registry,"xdroRegistry","sync reg")
    time2 = time.time()
    print "time for to sync = {0}".format(time2-time1)
    time1 - time.time()
    for fac in facilities:
        if fac.registry.hasRegistry("xdroRegistry"):
            fac.registry.transferFromOtherRegistry(syncRegistry,"sync reg","xdroRegistry")
    time2 = time.time()
    print 'time for from syn = {0}'.format(time2-time1)
            
def importConstants(valuePath, schemaPath, pathLookupDict=None):
    with open(pathTranslate(valuePath, pathLookupDict), 'rU') as f:
        cJSON = yaml.safe_load(f)
        if os.name != 'nt':
            validator = schemautils.getValidator(pathTranslate(schemaPath, pathLookupDict))
            try:
                nErrors = sum([1 for e in validator.iter_errors(cJSON)])  # @UnusedVariable
                if nErrors:
                    for e in validator.iter_errors(cJSON):
                        logger.error('Schema violation: %s: %s' %
                                     (' '.join([str(word) for word in e.path]), e.message))
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
