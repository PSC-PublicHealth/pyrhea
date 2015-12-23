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

_rhea_svn_id_="$Id$"

import os.path
import json
import yaml
import jsonschema
import urlparse
import urllib

import phacsl.utils.formats.yaml_ordered as yaml_ordered
yaml_ordered.install()


_schemaBasePath = None


def fileToJSON(fname):
    if fname.lower().endswith('.json') or fname.lower().endswith('.jsn'):
        with open(fname, "r") as f:
            tjson = json.load(f)
    else:
        assert fname.lower().endswith('.yaml') or fname.lower().endswith('.yml'), \
            "File type of %s is not understood" % fname
        with open(fname, "r") as f:
            tjson = yaml.safe_load(f)
            # print tjson
    return tjson


def fileURIHandler(uri):
    p = urlparse.urlsplit(uri)
    assert p.scheme == 'file', 'This handler only takes file URIs; got %s' % uri
    ext = os.path.splitext(p.path)[1]
    if ext in ['.yaml', '.yml']:
        with open(p.path, "rU") as f:
            result = yaml.safe_load(f)
    elif ext in ['.jsn', '.json']:
        result = json.loads(urllib.urlopen(uri).read().decode("utf-8"))
    else:
        raise RuntimeError('Unrecognized extension in the file URI %s' % uri)
    return result


def setSchemaBasePath(basePath):
    global _schemaBasePath
    _schemaBasePath = os.path.abspath(basePath)


def getValidator(schemaURI):
    p = urlparse.urlsplit(schemaURI)
    if p.scheme == '':
        # blank filename
        if p.path.startswith('/'):
            basePath = os.path.dirname(os.path.abspath(p.path))
            schema = fileToJSON(p.path)
        elif _schemaBasePath is None:
            basePath = os.path.dirname(os.path.abspath(p.path))
            schema = fileToJSON(p.path)
        else:
            basePath = _schemaBasePath
            schema = fileToJSON(os.path.join(basePath, p.path))
        resolver = jsonschema.RefResolver('file://' + basePath + '/', schema,
                                          handlers={'file': fileURIHandler})
    else:
        schema = json.loads(urllib.urlopen(schemaURI).read().decode("utf-8"))
        resolver = jsonschema.RefResolver(schemaURI, None,
                                          handlers={'file': fileURIHandler})
    validator = jsonschema.validators.validator_for(schema)(schema=schema,
                                                            resolver=resolver)
    return validator
