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

import sys
import os.path
import json
import yaml
from jsonschema import validate
import optparse

jsonSchemaFile = './sample_schema_jsonschema.json'
yamlSchemaFile = './facility_schema.yaml'
tfile = './data/GREE.yaml'



def fileToJSON(fname):
    if fname.lower().endswith('.json') or fname.lower().endswith('.jsn'):
        with open(fname, "r") as f:
            tjson = json.load(f)
    else:
        assert fname.lower().endswith('.yaml') or fname.lower().endswith('.yml'), \
            "File type of %s is not understood" % fname
        with open(fname, "r") as f:
            tjson = yaml.load(f)
    return tjson


def main(myargv=None):
    "Provides a few test routines"

    if myargv is None:
        myargv = sys.argv

    parser = optparse.OptionParser(usage="""
    %prog [-R] schemaFile checkFile

    The schema and the file(s) to be checked can be JSON or YAML.
    """)

    parser.add_option("-R", action="store_true",
                      help="If checkFile is a directory, check files recursively")

    opts, args = parser.parse_args()

    if len(args) == 2:
        schema = fileToJSON(args[0])
        if opts.R and os.path.isdir(args[1]):
            for root, dirs, files in os.walk(args[1]):  # @UnusedVariable
                for fName in files:
                    if fName.lower().endswith(('.yaml', '.yml', '.json', '.jsn')):
                        path = os.path.join(root, fName)
                        try:
                            validate(schema, fileToJSON(path))
                            print ('PASS: %s' % path)
                        except Exception, e:  # @UnusedVariable
                            print ('FAIL: %s' % path)
        else:
            validate(schema, fileToJSON(args[1]))
    else:
        parser.print_help()


############
# Main hook
############

if __name__ == "__main__":
    main()
