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

"""
This utility is used to compare directories of YAML files, for example to compare an
updated Facilities directory with the original.  Each YAML file must contain a unique
entry 'abbrev' which is used to match files.  (That is, they are not matched by filename!)

Lists are compared in unordered fashion.  For example, two yaml files of the same name
will match even if the list ['A', 'B', 'C'] in one matches against the list
['A', 'C', 'B'] in the other.
"""

import os.path
import signal
import optparse
import types

import phacsl.utils.formats.yaml_tools as yaml_tools


def innerTermCompare(facDict, updatedFacDict, oldTerm, newTerm, preamble=None):
    if preamble is None:
        preamble = ''
    # print '%s %s -> %s' % (preamble, oldTerm, newTerm)
    if isinstance(oldTerm, types.DictType):
        if isinstance(newTerm, types.DictType):
            allKeys = set(newTerm.keys() + oldTerm.keys())
            for key in allKeys:
                if key in oldTerm:
                    if key in newTerm:
                        if newTerm[key] != oldTerm[key]:
                            newPreamble = '%s %s: ' % (preamble, key)
                            innerTermCompare(facDict, updatedFacDict,
                                             oldTerm[key], newTerm[key], preamble=newPreamble)
                    else:
                        print '%s %s: this key was deleted' % (preamble, key)
                elif key in newTerm:
                    print '%s %s: this key was added' % (preamble, key)
                else:
                    raise RuntimeError('Vanishing key %s' % key)
        else:
            print '%s%s -> %s' % (preamble, oldTerm, newTerm)
    elif isinstance(oldTerm, types.ListType):
        if isinstance(newTerm, types.ListType):
            oldL = oldTerm[:]
            oldL.sort()
            newL = newTerm[:]
            newL.sort()
            idx = 0
            while oldL and newL:
                innerTermCompare(facDict, updatedFacDict,
                                 oldL.pop(0), newL.pop(0), preamble='%s %s:'%(preamble, idx))
                idx += 1
            if oldL:
                print '%s: lost %s in update' % (preamble, oldL)
            elif newL:
                print '%s: gained %s in update' % (preamble, newL)
            else:
                pass  # lists have been compared
        else:
            print '%s%s -> %s' % (preamble, oldTerm, newTerm)
    else:
        if oldTerm != newTerm:
            print '%s%s -> %s' % (preamble, oldTerm, newTerm)
    
def termCompare(facDict, updatedFacDict, abbrev, key):
    if abbrev in facDict:
        if abbrev in updatedFacDict:
            if key in facDict[abbrev]:
                if key in updatedFacDict[abbrev]:
                    oldTerm = facDict[abbrev][key]
                    newTerm = updatedFacDict[abbrev][key]
                    if newTerm != oldTerm:
                        innerTermCompare(facDict, updatedFacDict,
                                         oldTerm, newTerm, '%s: %s: ' % (abbrev, key))
                else:
                    pass
            else:
                if key in updatedFacDict[abbrev]:
                    print '%s: %s field was added in update' % (abbrev, key)
                else:
                    pass # no such abbrev key pair exists
        else:
            pass  # This abbrev was not updated
    elif abbrev in updatedFacDict:
        print '%s record was newly created in the update' % abbrev
    else:
        pass  # no such entry exists
        

def main():
    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)

    parser = optparse.OptionParser(usage="""
    %prog yamlDir1 yamlDir2
    """)
#     parser.add_option("-v", "--verbose", action="store_true",
#                       help="verbose output")

    opts, args = parser.parse_args()  # @UnusedVariable

#     if opts.verbose:
#         verbose = True

    if len(args) == 2:
        dir1Path = args[0]
        dir2Path = args[1]
    else:
        parser.error("Two directories containing YAML collections must be specified.")

    allKeySet, recs = yaml_tools.parse_all(dir1Path)  # @UnusedVariable
    facDict = {r['abbrev']:r for r in recs}
    
    kL, uRecs = yaml_tools.parse_all(dir2Path)  # @UnusedVariable
    updatedFacDict= {r['abbrev']:r for r in uRecs}

    for abbrev, rec in facDict.items():
        if abbrev in updatedFacDict:
            uRec = updatedFacDict[abbrev]
            allKeys = set(rec.keys() + uRec.keys())
            for key in allKeys:
                if key in rec:
                    if key in uRec:
                        termCompare(facDict, updatedFacDict, abbrev, key)
                    else:
                        print '%s: no %s in updated rec' % (abbrev, key)
                else:
                    print '%s: %s key was added' % (abbrev, key)
    
    print '----DONE----'


if __name__ == "__main__":
    main()

