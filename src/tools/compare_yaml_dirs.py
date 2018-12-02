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

import signal
import optparse
import types
import os
import yaml

import phacsl.utils.formats.yaml_tools as yaml_tools

_VERBOSE = False

def cleanDicts(thing):
    """
    This is used during printing to remove unsightly OrderedDicts from output
    """
    if isinstance(thing, types.ListType):
        return [cleanDicts(elt) for elt in thing]
    elif isinstance(thing, types.DictType):
        return {key: cleanDicts(val) for key, val in thing.items()}
    else:
        return thing

def innerTermCompare(facDict, updatedFacDict, oldTerm, newTerm, preamble=None):
    """
    Compare two terms, recursively.
    """
    if preamble is None:
        preamble = ''
    # print '%s %s -> %s' % (preamble, oldTerm, newTerm)
    if isinstance(oldTerm, types.DictType):
        oldTerm = {key:val for key, val in oldTerm.items()} # avoid OrderedDicts
        if isinstance(newTerm, types.DictType):
            newTerm = {key:val for key, val in newTerm.items()} # avoid OrderedDicts
            allKeys = set(newTerm.keys() + oldTerm.keys())
            for key in allKeys:
                if key in oldTerm:
                    if key in newTerm:
                        if newTerm[key] != oldTerm[key]:
                            newPreamble = '%s <%s>: ' % (preamble, key)
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
                                 oldL.pop(0), newL.pop(0), preamble='%s [%s]:'%(preamble, idx))
                idx += 1
            if oldL:
                print '%s: lost %s in update' % (preamble, cleanDicts(oldL))
            elif newL:
                print '%s: gained %s in update' % (preamble, cleanDicts(newL))
            else:
                pass  # lists have been compared
        else:
            print '%s%s -> %s' % (preamble, oldTerm, newTerm)
    else:
        if oldTerm != newTerm:
            try:
                print '%s%s -> %s' % (preamble, oldTerm, newTerm)
            except UnicodeEncodeError:
                print '%s%s -> %s' % (preamble, oldTerm.encode('utf-8'), newTerm.encode('utf-8'))


def termCompare(facDict, updatedFacDict, abbrev, key):
    """
    Compare two top-level dictionary entries.
    """
    if _VERBOSE:
        print 'comparing %s for %s' % (key, abbrev)
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

def loadByFName(dirPath):
    resultD = {}
    keyS = set()
    for nm in os.listdir(dirPath):
        if nm.endswith('.yaml'):
            fullPath = os.path.join(dirPath, nm)
            try:
                with open(fullPath, 'r') as f:
                    newD = yaml.safe_load(f)
                    if not isinstance(newD, dict):
                        print "%s is not a dict" % fullPath
                        continue
                    keyS.update(newD.keys())
                    resultD[nm] = newD
            except IOError, e:
                print 'Could not open %s: %s' % (fullPath, e)
            except Exception, e:
                print 'Failed to parse %s: %s' % (fullPath, e)
    return list(keyS), resultD

def main():
    global _VERBOSE
    
    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)

    parser = optparse.OptionParser(usage="""
    %prog [-v] [--byfilename] yamlDir1 yamlDir2
    """)
    parser.add_option("-v", "--verbose", action="store_true",
                      help="verbose output")
    parser.add_option("--byfilename", action="store_true",
                      help="match on filename rather than abbrev")

    opts, args = parser.parse_args()  # @UnusedVariable

    if opts.verbose:
        _VERBOSE = True

    if len(args) == 2:
        dir1Path = args[0]
        dir2Path = args[1]
    else:
        parser.error("Two directories containing YAML collections must be specified.")

    if opts.byfilename:
        dummy, facDict = loadByFName(dir1Path)
        dummy, updatedFacDict = loadByFName(dir2Path)
    else:
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

