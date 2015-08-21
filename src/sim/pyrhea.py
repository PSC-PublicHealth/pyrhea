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

_rhea_svn_id_ = "$Id$"

import sys
import os
import yaml
from imp import load_source
import jsonschema
from random import seed
from collections import defaultdict
import optparse

import patches
import yaml_tools

inputSchema = 'rhea_input_schema.yaml'


def loadFacilityImplementations(implementationDir):
    print 'Loading facility implementations'
    implDict = {}
    sys.path.append(implementationDir)
    for fname in os.listdir(implementationDir):
        name, ext = os.path.splitext(fname)
        if ext == '.py':
            newMod = load_source(name, os.path.join(implementationDir, fname))
            for requiredAttr in ['category', 'generateFull']:
                if not hasattr(newMod, requiredAttr):
                    raise RuntimeError('The facility implemenation in %s has no %s' %
                                       (os.path.join(implementationDir, fname),
                                        requiredAttr))
            assert newMod.category not in implDict, \
                "Redundant definitions for category %s" % newMod.category
            implDict[newMod.category] = newMod
            print 'Loaded %s implementing %s' % (fname, newMod.category)
        elif ext.startswith('.py'):
            pass  # .pyc files for example
        else:
            print 'Skipping non-python file %s' % fname
    sys.path.pop()  # drop implementaionDir
    return implDict


def loadFacilityDescriptions(dirList, facImplDict):
    """
    This loads facilities, including communities, checking them against a schema.
    """
    print 'Parsing facilities descriptions'
    allKeySet = set()
    rawRecs = []
    for d in dirList:
        kS, recs = yaml_tools.parse_all(d)
        rawRecs.extend(recs)
        allKeySet.update(kS)
    facRecs = []
    for rec in rawRecs:
        assert 'abbrev' in rec, "Facility description has no 'abbrev' field: %s" % rec
        assert 'category' in rec, \
            "Facility description for %(abbrev) has no 'category' field" % rec
        nErrors = facImplDict[rec['category']].checkSchema(rec)
        if nErrors:
            print 'dropping %s; %d schema violations' % (rec['abbrev'], nErrors)
        else:
            facRecs.append(rec)
    return facRecs


def distributeFacilities(comm, facDirs, facImplDict):
    """
    Estimate work associated with each facility and attempt to share it out fairly
    """
    if comm.rank == 0:
        facRecs = loadFacilityDescriptions(facDirs, facImplDict)
        assignments = defaultdict(list)
        tots = {k: 0 for k in xrange(comm.size)}
        pairL = [(facImplDict[rec['category']].estimateWork(rec), rec) for rec in facRecs]
        pairL.sort(reverse=True)
        for newWork, rec in pairL:
            leastWork, leastBusy = min([(v, k) for k, v in tots.items()])  # @UnusedVariable
            assignments[leastBusy].append(rec)
            tots[leastBusy] += newWork
        print 'Estimated work by rank: %s' % tots
        for targetRank in xrange(comm.size):
            if targetRank == comm.rank:
                myFacList = assignments[targetRank]
            else:
                comm.send(assignments[targetRank], dest=targetRank)
    else:
        myFacList = comm.recv(source=0)
    if comm.rank == 0:
        print 'Finished distributing facilities'
    return myFacList


class TweakedOptParser(optparse.OptionParser):
    def setComm(self, comm):
        self.comm = comm

    def exit(self, code, msg):
        print msg
        if hasattr(self, 'comm'):
            self.comm.Abort(code)
        else:
            sys.exit(code)


def checkInputFileSchema(fname, schemaFname, comm):
    try:
        with open(fname, 'rU') as f:
            inputJSON = yaml.safe_load(f)
        with open(os.path.join(os.path.dirname(__file__), schemaFname), 'rU') as f:
            schemaJSON = yaml.safe_load(f)
        validator = jsonschema.validators.validator_for(schemaJSON)(schema=schemaJSON)
        nErrors = sum([1 for e in validator.iter_errors(inputJSON)])  # @UnusedVariable
        if nErrors:
            print 'Input file violates schema:'
            for e in validator.iter_errors(inputJSON):
                print e
            comm.Abort(2)
        else:
            return inputJSON
    except Exception, e:
        print 'Error checking input against its schema: %s' % e
        comm.Abort(2)


def main():

    comm = patches.getCommWorld()

    if comm.rank == 0:
        parser = TweakedOptParser(usage="""
        %prog [-v][-d][-t][-D] input.yaml
        """)
        parser.setComm(comm)
        parser.add_option("-v", "--verbose", action="store_true",
                          help="verbose output")
        parser.add_option("-d", "--debug", action="store_true",
                          help="debugging output")
        parser.add_option("-t", "--trace", action="store_true",
                          help="enable greenlet thread tracing")
        parser.add_option("-D", "--deterministic", action="store_true",
                          help="deterministic mode")

        opts, args = parser.parse_args()
        clData = {'verbose': opts.verbose,
                  'debug': opts.debug,
                  'trace': opts.trace,
                  'deterministic': opts.deterministic
                  }
        if len(args) == 1:
            clData['input'] = checkInputFileSchema(args[0], inputSchema, comm)
        else:
            parser.error("A YAML-format file specifying run parameters must be specified.")
        parser.destroy()
    else:
        clData = None
    clData = comm.bcast(clData, root=0)

    verbose = clData['verbose']  # @UnusedVariable
    debug = clData['debug']  # @UnusedVariable
    trace = clData['trace']
    deterministic = clData['deterministic']
    inputDict = clData['input']

    if deterministic:
        seed(1234 + comm.rank)  # Set the random number generator seed

    facImplDict = loadFacilityImplementations(inputDict['facilityImplementationDir'])

    myFacList = distributeFacilities(comm, inputDict['facilityDirs'], facImplDict)
    print 'Rank %d has facilities %s' % (comm.rank, [rec['abbrev'] for rec in myFacList])

    patchGroup = patches.PatchGroup(comm, trace=trace, deterministic=deterministic)
    # Only one patch per rank for now
    patch = patchGroup.addPatch(patches.Patch(patchGroup))

    allIter = []
    allAgents = []
    for facDescription in myFacList:
        if facDescription['category'] in facImplDict:
            facilities, wards, patients = \
                facImplDict[facDescription['category']].generateFull(facDescription, patch)
            allIter.extend([fac.reqQueue for fac in facilities])
            allIter.extend([fac.holdQueue for fac in facilities])
            allIter.extend(wards)
            allAgents.extend([fac.manager for fac in facilities])
            allAgents.extend(patients)
        else:
            raise RuntimeError('Facility %(abbrev)s category %(category)s has no implementation' %
                               facDescription)

    patch.addInteractants(allIter)
    patch.addAgents(allAgents)
    print 'Rank %d: %d interactants, %d agents' % (comm.rank, len(allIter), len(allAgents))

    #patchGroup.start()
    print '%s all done (from main)' % patchGroup.name


############
# Main hook
############

if __name__ == "__main__":
    main()
