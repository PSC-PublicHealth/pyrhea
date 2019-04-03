import os
import sys
import yaml
from shutil import copytree

import pyrheautils
from lockserver import Lock

CACHE_TO_LOCAL_DISK = False

def mvCache():
    cache = pyrheautils.pathTranslate('$(COMMUNITYCACHEDIR)')
    if CACHE_TO_LOCAL_DISK:
        if 'SLURM_ARRAY_JOB_ID' in os.environ:
            lockName = '%s:%s' % (os.environ['SLURM_ARRAY_JOB_ID'], os.environ['HOSTNAME'])
            infoFName = '/tmp/info_%s' % os.environ['SLURM_ARRAY_JOB_ID']
        else:
            lockName = '%s:%s' % (os.environ['SLURM_JOB_ID'], newCache)
            infoFName = '/tmp/info_%s' % os.environ['SLURM_JOB_ID']
        newCache = None
        with Lock(lockName):
            doTheCopy = False
            if os.path.exists(infoFName):
                print 'reading ', infoFName
                with open(infoFName, 'rU') as f:
                    newCache = f.readline().strip()
                    print 'newCache is ', newCache
                if not os.path.exists(newCache):
                    # old copy is bad
                    print 'this newCache is empty; trying again'
                    doTheCopy = True
            else:
                doTheCopy = True
                print 'no info file found'

            if doTheCopy:
                newCache = os.path.join(os.environ['LOCAL'], os.path.split(cache)[1])
                with open(infoFName, 'w') as f:
                    f.write(newCache)
                print 'wrote new cache location %s to %s' % (newCache, infoFName)
                copytree(cache, newCache)
                print 'finished copying'
            else:
                print 'Someone else copied the cache to %s' % newCache
    else:
        newCache = os.path.join('/pylon5/pscstaff/welling/pyrhea/caches',
                                str(os.environ['SLURM_JOBID']), os.path.split(cache)[1])
        print "copying cache from %s to %s"%(cache, newCache)
        copytree(cache, newCache)
        print "finished making copy"

    pyrheautils.PATH_STRING_MAP['COMMUNITYCACHEDIR'] = newCache
                
    # And tell pyrhea that we moved this:
    assert newCache is not None, 'Somehow we avoided setting newCache'
    pyrheautils.PATH_STRING_MAP['COMMUNITYCACHEDIR'] = newCache

mvCache()

with open('scenario_custom.yaml', 'rU') as f:
    allData = yaml.load(f)
assert 'constantsReplacementData' in allData, 'scenario custom yaml has no constantsReplacementData'
assert 'facilitiesReplacementData' in allData, 'scenario custom yaml has no facilitiesReplacementData'

constantsReplacementData = allData['constantsReplacementData']
facilitiesReplacementData = allData['facilitiesReplacementData']

print 'constantsReplacementData: ', constantsReplacementData
print 'facilitiesReplacementData: ', facilitiesReplacementData

