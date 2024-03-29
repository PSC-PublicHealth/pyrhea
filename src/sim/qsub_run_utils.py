import sys
import os
import os.path
import yaml
import pickle
import time

from socket import gethostname
from getpass import getuser



def getDefaultPyrheaDir():
    host = gethostname()
    user = getuser()

    if user == 'jim' and host == 'nandi':
        return "/home/jim/projects/pyrhea/repo/pyrhea/"
    if user == 'jleonard':
        return "/home/jleonard/pyrhea/pyrhea/"

    raise RuntimeError("Can't determine pyrhea directory")


def getDefaultRunsDir():
    host = gethostname()
    user = getuser()

    if user == 'jim' and host == 'nandi':
        return "/home/jim/projects/pyrhea/repo/pyrhea/src/sim/runs"
    if user == 'jleonard':
        return "/home/jleonard/pyrhea/pyrhea/src/sim/runs"

    return "runs"



def replaceFile(base, out, repl):
    """
    creates a new file from a base file based on keyword substitution.
    base: input file name
    out: output file name
    repl: list of tuples of the form ("keyword", "substitution")
    """
    with open(base) as inf:
        with open(out, "w") as outf:
            for line in inf:
                for s,r in repl:
                    line = line.replace(s,r)
                outf.write(line)


def yamlReplaceFile(base, out, repl):
    """
    creates a new file based on yaml substitution.
    base: input file name
    out: output file name
    repl: list of tuples of the form (("path", 2, "yaml", "key"), "substitution")

    for lists of dictionaries where there is a key and a value in the dictionary, instead of needing
    to follow the index into the array use the following:
        '#KeyVal', <keyKey>, <keyValue>
    probably followed by 
        , <valueKey>, <valueValue>
    The first three will find the appropriate array index and traverse to that
    This will then traverse to <valueKey> (which is allowed to be the same as keyKey if you wish to change the key)
    and normally this would be followed by the value you wish to place in valueValue

    for instance in the main pyrhea yaml file 'pathTranslations' is an array of dicts and I want to change
    'CONSTANTS' value to something else (lets say 'new/constantsdir'):
    <keyKey> is 'key' <keyValue> is 'CONSTANTS', <valueKey> is 'value' and finally <valueValue> is 'new/constantsdir' so 
    I would use:
        '#KeyVal', 'key', 'CONSTANTS', 'value', 'new/constantsdir'
    """
    with open(base) as f:
        yData = yaml.load(f)

    for path, val in repl:
        pointer = yData
        idx = 0
        while idx < len(path) - 2:
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
                pointer = pointer[step]
                idx += 1
        pointer[path[-1]] = val
    with open(out, 'w') as f:
        yaml.dump(yData, f)


def writeReplPairs(runDir, rp):
    """
    takes a structure listing the replacement pairs and writes them to the run directory.
    This isn't used for anything except documenting after the fact what was done for this scenario.
    """

    rpYamlFile = os.path.join(runDir, 'cre_repl.yaml')
    with open(rpYamlFile, 'w') as f:
        f.write(yaml.dump(rp))
    


class RunEnvironment:
    def __init__(self, model, expName, expNum, instNum, baseConfFile,
                 replPairs=None, finalEditsFn=None, pyrheaDir=None, runsDir=None,
                 pyrheaOpts=None, pyrheaPrefix=None):
        if pyrheaDir is None:
            self.pyrheaDir = getDefaultPyrheaDir()
        else:
            self.pyrheaDir = pyrheaDir

        if runsDir is None:
            self.runsDir = getDefaultRunsDir()
        else:
            self.runsDir = runsDir

        self.model = model
        self.expName = expName
        self.expNum = expNum
        self.instNum = instNum
        self.baseConfFile = baseConfFile
        self.replPairs = replPairs
        self.finalEditsFn = finalEditsFn
        if pyrheaOpts is None:
            pyrheaOpts = ""
        self.pyrheaOpts = pyrheaOpts
        if pyrheaPrefix is None:
            pyrheaPrefix = ""
        self.pyrheaPrefix = pyrheaPrefix

    def buildEverything(self):
        self.setDirectories()
        self.buildRunEnvironment()

    def runSim(self, reallyRun=True):
        print 'Starting runSim'
        notesFile = os.path.join(self.runDir, "notes_%03d.pkl"%self.instNum)
        if os.path.isfile(notesFile):
            try:
                with open(notesFile) as f:
                    temp = pickle.load(f)
                return
            except:
                pass
        seed = self.expNum * 1000 + self.instNum

        runStr = "%s python pyrhea.py -o %s %s --seed %s %s"%(self.pyrheaPrefix, notesFile, self.pyrheaOpts, seed, self.confFile)
        print runStr
        if reallyRun:
            os.system(runStr)


    def setDirectories(self):
        self.modelDir = os.path.join(self.pyrheaDir, "models", self.model)
        self.baseConstDir = os.path.join(self.modelDir, "constants")

    def buildRunEnvironment(self): #, expNum, instNum, baseConfFile):
        runDir = os.path.join(self.runsDir, "%s_%06d"%(self.expName, self.expNum))
        self.runDir = runDir
        constantsDir = os.path.join(runDir, "constants")
        builtFile = os.path.join(runDir, "builtOk.txt")

        try:
            os.makedirs(constantsDir)
        except:
            for i in xrange(300):  # only wait up to 300 seconds for this to be built.
                if os.path.isfile(builtFile):
                    break
                time.sleep(1)
            # if we fall through this, some of the other stuff will probably crash.
            # I'd rather see a crash than an infinite loop, and there's probably not a great way of handling
            # some of this without lots of work.
                              

        self.confFile = os.path.join(runDir, os.path.basename(self.baseConfFile))

        if not os.path.isfile(builtFile):
            rp = self.replPairs
            if rp is None:
                rp = {}

            for root, dirs, files in os.walk(self.baseConstDir, followlinks=True):
                relPath = os.path.relpath(root, self.baseConstDir)
                if relPath == '.':
                    relPath = ''
                curConstDir = os.path.join(constantsDir, relPath)
                for d in dirs:
                    os.makedirs(os.path.join(curConstDir, d))

                for file in files:
                    f = os.path.join(relPath, file)
                    fullF = os.path.join(self.baseConstDir, f)
                    symF = os.path.join(constantsDir, f)
                    if f in rp:
                        yamlReplaceFile(fullF, symF, rp[f])
                    else:
                        os.system("cp %s %s"%(fullF, symF))
                        #os.symlink(fullF, symF)
                    
#            baseConfReplace = [[['pathTranslations', 0, 'value'], constantsDir]]
            baseConfReplace = [[['pathTranslations', '#KeyVal', 'key', 'CONSTANTS', 'value'], constantsDir]]
            yamlReplaceFile(self.baseConfFile, self.confFile, baseConfReplace)

            if self.finalEditsFn is not None:
                self.finalEditsFn(expNum, runDir, constantsDir)
            
            # tell any others waiting that we've built the directory
            with open(builtFile, "w") as f:
                f.write("done\n")


        

