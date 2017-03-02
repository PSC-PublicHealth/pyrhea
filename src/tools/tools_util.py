import os.path
import sys
cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))

#import map_transfer_matrix as mtm
import types
import schemautils
import pyrheautils
import phacsl.utils.formats.yaml_tools as yaml_tools
import yaml
import cPickle as pickle

SCHEMA_DIR = '../schemata'
INPUT_SCHEMA = 'rhea_input_schema.yaml'


def checkInputFileSchema(fname, schemaFname):
    try:
        with open(fname, 'rU') as f:
            inputJSON = yaml.safe_load(f)
        if os.name != "nt":
            validator = schemautils.getValidator(os.path.join(os.path.dirname(__file__), schemaFname))
            nErrors = sum([1 for e in validator.iter_errors(inputJSON)])  # @UnusedVariable
            if nErrors:
                print 'Input file violates schema:'
                for e in validator.iter_errors(inputJSON):
                    print ('Schema violation: %s: %s' %
                           (' '.join([str(word) for word in e.path]), e.message))
                sys.exit('Input file violates schema')
            else:
                return inputJSON
        else:
            return inputJSON
    except Exception, e:
        sys.exit('Error checking input against its schema: %s' % e)


def parseFacilityData(fnameOrNameList):
    if not isinstance(fnameOrNameList, types.ListType):
        fnameOrNameList = [fnameOrNameList]
    facDict = {}
    for fn in fnameOrNameList:
        junkKeys, facRecs = yaml_tools.parse_all(fn)  # @UnusedVariable
        for r in facRecs:
            assert r['abbrev'] not in facDict, 'Redundant definitions for %s' % r['abbrev']
            facDict[r['abbrev']] = r
    return facDict


def getFacDict(inputDict):
    facDirList = [pyrheautils.pathTranslate(pth) for pth in inputDict['facilityDirs']]
    facDict = parseFacilityData(facDirList)
    return facDict

def readModelInputs(runDesc):
    schemautils.setSchemaBasePath(SCHEMA_DIR)
    inputDict = checkInputFileSchema(runDesc,
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    modelDir = inputDict['modelDir']
    pyrheautils.PATH_STRING_MAP['MODELDIR'] = modelDir
    implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])
    pyrheautils.PATH_STRING_MAP['IMPLDIR'] = implDir

    return inputDict



def getNotesFileList(notesPathList, globFlag=False):
    if globFlag:
        newNotesPathList = []
        for notesFName in notesPathList:
            print '%s yields %s' % (notesFName, glob.glob(notesFName))
            newNotesPathList += glob.glob(notesFName)
        newNotesPathList.sort()
        notesPathList = newNotesPathList

    return notesPathList

def readNotesFiles(notesFileList):
    notesDataList = []
    for n in notesFileList:
        with open(n) as f:
            notesDataList.append(pickle.load(f))

    return notesDataList

def getNotesData(notesPathList, globFlag=False):
    notesFileList = getNotesFileList(notesPathList, globFlag)
    notesData = readNotesFiles(notesFileList)

    return notesData
