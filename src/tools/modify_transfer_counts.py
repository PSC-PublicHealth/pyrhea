import os.path
import sys

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))


import re
import yaml
from optparse import OptionParser
import tools_util as tu
import copy


def getLocalOccupancyDataList(notesDataList):
    ret = []
    for n in notesDataList:
        for k in n.keys():
            if k.startswith('Patch_'):
                ret.append(n[k]['localoccupancy'])

    return ret


def getMaxDays(dataList):
    ret = -1
    for d in dataList:
        ret = max(len(d), ret)

    return ret

def getSampleDayList(maxDay):
    "get a list of 5 days spaced evenly over the last third of the run"
    dayList = []
    dayList.append(maxDay * 2 / 3 - 1)    # 8/12 of the way through
    dayList.append(maxDay * 3 / 4 - 1)    # 9/12
    dayList.append(maxDay * 5 / 6 - 1)    # 10/12
    dayList.append(maxDay * 11 / 12 - 1)  # 11/12
    dayList.append(maxDay - 1)            # 12/12
    return dayList

def prepOccupancyData(notesDataList):
    lodl = getLocalOccupancyDataList(notesDataList)
    maxDay = getMaxDays(lodl)
    dayList = getSampleDayList(maxDay)

    return lodl, dayList

def getDayStats(fac, day, dataList):
    ret = []
    for d in dataList:
        try:
            val = d[day][fac]
            ret.append(val)
        except:
            pass

    return ret

def getPopData(fac, facDict, lodl, dayList):
    popVals = []
    
    for day in dayList:
        popVals.extend(getDayStats(fac, day, lodl))
    try:
        aPop = float(sum(popVals)) / len(popVals) # mean actual population
    except:
        aPop = None

    try:
        ePop = facDict[fac]['meanPop']['value'] # expected population
    except:
        ePop = None

    print "%s: expected %s, actual %s"%(fac, ePop, aPop)
    return ePop, aPop


def loadTransferCounts():
    with open("transfer_counts.yaml") as f:
        d = yaml.load(f)
    return d

def adjustCounts(fac, ratio, xd):
    for loc,table in xd.items():
        for dest, count in table.items():
            if dest == fac:
                table[dest] = float(count) * ratio

def main():
    """
    main
    """
#    import sys
#    print sys.argv
    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl] [--glob] [--out outname.yaml] run_descr.yaml
    """)
    parser.add_option('-n', '--notes', action='append', type='string',
                      help="Notes filename - may be repeated")
    
    #parser.add_option('-o', '--out', action='store', type='string',
    #                  help="Output file (defaults to %s)" % DEFAULT_OUT_FILE)
    
    parser.add_option('--glob', action='store_true',
                      help=("Apply filename globbing for notes files."
                            "  (Remember to protect the filename string from the shell!)"))
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    if not opts.notes:
        parser.error('At least one --notes option is required')

    parser.destroy()
    runDesc = args[0]

    inputDict = tu.readModelInputs(runDesc)
    #print inputDict.keys()
    facDict = tu.getFacDict(inputDict)

    notesDataList = tu.getNotesData(opts.notes, opts.glob)
    xferCounts = loadTransferCounts()
    xferData = []
    for i in xrange(10):
        xferData.append(copy.deepcopy(xferCounts))

    lodl, dayList = prepOccupancyData(notesDataList)
    facList = inputDict['trackedFacilities']
    for fac in facList:
        ePop, aPop = getPopData(fac, facDict, lodl, dayList)
        try:
            popRatio = ePop / aPop
        except:
            continue

        print popRatio

        adjRatio = 0.1
        for xd in xferData:
            facAdj = (1.0 - popRatio) * adjRatio + 1.0
            adjustCounts(fac, facAdj, xd)
            adjRatio += 0.1

    for i,xd in enumerate(xferData):
        with open("transfer_counts_%s.yaml"%i, "w") as f:
            yaml.dump(xd, f)

if __name__ == "__main__":
    main()
