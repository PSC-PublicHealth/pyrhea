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

import sys
import os.path
import signal
import optparse
import yaml
import cPickle as pickle
import matplotlib.cm as pltcm

import schemautils
import pyrheautils
from map_transfer_matrix import parseFacilityData
import phacsl.utils.formats.csv_tools as csv_tools
import phacsl.utils.formats.yaml_tools as yaml_tools
import geomap
import geodatamanager

# http://stackoverflow.com/questions/18229628/python-profiling-using-line-profiler-clever-way-to-remove-profile-statements
try:
    profile = __builtins__.profile
except AttributeError:
    # No line profiler, provide a pass-through version
    def profile(func):
        """
        Dummy substitute for __builtin__.profile line profiler
        """
        return func


LTRED = '#cc6666'
RED = '#cc0000'
LTMAGENTA = '#cc66cc'
MAGENTA = '#cc00cc'
LTBLUE = '#6666cc'
BLUE = '#0000cc'
ALTBLUE = '#6699cc'
LTCYAN = '#66cccc'
CYAN = '#00cccc'
LTGREEN = '#66cc66'
GREEN = '#00cc00'
LTYELLOW = '#cccc66'
YELLOW = '#cccc00'

SCHEMA_DIR = '../schemata'
INPUT_SCHEMA = 'rhea_input_schema.yaml'

DEFAULT_NOTES_FNAME = 'notes.pkl'

def categoryToColor(cat):
    return {'HOSPITAL': 'red', 'LTAC': 'green', 'NURSINGHOME': 'blue', 'COMMUNITY': 'pink',
            'SNF': 'blue', 'VSNF': 'purple'}[cat]


def loadCountyCodes(fname):
    with open(fname, 'rU') as f:
        keys, recs = csv_tools.parseCSV(f)
    for k in ['STATE', 'STATEFP', 'COUNTYNAME', 'COUNTYFP']:
        assert k in keys, "%s is missing the required column %s" % k
    result = {}
    for r in recs:
        countyNm = r['COUNTYNAME'].strip()
        if countyNm.lower().endswith('county'):
            countyNm = countyNm[:-len('county')].strip().lower()
        stateNm = r['STATE'].strip().upper()
        key = (stateNm, countyNm)
        result[key] = (r['STATEFP'], r['COUNTYFP'])
    return result


def drawMap(geoDataPathList, locRecs, stateFIPSRE, countyFIPSRE, countySet):

    sumLon = 0.0
    sumLat = 0.0
    for r in locRecs:
        if len(r['name']) > 0:
            try:
                sumLon += float(r['longitude'])
                sumLat += float(r['latitude'])
            except Exception, e:
                print 'Got exception %s on %s' % (e, r)

    ctrLon = sumLon / len(locRecs)
    ctrLat = sumLat / len(locRecs)
    ctrLatRad = 0.5*(41.615453 + 41.419047)
    ctrLonRad = 0.5*(-87.575672 + -87.180045)

    myMap = geomap.Map(geoDataPathList,
                       stateFIPSRE,
                       countyFIPSRE,
                       ctrLon, ctrLat,
                       annotate=False,
                       annotateTracts=False,
                       nameMap={'STATE': 'STATE', 'COUNTY': 'COUNTY',
                                'TRACT': 'TRACT', 'NAME': 'NAME',
                                'GEOKEY': 'GEO_ID'},
                       regexCodes=True,
                       mrkSzDict={'*': 15, '+': 15, 'o': 10})
    clrTupleSeq = [(LTRED, RED), (LTMAGENTA, MAGENTA), (LTBLUE, BLUE),
                   (LTCYAN, CYAN), (LTGREEN, GREEN), (LTYELLOW, YELLOW)]
    censusTractFacD = {r['censusTract']:r for r in locRecs if 'censusTract' in r}
    minPopDensity = None
    maxPopDensity = None
    for geoID in myMap.tractIDList():
        propRec = myMap.getTractProperties(geoID)
        area = propRec['CENSUSAREA']
        if propRec['TRACT'] in censusTractFacD and area >= 0.1:  # exclude tiny tracts
            popDens = censusTractFacD[propRec['TRACT']]['meanPop']['value']/propRec['CENSUSAREA']
            if not minPopDensity or minPopDensity > popDens:
                minPopDensity = popDens
            if not maxPopDensity or maxPopDensity < popDens:
                maxPopDensity = popDens
    cmap = pltcm.get_cmap('plasma')
    for geoID in myMap.tractIDList():
        propRec = myMap.getTractProperties(geoID)
        countyKey = (int(propRec['STATE']), int(propRec['COUNTY']))
        if countyKey in countySet:
            # This block controls coloring of census tracts
            # To color according to the color tuple sequence:
            # clr, dummy = clrTupleSeq[int(propRec['COUNTY']) % len(clrTupleSeq)]
            #
            # To color according to population density:
            if propRec['TRACT'] in censusTractFacD:
                popDens = censusTractFacD[propRec['TRACT']]['meanPop']['value']/propRec['CENSUSAREA']
                clr = cmap(popDens/maxPopDensity)
            else:
                clr = cmap(0.0)
            # End of block controlling coloring of census tracts
                
            myMap.plotTract(geoID, clr)

    for r in locRecs:
        name = r['abbrev']
        if 'censusTract' in r:
            pass  # Already drew these while plotting tracts above
        else:
            try:
                mrk = {'SNF': '^', 'VSNF': '^', 'NURSINGHOME': '^', 'HOSPITAL': '*',
                       'LTAC': 'o', 'LTACH': 'o', 'CHILDREN': 'D'}[r['category']]
            except KeyError:
                print 'No marker type for %s' % r['category']
                mrk = 'x'
            if name:
                myMap.plotMarker(float(r['longitude']), float(r['latitude']), mrk, name, 'black',
                                 sz=0.1*r['meanPop']['value'])

    myMap.draw(mrkLegendDict={'*': 'Hospitals', '^': 'SNFs, VSNFs', 'o': 'LTACHs', 'D': "Children's"})


def checkInputFileSchema(fname, schemaFname):
    try:
        with open(fname, 'rU') as f:
            inputJSON = yaml.safe_load(f)
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
    except Exception, e:
        sys.exit('Error checking input against its schema: %s' % e)


def importNotes(fname):
    with open(fname, 'r') as f:
        stuff = pickle.load(f)
    return stuff


def main():
    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)
    

    parser = optparse.OptionParser(usage="""
    %prog [-n notes.pkl] run_descr.yaml
    """)

    parser.add_option('-n', '--notes', action='store', type='string',
                      help="Notes filename (overrides any name in the run description)")
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error('A YAML run description is required')

    parser.destroy()

    inputDict = checkInputFileSchema(args[0],
                                     os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    if 'modelDir' in inputDict:
        pyrheautils.PATH_STRING_MAP['MODELDIR'] = pyrheautils.pathTranslate(inputDict['modelDir'])
    implDir = inputDict['facilityImplementationDir']

    if opts.notes:
        notesFName = opts.notes
    elif 'notesFileName' in inputDict:
        notesFName = pyrheautils.pathTranslate(inputDict['notesFileName'])
    else:
        notesFName = DEFAULT_NOTES_FNAME

    notesDict = importNotes(notesFName)

    gDM = geodatamanager.GeoDataManager()

    facDirList = [pyrheautils.pathTranslate(fD) for fD in inputDict['facilityDirs']]
    facDict = parseFacilityData(facDirList)
    countySet = set()
    for rec in facDict.values():
        if 'FIPS' in rec:
            fipsCode = rec['FIPS']
            countySet.add((int(fipsCode[:2]), int(fipsCode[2:])))
            gDM.addFIPS(fipsCode)
        else:
            line = rec['address']
            strs = line.split(',')
            if 'USA' in strs[-1]:
                strs = strs[:-1]
            if 'county' in strs[-2].lower():
                stateStr = strs[-1].split()[0].upper()
                countyStr = strs[-2].split()[0].lower()
                try:
                    gDM.addCounty(stateStr, countyStr)
                except:
                    print 'Unknown county %s %s' % (countyStr, stateStr)
    countyFIPSRE = r'.*'
    stateFIPSRE = r'.*'


#     if len(args) == 1:
#         inputPath = args[0]
#     else:
#         parser.error("A single CSV file containing coordinate data must be specified.")



#     for cg, mrk, sz in {('HOSPITAL', '*', 15), ('LTAC', '+', 15), ('NURSINGHOME', 'o', 10)}:
#         coordList = [(mapProj.project(rec['longitude'], rec['latitude']),
#                      rec['abbrev'])
#                      for rec in facDict.values() if rec['category'] == cg]
#         coordList = [(x, y, abbrev) for (x, y), abbrev in coordList]
#         ax.scatter([x for x, y, abbrev in coordList], [y for x, y, abbrev in coordList],
#                    c=categoryToColor(cg), marker=mrk)
#         for x, y, abbrev in coordList:
#             ax.annotate(abbrev, (x, y), (x, y), ha='left', va='bottom', fontsize=sz)


#     with open(inputPath, 'rU') as f:
#         keys, recs = parseCSV(inputPath)
#     for k in ['location ID', 'latitude', 'longitude']:
#         assert k in keys, "Required CSV column %s is missing" % k
#         
#     geoDataPath = '/home/welling/geo/India/india.geojson'
#     stateCode = '06'
#     countyCode = '059'

    drawMap(gDM.getGeoDataFileList(), facDict.values(), stateFIPSRE, countyFIPSRE, countySet)

if __name__ == "__main__":
    main()
