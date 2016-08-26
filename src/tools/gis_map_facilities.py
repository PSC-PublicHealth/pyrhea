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

import os.path
import signal
import optparse

from map_transfer_matrix import parseFacilityData
import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.formats.csv_tools as csv_tools


import geomap

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
    print result
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

    myMap = geomap.Map(geoDataPathList,
                       stateFIPSRE,
                       countyFIPSRE,
                       ctrLon, ctrLat,
                       annotate=True,
                       nameMap={'STATE': 'STATE', 'COUNTY': 'COUNTY',
                                'TRACT': 'TRACT', 'NAME': 'NAME',
                                'GEOKEY': 'GEO_ID'},
                       regexCodes=True,
                       mrkSzDict={'*': 15, '+': 15, 'o': 10})
    clrTupleSeq = [(LTRED, RED), (LTMAGENTA, MAGENTA), (LTBLUE, BLUE),
                   (LTCYAN, CYAN), (LTGREEN, GREEN), (LTYELLOW, YELLOW)]
    for geoID in myMap.tractIDList():
        propRec = myMap.getTractProperties(geoID)
        countyKey = (int(propRec['STATE']), int(propRec['COUNTY']))
        if countyKey in countySet:
            clr, dummy = clrTupleSeq[int(propRec['COUNTY']) % len(clrTupleSeq)]
            myMap.plotTract(geoID, clr)

    for r in locRecs:
        name = r['abbrev']
        try:
            mrk = {'SNF': '+', 'VSNF': '+', 'NURSINGHOME': '+', 'HOSPITAL': '*',
                   'LTAC': 'o', 'LTACH': 'o', 'CHILDREN': 'x'}[r['category']]
        except KeyError:
            print 'No marker type for %s' % r['category']
            mrk = 'x'
        if name:
            myMap.plotMarker(float(r['longitude']), float(r['latitude']), mrk, name, 'black')

    myMap.draw()


def main():
    verbose = False
    debug = False
    # Thanks to http://stackoverflow.com/questions/25308847/attaching-a-process-with-pdb for this
    # handy trick to enable attachment of pdb to a running program
    def handle_pdb(sig, frame):
        import pdb
        pdb.Pdb().set_trace(frame)
    signal.signal(signal.SIGUSR1, handle_pdb)

    geoDataDir = '/home/welling/geo/USA'
    geoDataPathList = [os.path.join(geoDataDir, 'IL/tigr_2010_17.geojson'),
                       os.path.join(geoDataDir, 'IN/tigr_2010_18.geojson'),
                       os.path.join(geoDataDir, 'IA/tigr_2010_19.geojson'),
                       os.path.join(geoDataDir, 'WI/tigr_2010_55.geojson'),
                       ]
    countyCodePath = os.path.join(geoDataDir, 'fips_county_codes.csv')

    modelDir = os.path.join(os.path.dirname(__file__), '../../models/ChicagoLand/')
    facDict = parseFacilityData([
                                 os.path.join(modelDir, 'facilityfacts'),
#                                  os.path.join(modelDir, 'synthCommunities')
                                 ])
    countyCodeDict = loadCountyCodes(countyCodePath)
    cL = []
    for rec in facDict.values():
        line = rec['address']
        strs = line.split(',')
        stateStr = strs[3].split()[0].upper()
        countyStr = strs[2].split()[0].lower()
        cL.append((stateStr, countyStr))
    countiesInPlay = list(set(cL))
    print countiesInPlay
    countySet = set()
    for sname, cname in countiesInPlay:
        if (sname, cname) in countyCodeDict:
            countySet.add(countyCodeDict[(sname, cname)])
        else:
            print 'No code for <%s><%s>' % (sname, cname)
#     countyFIPSList = list(countySet)
#     countyFIPSRE = r'((' + r')|('.join(["%03d" % cF for cF in countyFIPSList]) + r'))'
#     stateFIPSRE = r'17'
    countyFIPSRE = r'.*'
    stateFIPSRE = r'.*'

    parser = optparse.OptionParser(usage="""
    %prog [-v][-d] input.csv
    """)
    parser.add_option("-v", "--verbose", action="store_true",
                      help="verbose output")
    parser.add_option("-d", "--debug", action="store_true",
                      help="debugging output")

    opts, args = parser.parse_args()

    if opts.debug:
        debug = True
    elif opts.verbose:
        verbose = True

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

    drawMap(geoDataPathList, facDict.values(), stateFIPSRE, countyFIPSRE, countySet)

if __name__ == "__main__":
    main()
