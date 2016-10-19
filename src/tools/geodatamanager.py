#! /usr/bin/env python

import os.path

import phacsl.utils.formats.csv_tools as csv_tools

GEO_DATA_DIR = '/home/welling/geo/USA'
COUNTY_CODE_PATH = os.path.join(GEO_DATA_DIR, 'fips_county_codes.csv')

class GeoDataManager(object):
    def __init__(self):
        self.stateSet = set()
        self.namesToCodes, self.fipsToNames = self._loadCountyCodes(COUNTY_CODE_PATH)

    def _loadCountyCodes(self, fname):
        with open(fname, 'rU') as f:
            keys, recs = csv_tools.parseCSV(f)
        for k in ['STATE', 'STATEFP', 'COUNTYNAME', 'COUNTYFP']:
            assert k in keys, "%s is missing the required column %s" % k
        index = {}
        reverseIndex = {}
        for r in recs:
            countyNm = r['COUNTYNAME'].strip()
            if countyNm.lower().endswith('county'):
                countyNm = countyNm[:-len('county')].strip().lower()
            stateNm = r['STATE'].strip().upper()
            key = (stateNm, countyNm)
            fipsTuple = (r['STATEFP'], r['COUNTYFP'])
            index[key] = fipsTuple
            fipsStr = '%02d%03d' % fipsTuple
            reverseIndex[fipsStr] = key
        return index, reverseIndex

    def addFIPS(self, fipsCodeStr):
        try:
            stateNm, countyNm = self.fipsToNames[fipsCodeStr]
            stateFP, countyFP = self.namesToCodes[(stateNm, countyNm)]  # @UnusedVariable
            self.stateSet.add((stateNm, stateFP))
        except KeyError:
            raise RuntimeError('Invalid FIPS code %s' % fipsCodeStr)
        
    def addCounty(self, stateStr, countyName):
        try:
            stateFP, countyFP = self.namesToCodes[(stateStr, countyName)]
            self.stateSet.add((stateStr, stateFP))
        except KeyError:
            raise RuntimeError('Invalid state abbrev %s or county name %s' %
                               (stateStr, countyName))
    
    def getGeoDataFileList(self):
        result = []
        for abbrev, code in self.stateSet:
            result.append(os.path.join(GEO_DATA_DIR, '%s' % abbrev,
                                       'tigr_2010_%02d.geojson' % int(code)))
        return result