#! /usr/bin/env python

"""
This assumes that lat/lon is present the yaml files, and uses Google's distancematrix
API to grab a section of the travel time/ distance matrix.
"""

import sys
import os
import urllib
import json
import time
import phacsl.utils.formats.yaml_tools as yaml_tools
import phacsl.utils.formats.csv_tools as csv_tools
import datetime


class PST(datetime.tzinfo):
    def utcoffset(self, dt):
        return datetime.timedelta(hours=-8)

    def dst(self, dt):
        return datetime.timedelta(0)


class UTC(datetime.tzinfo):
    def utcoffset(self, dt):
        return datetime.timedelta(hours=0)

    def dst(self, dt):
        return datetime.timedelta(0)

targetDateTime = datetime.datetime(2015, 7, 11, 0, tzinfo=PST())  # 11AM on 7/11/2015
epochDateTime = datetime.datetime(1970, 1, 1, tzinfo=UTC())  # Unix epoch
targetTimestamp = (targetDateTime - epochDateTime).total_seconds()

nRequestsThisBatch = 52  # Google API permits 100, but this keeps url length down

facKeys, facRecs = yaml_tools.parse_all('/home/welling/workspace/pyRHEA/models/OrangeCounty2013/facilityfactsCurrent2013')
facDict = {r['abbrev']: r for r in facRecs}

try:
    with open('/home/welling/workspace/pyRHEA/models/OrangeCounty2013/transitmatrix.csv', 'rU') as f:
        transitKeys, transitRecs = csv_tools.parseCSV(f)  #
except Exception, e:
    print 'Could not parse input transitmatrix'
    transitKeys = ['From', 'To', 'Seconds', 'Meters']
    transitRecs = []
nInputRecs = len(transitRecs)

transitDict = {}
for r in transitRecs:
    if r['From'] not in transitDict:
        transitDict[r['From']] = {}
    transitDict[r['From']][r['To']] = {'Seconds': r['Seconds'], 'Meters': r['Meters']}

newRecs = []

fromLoc = None
for k in facDict.keys():
    if k not in transitDict:
        fromLoc = k
        transitDict[fromLoc] = {}
        break
    elif len(transitDict[k]) < len(facDict):
        fromLoc = k
        break
if fromLoc is None:
    print "All transit data is already present"
else:
    toLocList = []
    nFetch = 0
    for loc in facDict.keys():
        if loc not in transitDict[fromLoc]:
            toLocList.append(loc)
            nFetch += 1
            if nFetch >= nRequestsThisBatch:
                break
print 'Fetching data from %s' % fromLoc
print '...to %s' % toLocList

for k in ['latitude', 'longitude']:
    assert k in facDict[fromLoc], '%s is missing %s data' % (fromLoc, k)

try:
    query = {'origins': '%f,%f' % (facDict[fromLoc]['latitude'], facDict[fromLoc]['longitude'])}
    dstStr = ""
    for toLoc in toLocList:
        for k in ['latitude', 'longitude']:
            assert k in facDict[toLoc], '%s is missing %s data' % (toLoc, k)
        dstStr += '|%f,%f' % (facDict[toLoc]['latitude'], facDict[toLoc]['longitude'])
    dstStr = dstStr[1:]
    query = {'arrival_time': int(targetTimestamp),
             'origins': '%f,%f' % (facDict[fromLoc]['latitude'], facDict[fromLoc]['longitude']),
             'destinations': dstStr}
    url = ('http://maps.googleapis.com/maps/api/distancematrix/json?%s' %
           urllib.urlencode(query))
except Exception, e:
    print 'bad encode: %s' % e
    print 'query: %s' % query
    print 'url: %s' % url

try:
    data = None
    time.sleep(1.0)
    response = urllib.urlopen(url)
    data = json.loads(response.read())
except Exception, e:
    print 'bad transaction with googleapis: %s' % e
    print 'returned JSON: %s' % data

for d, toLoc in zip(data['rows'][0]['elements'], toLocList):
    if d['status'] == 'OK':
        transitDict[fromLoc][toLoc] = {'Seconds': d['duration']['value'],
                                       'Meters': d['distance']['value']}
    else:
        print 'Status = %s for %s -> %s' % (d['status'], fromLoc, toLoc)

outTransitRecs = []
for fromLoc, fromDict in transitDict.items():
    for toLoc, d in fromDict.items():
        outRec = {'From': fromLoc, 'To': toLoc}
        outRec.update(d)
        outTransitRecs.append(outRec)

print 'Added %d recs; %d total' % (len(outTransitRecs) - nInputRecs, len(outTransitRecs))
try:
    os.unlink('/home/welling/workspace/pyRHEA/models/OrangeCounty2013/transitmatrix.sav.csv')
except:
    pass
try:
    os.rename('/home/welling/workspace/pyRHEA/models/OrangeCounty2013/transitmatrix.csv',
              '/home/welling/workspace/pyRHEA/models/OrangeCounty2013/transitmatrix.sav.csv')
except:
    pass
with open('/home/welling/workspace/pyRHEA/models/OrangeCounty2013/transitmatrix.csv', 'w') as f:
    csv_tools.writeCSV(f, transitKeys, outTransitRecs)
print 'Wrote %s' % '/home/welling/workspace/pyRHEA/models/OrangeCounty2013/transitmatrix.csv'
