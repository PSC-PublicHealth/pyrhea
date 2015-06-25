#! /usr/bin/env python

"""
This tool generates ready-to-use facility information in yaml format from facilityfacts data,
also in yaml format.
"""

import types
import math
import yaml_tools

"""An empirical estimate"""
bedsPerWard = 20
"""Another empirical estimate"""
bedsPerICUWard = 12


def estimateHospitalWards(meanPop, fracAdultPatientDaysICU):
    meanICUPop = meanPop * fracAdultPatientDaysICU
    meanHOSPPop = meanPop - meanICUPop
    nICUWards = int(math.ceil(meanICUPop / bedsPerICUWard))
    icuHeadroom = (nICUWards * bedsPerICUWard) - meanICUPop
    nHOSPWards = int(math.ceil(meanHOSPPop / bedsPerWard))
    hospHeadroom = (nHOSPWards * bedsPerWard) - meanHOSPPop
    print 'headrooms: %s %s' % (hospHeadroom, icuHeadroom)
    return nHOSPWards, nICUWards


def createWard(abbrev, tier, nBeds, counter):
    return {'name': 'Ward_%s_%s_%d' % (abbrev, tier, counter),
            'nBeds': nBeds,
            'tier': tier
            }

allKeySet, factRecs = yaml_tools.parse_all('facilityfacts7')

facilityRecs = []
for factR in factRecs:
    try:
        facility = {'name': factR['name'],
                    'abbrev': factR['abbrev'],
                    'category': factR['category'],
                    'initialPop': factR['meanPop']}

        if factR['category'] == 'HOSPITAL':
            assert 'meanPop' in factR and isinstance(factR['meanPop'], types.FloatType), \
                'meanPop is missing or invalid'
            nWards, nICUWards = estimateHospitalWards(factR['meanPop'],
                                                      factR['fracAdultPatientDaysICU'])
            facility['wards'] = []
            for i in xrange(nWards):
                facility['wards'].append(createWard(factR['abbrev'], 'HOSP', bedsPerWard, i))
            for i in xrange(nICUWards):
                facility['wards'].append(createWard(factR['abbrev'], 'ICU', bedsPerICUWard, i))

        elif factR['category'] == 'NURSINGHOME':
            assert 'nBeds' in factR and isinstance(factR['nBeds'], types.IntType), \
                'nBeds is missing or invalid'
            facility['wards'] = [createWard(factR['abbrev'], 'NURSING', factR['nBeds'], 0)]

        elif factR['category'] == 'LTAC':
            # LTACs are like HOSPITALS but are treated as having no ICUs
            assert 'meanPop' in factR and isinstance(factR['meanPop'], types.FloatType), \
                'meanPop is missing or invalid'
            nWards = int(math.ceil(factR['meanPop']/float(bedsPerWard)))
            facility['wards'] = []
            for i in xrange(nWards):
                facility['wards'].append(createWard(factR['abbrev'], 'HOSP', bedsPerWard, i))
            headRoom = bedsPerWard*nWards - factR['meanPop']
            print '%s headroom: %s' % (factR['abbrev'], headRoom)

        else:
            raise RuntimeError('Unknown category %s' % factR['category'])

        facilityRecs.append(facility)
    except Exception, e:
        print 'Excluding %s: <%s>' % (factR['abbrev'], e)

yaml_tools.save_all('facilities', facilityRecs)
