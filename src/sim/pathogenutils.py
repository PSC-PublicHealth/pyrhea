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

import logging
from facilitybase import CareTier

logger = logging.getLogger(__name__)

def valFromCategoryEntry(key, ctg, constantsJSON):
    if key not in constantsJSON:
        raise RuntimeError('Constant list for %s was not found' % key)
    for val in constantsJSON[key]:
        if val['category'] == ctg:
            return val['frac']['value']
    raise RuntimeError('Constant entry for category %s was not found' % ctg)

def parseValByTier(fieldStr, constants):
    topD = {}
    for elt in constants[fieldStr]:
        topD[elt['tier']] = elt['value']
    return topD

def parseFracByTierByFacilityByCategory(fieldStr, constants):
    topD = {}
    for elt in constants[fieldStr]:
        cat = elt['category']
        facilities = elt['facilities']
        if cat not in topD:
            topD[cat] = {}
        for facElt in facilities:
            abbrev = facElt['abbrev']
            tiers = facElt['tiers']
            if abbrev not in topD[cat]:
                topD[cat][abbrev] = {}
            for tierElt in tiers:
                tier = tierElt['tier']
                val = tierElt['frac']['value']
                assert tier not in topD[cat][abbrev], ('Redundant %s for %s %s' %
                                                      (fieldStr, abbrev, CareTier.names[tier]))
                topD[cat][abbrev][tier] = val
    return topD


def parseFracByTierByCategory(fieldStr, constants, eltKey='frac'):
    topD = {}
    for elt in constants[fieldStr]:
        cat = elt['category']
        tiers = elt['tiers']
        if cat not in topD:
            topD[cat] = {}
        for tierElt in tiers:
            tier = tierElt['tier']
            val = tierElt[eltKey]['value']
            assert tier not in topD[cat], ('Redundant %s for %s %s' %
                                           (fieldStr, cat, CareTier.names[tier]))
            topD[cat][tier] = val
    return topD

def parseScaleByTierByCategory(fieldStr, constants):
    return parseFracByTierByCategory(fieldStr, constants, eltKey='scale')

def parseValByTierByCategory(fieldStr, constants):
    return parseFracByTierByCategory(fieldStr, constants, eltKey='value')

def parseTierTierScaleList(fieldStr, constants):
    result = {}
    for elt in constants[fieldStr]:
        fmTier = CareTier.__dict__[elt['tierFrom']]
        toTier = CareTier.__dict__[elt['tierTo']]
        value = elt['scale']['value']
        result[(fmTier, toTier)] = value
    return result

def getValByTierByCategory(tbl, tblNameForErr, ward, wardCategory, overrideTbl=None,
                            default=None):
    """
    tbl is expected to have the structure tbl[category][tier] -> value
    overrideTbl is expected to have the structure overrideTbl[category][abbrev][tier] -> value
    """
    tierStr = CareTier.names[ward.tier]
    if overrideTbl:
        # Potential override values are stored by [category][abbrev][tierName]
        abbrev = ward.fac.abbrev
        if (wardCategory in overrideTbl
                and abbrev in overrideTbl[wardCategory]
                and tierStr in overrideTbl[wardCategory][abbrev]):
            return overrideTbl[wardCategory][abbrev][tierStr]
    if wardCategory in tbl and tierStr in tbl[wardCategory]:
        return tbl[wardCategory][tierStr]
    else:
        if default is not None:
            return default
        else:
            raise RuntimeError('No way to set %s for %s tier %s' %
                               (tblNameForErr, ward.fac.abbrev, CareTier.names[ward.tier]))

def getValByTier(tbl, tblNameForErr, ward, overrideTbl=None, default=None):
    """
    tbl is expected to have the structure tbl[tier] -> value
    overrideTbl is expected to have the structure overrideTbl[abbrev][tier] -> value
    """
    tierStr = CareTier.names[ward.tier]
    if overrideTbl:
        # Potential override values are stored by [category][abbrev][tierName]
        abbrev = ward.fac.abbrev
        if abbrev in overrideTbl and tierStr in overrideTbl[abbrev]:
            return overrideTbl[abbrev][tierStr]
    if tierStr in tbl:
        return tbl[tierStr]
    else:
        if default is not None:
            return default
        else:
            raise RuntimeError('No way to set %s for %s tier %s' %
                               (tblNameForErr, ward.fac.abbrev, CareTier.names[ward.tier]))

