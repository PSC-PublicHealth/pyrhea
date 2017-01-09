#! /usr/bin/env python

"""
This tool reads two directories full of YAML facility descriptions, finds
the straight-line distances between each location in the first and all
locations in the second, and generates a table satisfying transfermatrix_schema
giving the distances to the nearest N locations.
"""

import math
import optparse
import yaml

import phacsl.utils.formats.yaml_tools as yaml_tools

DEFAULT_N_NEAREST = 10
DEFAULT_OUT_FNAME = "separation_table.yaml"


def longitudeLatitudeSep(lon1, lat1, lon2, lat2):
    "Inputs are in floating point degrees.  Returns separation in kilometers"
    scale = lat1r = lat2r = lon1r = lon2r = a = b = 0.0
    try:
        scale = math.pi / 180.
        lat1r = lat1*scale
        lon1r = lon1*scale
        lat2r = lat2*scale
        lon2r = lon2*scale

        a = math.sin(lat1r)*math.sin(lat2r)
        b = math.cos(lat1r)*math.cos(lat2r) * math.cos(lon2r-lon1r)
        apb = a + b
        if apb > 1.0:
            apb = 1.0  # avoid rounding error
        if apb < -1.0:
            apb = -1.0  # avoid rounding error
        c = math.acos(apb)
    except Exception, e:
        print ('longitudeLatitudeSep: <%s> <%s> <%s> <%s> -> %s %s -> %s %s'
               % (lon1, lat1, lon2, lat2, lat1r, lat2r, a, b))
        raise e
    R = 6378.  # radius of earth in km; in miles it's 3963.189
    return R*c


def main():
    """
    main
    """
    parser = optparse.OptionParser(usage="""
    %prog [-n NNearest][--category CAT] srcYamlDir dstYamlDir
    """)
    parser.add_option('-n', '--nnearest', action='store', type='int',
                      default=DEFAULT_N_NEAREST,
                      help="How many nearby locations to save")
    parser.add_option('--inv', action='store_true', default=False,
                      help="Store the inverse of the separation distance")
    parser.add_option('-o', '--out', action='store', default=DEFAULT_OUT_FNAME,
                      help="Output yaml file name")
    parser.add_option('--category', action='store', type='string',
                      help="Include only destinations belonging to this category (e.g. HOSPITAL)")
    opts, args = parser.parse_args()
    if len(args) != 2:
        parser.error('Two directory names are required')

    nToSave = opts.nnearest
    saveInverse = opts.inv
    outFName = opts.out
    parser.destroy()

    srcKeySet, srcRecs = yaml_tools.parse_all(args[0])  # @UnusedVariable
    dstKeySet, dstRecs = yaml_tools.parse_all(args[1])  # @UnusedVariable
    if opts.category:
        dstRecs = [rec for rec in dstRecs if rec['category'] == opts.category]

    outTbl = {}
    for sR in srcRecs:
        srcNm = sR['abbrev']
        if srcNm in outTbl:
            raise RuntimeError("Duplicate records for source %s" % srcNm)
        else:
            pairL = []
            srcLon = sR['longitude']
            srcLat = sR['latitude']
            for dR in dstRecs:
                dstNm = dR['abbrev']
                dstLon = dR['longitude']
                dstLat = dR['latitude']
                pairL.append((longitudeLatitudeSep(srcLon, srcLat, dstLon, dstLat), dstNm))
            pairL.sort()
            pairL = pairL[:nToSave]
            if saveInverse:
                pairL = [(1.0/sep, nm) for sep, nm in pairL]
            outTbl[srcNm] = {nm: wt for wt, nm in pairL}

    with open(outFName, 'w') as f:
        yaml.safe_dump(outTbl, f, default_flow_style=True, indent=4,
                       encoding='utf-8', width=130, explicit_start=True)

############
# Main hook
############

if __name__ == "__main__":
    main()
