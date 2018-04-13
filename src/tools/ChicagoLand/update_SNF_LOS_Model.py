import csv
import yaml
import os.path
from phacsl.utils.misc.util import isiterable

def readDictFromCsv(fileName, key, val):
    ret = {}
    with open(fileName) as f:
        reader = csv.DictReader(f)
        for row in reader:
            k = row[key]
            if isiterable(val):
                try:
#                    ret[row[key]] = [row[h] for h in val]
                    ret[row[key]] = [float(row[h]) for h in val]
                except:
                    ret[row[key]] = [.800089, 3.1730180781906, 1.0497094393532238, 0.0002451]
            else:
                ret[row[key]] = row[val]

    return ret

fitParms = readDictFromCsv("../../../models/ChicagoLand/los_model_fit_parms_fixedmean.csv", "abbrev", ["k", "mu", "sigma", "lmda"])
#provenance = "Some BS Joel cooked up"
provenance = "los_model_fits_fixedmean.py f9405f1d applied to Histogram_09292016.csv"

facilityfactsDir = "../../../models/ChicagoLand/facilityfacts"
facilityfactsUpdateDir = "../../../models/ChicagoLand/facilityfacts_new"
for site, parms in fitParms.items():
    if not site.endswith("_S"):
        continue
    
    siteFName = "%s.yaml"%site
    fname = os.path.join(facilityfactsDir, siteFName)
    updateFname = os.path.join(facilityfactsUpdateDir, siteFName)
    try:
        with open(fname) as f:
            yData = yaml.load(f)
    except:
        print "Can't find file %s"%fname
        continue

    yData['losModel']['parms'] = parms
    yData['losModel']['prov'] = provenance

#    print yData
    with open(updateFname, 'w') as f:
        yaml.dump(yData, f)


