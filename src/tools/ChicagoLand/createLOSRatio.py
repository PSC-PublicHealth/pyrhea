import csv
import yaml
import os.path

def readDictFromCsv(fileName, key, val):
    ret = {}
    with open(fileName) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ret[row[key]] = row[val]

    return ret

data_1023 = readDictFromCsv("../../../models/ChicagoLand/Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016_BedLOS.csv", "UNIQUE_ID", "LOS")
data_0929 = readDictFromCsv("../../../models/ChicagoLand/Matrices_LOS_09292016_Facilities_LOS.csv", "UNIQUE_ID", "Mean")

ratios = {}
for site, los in data_1023.items():
    try:
        los929 = data_0929[site]
        ratios[site] = float(los) / float(los929)

        if site.endswith('L'):
            print "%s: los1023: %s, los929: %s, ratio: %s"%(site, los, los929, ratios[site])
        
    except:
        pass


facilityfactsDir = "../../../models/ChicagoLand/facilityfacts"

provenance = "Geocoded_LOSCMSAHQ_2010Cohort_LOS_Bedsize_10-23-2016_BedLOS.csv mean LOS vs Matrices_LOS_09292016_Facilities_LOS.csv LOS"
    
for site, ratio in ratios.items():
    siteFName = "%s.yaml"%site
    fname = os.path.join(facilityfactsDir, siteFName)
    try:
        with open(fname) as f:
            yData = yaml.load(f)
    except:
        print "Can't find file %s"%fname
        continue

    yData['scaleLengthOfStay'] = {'prov': provenance, 'value': ratio}

#    print yData
    with open(fname, 'w') as f:
        yaml.dump(yData, f)


