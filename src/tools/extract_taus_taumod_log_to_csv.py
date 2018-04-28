import re
import cPickle as pickle
from collections import defaultdict

numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
rx = re.compile(numeric_const_pattern, re.VERBOSE)

def floatExtract(starts, s):
    k,v = s.split(' ')
    if k != starts:
        raise RuntimeException("%s doesn't match %s"%(k, starts))
    return float(rx.match(v).group(0))

vsnfExceptions = ["GLEN_270_S", "PINN_2222_S", "SALE_1314_S", "THI_5400_S"]
catList = ["COMMUNITY",
           "SNF",
           "VSNF",
           "LTACH",
           "HOSPITAL"]

def whichCat(fac):
    if fac.endswith("_C"):
        return "COMMUNITY"
    elif fac.endswith("_S"):
        if fac in vsnfExceptions:
            return "VSNF"
        else:
            return "SNF"
    elif fac.endswith("_V"):
        return "VSNF"
    elif fac.endswith("_L"):
        return "LTACH"
    elif fac.endswith("_H"):
        return "HOSPITAL"
    raise RuntimeError("invalid facility name")



def readTaumodLog(fname, outname):
    tauDict = {}
    instances = defaultdict(int)

    with open(fname) as f:
        with open(outname, "w") as o:
            o.write("fac,tier,category,day,prevalence,expected,ratio,tau,newTau\n")
            for line in f.readlines():
                try:
                    facInfo, data = line.split(':')
                    data = data.strip()
                    if not data.startswith('prevalence'):
                        continue

                    fac, tier = facInfo.split(' ')
                    prevalence, expected, ratio, tau, newTau = data.split(', ')
                    prevalence = floatExtract('prevalence', prevalence)
                    expected = floatExtract('expected', expected)
                    ratio = floatExtract('ratio', ratio)
                    tau = floatExtract('tau', tau)
                    newTau = floatExtract('newTau', newTau)
                    instances[(fac,tier)] += 50

                    print "%s, %s: %s %s %s %s %s"%(fac, tier, prevalence, expected, ratio, tau, newTau)
                except:
                    continue

                o.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%(fac, tier, whichCat(fac), instances[(fac,tier)], prevalence,
                                                     expected, ratio, tau, newTau))
                
                
                
            

def main():
    readTaumodLog("/home/jim/mnt/olympus/WORK/pyrhea/pyrhea/src/sim/taumod.out", "tauInfo.csv")



main()
