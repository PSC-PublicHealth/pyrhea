import re
import cPickle as pickle

numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
rx = re.compile(numeric_const_pattern, re.VERBOSE)

def floatExtract(starts, s):
    k,v = s.split(' ')
    if k != starts:
        raise RuntimeException("%s doesn't match %s"%(k, starts))
    return float(rx.match(v).group(0))


def readTaumodLog(fname):
    tauDict = {}
    
    with open(fname) as f:
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

                print "%s, %s: %s %s %s %s %s"%(fac, tier, prevalence, expected, ratio, tau, newTau)
            except:
                pass

            if newTau > 0.9:
                newTau = 0.9
            tauDict[(fac, tier)] = newTau
            
    return tauDict

def main():
    tauDict = readTaumodLog("/home/jim/mnt/olympus/WORK/pyrhea/pyrhea/src/sim/taumod.out")

    with open("new_taus.pkl", "w") as f:
        pickle.dump(tauDict, f)


main()
