import sys
import re
from collections import defaultdict

"""
example lines:

[0]INFO:infectionTracking:PatientAgent_ICU_HOSPITAL_Patch_0_0_COPL_2000_H_ICU_0_5 newly colonized at 1 in fac COPL_2000_H, tier 
ICU, ward 0

[0]INFO:infectionTracking:PatientAgent_HOSP_HOSPITAL_Patch_0_0_SAIN_2875_H_HOSP_2_42 arriving in community at time 1 with colonization status COLONIZED

[0]INFO:infectionTracking:PatientAgent_COMMUNITY_Patch_0_0_C2128_0_birth arriving in community at time 1 with colonization statu
s CLEAR

[0]INFO:infectionTracking:PatientAgent_HOME_COMMUNITY_Patch_0_0_C649_823 leaving community at time 1 with colonization status CLEAR

"""

reColonize = re.compile(r"INFO:infectionTracking:(\S+)\snewly colonized at\s(\d+)\sin fac\s(\S+),\stier\s(\S+),\sward\s(\d+)")
reDecolonize = re.compile(r"INFO:infectionTracking:(\S+)\sdecolonized at\s(\d+)\sin fac\s(\S+),\stier\s(\S+),\sward\s(\d+)")
reArrive = re.compile(r"INFO:infectionTracking:(\S+)\sarriving in community at time\s(\d+)\swith colonization status\s(\S+)")
reDepart = re.compile(r"INFO:infectionTracking:(\S+)\sleaving community at time\s(\d+)\swith colonization status\s(\S+)")

colonizations = []  # agent, time, fac, tier, ward
decolonizations = []  # agent, time, fac, tier, ward
arrivals = []       # agent, time, status
departures = []     # agent, time, status

agentHistory = defaultdict(list)

def appendHistory(event, data):
    agent = data[0]
    agentHistory[agent].append(((event,) + data[1:]))


def printHistory(history):
    for h in history:
        print h
    
def processHistory():
    colonizationCountHist = defaultdict(int)
    colonizationCountHistByYear = defaultdict(lambda:defaultdict(int))
    daysPerYear = 100
    maxYears = 30
    
    for agent,history in agentHistory.items():
        lastDepArr = None
        lastDepArrStatus = None
        colonizationCount = 0
        colonizationCountByYear = defaultdict(int)
        for hItem in history:
            event = hItem[0]
            date = int(hItem[1])
            if event=="comArrival" or event=="comDepart":
                if lastDepArr == event:
                    print "!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?"
                    print "arrivals/departures wrong for agent %s"%agent
                    printHistory(history)
                    print "!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?"
                    print "**********************************"
                lastDepArr = event
                
            if event == "colonization":
                if lastDepArr == "comArrival":
                    print "!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?"
                    print "colonization event while in community for agent %s"%agent
                    printHistory(history)
                    print "!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?"
                    print "**********************************"
                    
                colonizationCount +=1
                year = int(date / daysPerYear)
                colonizationCountByYear[year] += 1
                
        if 0 and colonizationCount > 3:
            print "**********************************"
            print "agent with high colonization count: %s"%agent
            printHistory(history)
            print "**********************************"
            
        colonizationCountHist[colonizationCount] += 1
        for i in xrange(maxYears):
            colonizationCountHistByYear[i][colonizationCountByYear[i]] += 1

    print "colonization count histogram:"
    print colonizationCountHist
    print "by year"
    for y in xrange(maxYears):
        print colonizationCountHistByYear[y]

def appendBuggyHistory():
    appendHistory('comArrival', ("buggyAgent1", 5, 'CLEAR'))
    appendHistory('comArrival', ("buggyAgent1", 10, 'CLEAR'))
    appendHistory('comDepart', ("buggyAgent2", 5, 'CLEAR'))
    appendHistory('comDepart', ("buggyAgent2", 10, 'CLEAR'))
    appendHistory('comArrival', ("buggyAgent3", 10, 'CLEAR'))
    appendHistory('colonization', ("buggyAgent3", 15, 'fac', 'HOSP', 0))
    
def parse(line):
    line = line.strip()
    m = reColonize.search(line)
    if (m):
        data = m.group(1,2,3,4,5)
        appendHistory('colonization', data)
        colonizations.append(data)
        return

    m = reDecolonize.search(line)
    if (m):
        data = m.group(1,2,3,4,5)
        appendHistory('decolonization', data)
        decolonizations.append(data)
        return

    m = reArrive.search(line)
    if (m):
        data = m.group(1,2,3)
        appendHistory('comArrival', data)
        arrivals.append(data)
        return

    m = reDepart.search(line)
    if (m):
        data = m.group(1,2,3)
        appendHistory('comDepart', data)
        departures.append(data)
        return

    raise RuntimeError("unreachable!")

    x,it,remainder = line.partition(':infectionTracking:')
    agent,sp,remainder = remainder.partition(' ')
    if remainder.startswith('newly colonized'):
        x,y,remainder = remainder.partition('newly colonized at ')
        time,x,loc = remainder.partition(' in fac ')
        

def writeCsvs():
    with open("infTrack_colonizations.csv", "w") as f:
        f.write("agent, time, fac, tier, ward\n")
        for e in colonizations:
            f.write("%s,%s,%s,%s,%s\n"%(e[0], e[1], e[2], e[3], e[4]))

    with open("infTrack_decolonizations.csv", "w") as f:
        f.write("agent, time, fac, tier, ward\n")
        for e in decolonizations:
            f.write("%s,%s,%s,%s,%s\n"%(e[0], e[1], e[2], e[3], e[4]))

    with open("infTrack_arrivals.csv", "w") as f:
        f.write("agent, time, status\n")
        for e in arrivals:
            f.write("%s,%s,%s\n"%(e[0], e[1], e[2]))

    with open("infTrack_departures.csv", "w") as f:
        f.write("agent, time, status\n")
        for e in departures:
            f.write("%s,%s,%s\n"%(e[0], e[1], e[2]))

            
def main():
    inputLog = sys.argv[1]

    with open(inputLog) as f:
        for line in f:
            if "INFO:infectionTracking" in line:
                parse(line)

    if 0:
        appendBuggyHistory()

    if 0:
        writeCsvs()

    if 1:
        processHistory()
        

main()
