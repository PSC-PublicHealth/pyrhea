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

_rhea_svn_id_ = "$Id: pyrhea.py 30 2015-06-01 21:48:41Z welling $"

import sys
import math
import set_paths
from random import randint, random, shuffle, seed, choice
import patches
import pyrheabase
from pyrheautils import enum, namedtuple
import random

DiseaseStatus = enum('Susceptible',
                    'Exposed',
                    'Infectious',
                    'Asymptomatic',
                    'Recovered')

StatusList = ['Susceptible',
              'Exposed',
              'Infectious',
              'Asymptomatic',
              'Recovered']

PersonStatus = namedtuple('PersonStatus',
                          ['diseaseStatus',
                           'diseaseStatusCountDown'],
                          field_types=[DiseaseStatus,None])

R0 = 1.4
InfectiousPeriod = 4.1
LatentPeriod = 1.9
Beta = R0/(InfectiousPeriod)
AsymptProb = 0.33


class Community(patches.MultiInteractant):
    def __init__(self,name,patch,nAgents):
        patches.MultiInteractant.__init__(self,name,nAgents,patch)
        self.checkInterval = 1 # check the population once per day
        self.populationStatus = [0.0 for id,key in DiseaseStatus.names.items()]
        self.totalAgents = nAgents
        self.dailyInfectiousAttemps = []
        
    def printPopulationStatus(self,mainLoop,timeNow):
        returnString = "Community Status for Time {0}: ".format(timeNow)
        for i in xrange(len(self.populationStatus)):
            returnString += " {0}:{1} ".format(StatusList[i],self.populationStatus[i])
        
        print returnString
            
class Person(patches.Agent):
    def __init__(self,name,patch,timeNow=0,comIndex=0,debug=False):
        patches.Agent.__init__(self,name,patch,debug=debug)
        self._status = PersonStatus(DiseaseStatus.Susceptible,0)
        self.lastUpdateTime = timeNow
        self.comIndex = comIndex
    
    def setDiseaseStatus(self,status,decrement_=True):
        currentStatus =self._status.diseaseStatus
        newCountDown = 0;
        if status == 'Susceptible':
            newCountDown = 1;
        elif status == 'Exposed':
            newCountDown = LatentPeriod
        elif status == 'Infectious' or status == 'Asymptomatic':
            newCountDown = InfectiousPeriod
        elif status == "Recovered":
            newCountDown = 0;
        
        self._status = PersonStatus(getattr(DiseaseStatus,status),newCountDown)
                            
        if hasattr(self,'community'):
            #print "incrementing {0} by 1".format(self._status.diseaseStatus)
            if decrement_:
                self.community.populationStatus[currentStatus] -=1
            self.community.populationStatus[self._status.diseaseStatus] += 1
    
    def updateIfInfected(self,timeNow):
        if self._status.diseaseStatus == DiseaseStatus.Susceptible:
            if self.comIndex in self.community.dailyInfectiousAttemps:
                self._status = PersonStatus(DiseaseStatus.Exposed,LatentPeriod)
                self.community.populationStatus[DiseaseStatus.Susceptible] -=1
                self.community.populationStatus[DiseaseStatus.Exposed] +=1
        if self.comIndex in self.community.dailyInfectiousAttemps:
            del self.community.dailyInfectiousAttemps[self.community.dailyInfectiousAttemps.index(self.comIndex)]
        
    def updateDiseaseStatus(self,timeNow):
        dT = timeNow - self.lastUpdateTime
        #print "TimeNow = {0}".format(timeNow)
        if dT > 0:
            self.lastUpdateTime = timeNow
            ## Decrement the state counter
            newCountDown = self._status.diseaseStatusCountDown - dT
            newDiseaseStatus = self._status.diseaseStatus
            #if self._status.diseaseStatus == DiseaseStatus.Exposed:
                #print "New CountDown for Exposed: {0}".format(newCountDown)
            if newCountDown <= 0:
                if self._status.diseaseStatus == DiseaseStatus.Exposed:
                    #print "Updating Exposed Status: {0}".format(newCountDown)
                    if random.random() < AsymptProb:
                         newDiseaseStatus = DiseaseStatus.Asymptomatic
                    else:
                        newDiseaseStatus = DiseaseStatus.Infectious
                    timeDiff = newCountDown
                    newCountDown = InfectiousPeriod - timeDiff
                    #print "Exposed new countdown: {0}".format(newCountDown)
                elif self._status.diseaseStatus == DiseaseStatus.Infectious or self._status.diseaseStatus == DiseaseStatus.Asymptomatic:
                    newDiseaseStatus = DiseaseStatus.Recovered
                    newCountDown = 0
                elif self._status.diseaseStatus == DiseaseStatus.Susceptible:
                    newCountDown = 1.0
                elif self._status.diseaseStatus == DiseaseStatus.Recovered:
                    newCountDown = 0.0
                    pass
            
            if self._status.diseaseStatus == DiseaseStatus.Infectious or self._status.diseaseStatus == DiseaseStatus.Asymptomatic:
                ## determine the number of attempts to infect
                factor=1.0
                if self._status.diseaseStatus == DiseaseStatus.Asymptomatic:
                    factor = 0.5
                
                nAttempts = int(math.floor(Beta*factor))
                if random.random() < (Beta*factor % 1.0):
                    nAttempts += 1
                
                for i in xrange(nAttempts):
                    attIndex = random.randint(0,self.community.totalAgents)
                    if attIndex not in self.community.dailyInfectiousAttemps:
                        self.community.dailyInfectiousAttemps.append(attIndex)
                                
            self.community.populationStatus[self._status.diseaseStatus] -=1
            self.community.populationStatus[newDiseaseStatus] +=1
            self._status = PersonStatus(newDiseaseStatus,newCountDown)
            
            ## must check if this agent will get infected
                    ## get the total number of infectious people in this bin
    
    def run(self,startTime):
        timeNow =startTime 
        while True:
            #timeNow = self.community.lock(self)
            self.updateIfInfected(timeNow)
            self.updateDiseaseStatus(timeNow)
            timeNow = self.sleep(1)
            if timeNow > 200:
                sys.exit()
            
        

    def __getstate__(self):
        d = patches.Agent.__getstate__(self)
        d['status'] = self.status
        d['lastUpdateTime'] = self.lastUpdateTime
        
    def __setstate__(self,stateDict):
        patches.Agent.__setstate(self,stateDict)
        self.status = stateDict['status']
        self.lastUpdateTime = stateDict['lastUpdateTime']
        


def main():    

    trace = False
    verbose = False # @UnusedVariable
    debug = False
    deterministic = False
    
    for a in sys.argv[1:]:
        if a == '-v':
            verbose = True  # @UnusedVariable
        elif a == '-d':
            debug = True
        elif a == '-t':
            trace = True
        elif a == '-D':
            deterministic = True
        else:
            describeSelf()
            sys.exit('unrecognized argument %s' % a)
        
    comm = patches.getCommWorld()
    
    patchGroup = patches.PatchGroup(comm,trace=trace,deterministic=False)
    
    nPatches = 1
    nPeoplePerPatch = 10000
    nSeeds = int(nPeoplePerPatch*.01)
    
    for j in xrange(nPatches):
        patch = patchGroup.addPatch(patches.Patch(patchGroup))
        community = Community("Pittsburgh",patch,nPeoplePerPatch)
        patch.loop.addPerDayCallback(community.printPopulationStatus)
        allAgents = []
        for i in xrange(nPeoplePerPatch):
            person = Person('Person_{0}_{1}'.format(str(patch.tag),i),
                            patch,comIndex=i,debug=debug)
            community.lock(person)
            person.community = community
            #community.unlock(person)
            person.setDiseaseStatus('Susceptible',False)
            allAgents.append(person)
        
        seedAgents = []
        for i in xrange(nSeeds):
            if len(seedAgents)>0:
                seedIndex = seedAgents[-1]
                while seedIndex in seedAgents:
                    seedIndex = random.randint(0,nPeoplePerPatch-1)
                seedAgents.append(seedIndex)
            else:
                seedAgents.append(random.randint(0,nPeoplePerPatch-1))
        
        for iS in seedAgents:
            allAgents[iS].setDiseaseStatus('Exposed')
            #allAgents[iS]._status.diseaseStatusCountDown = LatentPeriod
        
        for agent in allAgents:
            print agent.name
            #community.populationStatus[agent._status.diseaseStatus] += 1.0
        
        print community.populationStatus
        patch.addAgents(allAgents)
        
    print allAgents[0].name
    print "Done Initializing"
    patchGroup.start()
    print "we are done!!!!"
                   
############
# Main hook
############

if __name__ == "__main__":
    main()       
            
    
                        
                
                    

