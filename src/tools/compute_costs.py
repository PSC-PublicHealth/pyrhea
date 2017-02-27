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

import os.path
import sys
import math
import yaml
import glob
import random
from optparse import OptionParser
from multiprocessing import Process,Manager

import numpy as np
import scipy.stats as st

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))
import pyrheautils
import schemautils
import pathogenbase as pth
from facilitybase import CareTier
from notes_plotter import readFacFiles, checkInputFileSchema
from notes_plotter import SCHEMA_DIR, INPUT_SCHEMA
from time_series_plotter import mergeNotesFiles, getTimeSeriesList
import print_xdro_counts

DEFAULT_OUT_FILE = 'costs_output.yaml'

#targetYear = 2017
#incidentCarriers = 393.41
#numXDROReg = 393.41
#numCREBund = 9682.83
#totCarriersGettingCRE = 17.95
#carrierCPdays = 11.21
#nonCarrierCPDays = 3350.67
fractionAttribMorts = [0.35]
infectionGivenCol = 0.25

class Constants:
    daysPerYear = 365.0
    minutesPerHour = 60.0
    workDayHours = 8.0
class Cost:
    def __init__(self,value_,year_):
        self.value = value_
        self.year = year_
        
    def discountedValue(self,targetYear,discountRate=0.03):
        yearDiff = targetYear - self.year
        if yearDiff == 0:
            return Cost(self.value, self.year)
        else:
            return Cost(self.value * (1.0 + discountRate)**yearDiff,targetYear)
    
    def valueInYear(self,targetYear,discountRate=0.03):
        yearDiff = targetYear - self.year
        if yearDiff == 0:
            return Cost(self.value, self.year)
        else:
            return Cost(self.value / (1.0 + discountRate)**yearDiff,targetYear)
    def __str__(self):
        return "{0} USD {1}".format(self.value, self.year)

    def __iadd__(self,other):
        self.value += other.discountedValue(self.year).value
        return self
    
    def __mul__(self,other):
        if isinstance(other,float):
            return Cost(self.value*other,self.year)
        else:
            raise RuntimeError("Doesn't make sense")
        
            
    
    
class Distribution(object):
    def __init__(self,type_,args_):
        self.type = type_
        self.args = args_
        
    def draw(self):
        if self.type == 'value':
            return self.args['value']
        elif self.type == 'gamma':
            mean = self.args['mean']
            stdev = self.args['stdev']
            
            alpha = (mean/stdev)**2
            beta = (stdev**2)/mean
            
            return np.random.gamma(alpha,beta)
        
        elif self.type == 'beta':
            mean = self.args['mean']
            stdev = self.args['stdev']
            
            alpha = ((1.0 - mean)/(stdev**2)-(1.0/mean))*(mean**2)
            beta = alpha*((1.0/mean)-1.0)
            
            return np.random.beta(alpha,beta)
        
        elif self.type == 'uniform':
            return np.random.uniform(self.args['low'],self.args['high'])
        
        else:
            raise RuntimeError("unknown distribution: {0}".fromat(self.type))
        

def calculateNPVValues(startAge_,annualWage_,targetYear_, params_, discountRate_=0.03):
    
    annualWageDiscounted = annualWage_.discountedValue(targetYear_)

    
    lifeExpAtAge = int(round(params_['lifetimeProdLossParams']['expectancyAtAge'][startAge_]))
    
    #print lifeExpAtAge-startAge_
    npvProdLoss = Cost(sum([annualWageDiscounted.valueInYear(targetYear_ + x).value for x in range(0,lifeExpAtAge)]),targetYear_)
    
    #print npvProdLoss
    
    qWeight = 1.0
    for i in params_['lifetimeProdLossParams']['npvQALYWeight']['startAges']:
        #print "{0}: {1}".format(startAge_,i)
        if i > startAge_:
            break
        qWeight = params_['lifetimeProdLossParams']['npvQALYWeight']['weights'][params_['lifetimeProdLossParams']['npvQALYWeight']['startAges'].index(i)]
    
    
    npvQ = 0.0
    for i in range(0,lifeExpAtAge):
        npvQ += Cost(qWeight,targetYear_).valueInYear(targetYear_+i).value
    
    npvQCost = Cost(npvQ,targetYear_)
    
    return npvProdLoss, npvQCost, Cost(qWeight,targetYear_)
 
def computeCostOfXDROReg(numPatients_,targetYear_, params_):
    
    tRegistryLogin = params_['interventionParameters']['personnelTimeXDROLoginMinutes']
    tRegistryPerPatient = params_['interventionParameters']['personnelTimeXDROSearchMinutesPerPatient']
    
    cInfectControlWageDict = params_['interventionParameters']['personnelXDROWage']
    cInfectControlWage = Cost(Distribution(cInfectControlWageDict['cost']['distribution']['type'],
                                           cInfectControlWageDict['cost']['distribution']['args']).draw(),
                              cInfectControlWageDict['cost']['year'])
    
    value = ((tRegistryLogin*Constants.daysPerYear) + (numPatients_*tRegistryPerPatient))\
          * (cInfectControlWage.discountedValue(targetYear_).value/Constants.minutesPerHour)

    return Cost(value,targetYear_)

def getNContactsPerDayByType(locationType_, params_):
    contactRates = params_['contacts']
    if locationType_ in ['LTACH','HOSPITAL']:
        locationType_ = 'GeneralWards'
    ### Update these with RHEA monikers, may need a translate
    return Distribution(contactRates['perDayIn{0}'.format(locationType_)]['distribution']['type'],
                        contactRates['perDayIn{0}'.format(locationType_)]['distribution']['args']).draw()
                        
def computeCostsOfContactPrecautions(contactPrecautionDays_, locationType_, rNurseWage_, cGloves_, cGowns_, targetYear_, params_):
    
    contactRate = getNContactsPerDayByType(locationType_,params_)
    tContactPrecautionDict = params_['interventionParameters']['timeDonDoffMinutes']
    tCPMinutes = Distribution(tContactPrecautionDict['distribution']['type'],
                              tContactPrecautionDict['distribution']['args']).draw()
    
    value = (contactPrecautionDays_ \
          * contactRate) \
          * ((rNurseWage_.discountedValue(targetYear_).value/60.0) \
          * tCPMinutes \
          + cGloves_.discountedValue(targetYear_).value \
          + cGowns_.discountedValue(targetYear_).value)
    
    return Cost(value,targetYear_)

def computeCostsOfBundles(numBundles_,numBaths_, targetYear_, params_):
    chgWipeDict = params_['interventionParameters']['chgWipesPerBath']['cost']
    cChgWipesPerBath = Cost(Distribution(chgWipeDict['distribution']['type'],
                                        chgWipeDict['distribution']['args']).draw(),
                           chgWipeDict['year'])
    
    screenDict = params_['interventionParameters']['screeningTotal']['cost']
    cScreening = Cost(Distribution(screenDict['distribution']['type'],
                                   screenDict['distribution']['args']).draw(),
                           screenDict['year'])
    
    #print cChgWipesPerBath.discountedValue(targetYear_)
    #print cScreening.discountedValue(targetYear_)
    value = numBundles_ * cScreening.discountedValue(targetYear_).value \
          + numBaths_ * cChgWipesPerBath.discountedValue(targetYear_).value
    
    return Cost(value, targetYear_)      
        
                                        
def determineOutcomes(nIncidence_,nCarriersCRE_, probInfect_, attribMort_, cGloves_, cGowns_, rNurseWage_, personAge_, annualWage_, hourlyWage_, targetYear_, params_):
    outcomesDict = params_['outcomes']
    outcomeCosts = {x:{} for x in outcomesDict.keys()}
    
    #nOutsDict = {x:Distribution(outcomesDict[x]['probability']['distribution']['type'],
    #                            outcomesDict[x]['probability']['distribution']['args']).draw() \
    #             * (nIncidence_+nCarriersCRE_) \
    #             * probInfect_ \
    #             for x in outcomesDict.keys()}
    therapyDict = params_['therapy']
    
#   npvInf = []
#     npQInf = []
#     baseQALYInf = []
#     for i in range(0,int(nIncidence)):
#         startAge =  int(round(Distribution('uniform',{'low':60.0,'high':80.0}).draw()))
#         npv,npQ,baseQALY = calculateNPVValues(startAge,  annualWage_, params_)
#         npvInf.append(npv)
#         npQInf.append(npQ)
#         baseQALYInf.append(baseQALY)  
#     
    ### Get the cost of tests
    testsDict = params_['tests']
    testsCosts = {}
    testsNums = {}
    for t,tD in testsDict.items():
        testsCosts[t] = Cost(Distribution(tD['cost']['distribution']['type'],
                            tD['cost']['distribution']['args']).draw(),
                     tD['cost']['year'])
        
        testsNums[t] = {}
        nTests = tD['numberOfTests']
        for nt,ntD in nTests.items():
            testsNums[t][nt] = Distribution(ntD['distribution']['type'],
                                            ntD['distribution']['args']).draw()
            
            
    ### Gets the cost of drugs 
    drugsDict = params_['drugs']
    drugCosts = {}
    for d,dD in drugsDict.items():
        drugCosts[d] = Cost(Distribution(dD['cost']['distribution']['type'],
                                         dD['cost']['distribution']['args']).draw(),
                            dD['cost']['year'])
    
                          
        
    for k,v in outcomesDict.items():
        #if k != 'nonVentilatorAssociatedPneumonia':
        #    continue
        if k not in ['pneumoniaAll']:
            nCases = Distribution(v['probability']['distribution']['type'],
                                  v['probability']['distribution']['args']).draw() \
                   * (nIncidence_+nCarriersCRE_) * probInfect_
            outcomeCosts[k]['nCases'] = nCases
            #print "{0}: {1}".format(k,nCases)
            mortProbDict = v['probabilityOfMortality']
            nDead = 0.0
            therapyProbs = {}
            therapyDrugReg = {}
            for mk,mD in mortProbDict.items():
                if mk == 'combination':
                    continue
                therapyProbs[mk] = Distribution(therapyDict[mk]['prob']['distribution']['type'],
                                                therapyDict[mk]['prob']['distribution']['args']).draw()
                
                therapyDrugReg[mk] = {}
                for d,dD in therapyDict[mk]['drugs'].items():
                    if d.find('_and_') > -1:
                        drugs = d.split("_and_")
                    else:
                        drugs = [d]
                    
                    for dr in drugs:
                        therapyDrugReg[mk][dr] = dD
                
            
                #print type(attribMort_)
                #print type(therapyProbs[mk])
                #print type(nCases)
                nDead += Distribution(mD['distribution']['type'],
                                      mD['distribution']['args']).draw() \
                       * nCases * attribMort_ * therapyProbs[mk]
            
                
                drugCostsForTherapy = {}
                for th,thD in therapyDrugReg.items():
                    drugCostsForTherapy[th] = 0.0
                    for d,dD in thD.items():
                        drugCost = drugCosts[d]
                        drugCostsForTherapy[th] += drugCost.discountedValue(targetYear_).value * dD
                
            outcomeCosts[k]['nDead'] = nDead    
            
            ## NPV and QALYs for live
            
            npvCases = Cost(0.0,targetYear_)
            npQCases = Cost(0.0,targetYear_)
            baseQALYCases = Cost(0.0,targetYear_)
            npvDead = Cost(0.0,targetYear_)
            npQDead = Cost(0.0,targetYear_)
            baseQALYDead = Cost(0.0,targetYear_)
            
            nCasesInt = math.floor(nCases)
            nCasesRem = nCases - nCasesInt
            nDeadInt = math.floor(nDead)
            nDeadRem = nDead - nDeadInt
            
            for i in range(0,int(nCasesInt)):
                startAge = int(round(Distribution('uniform',{'low':60.0,'high':80.0}).draw()))
                npv,npQ,baseQALY = calculateNPVValues(startAge,  annualWage_, targetYear_, params_)
                if i < nDeadInt:
                    #print npv
                    #print npvDead
                    npvDead += npv
                    npQDead += npQ
                    baseQALYDead += baseQALY
                npvCases += npv
                npQCases += npQ
                baseQALYCases += baseQALY
                
            #REMANDER
            startAge = int(round(Distribution('uniform',{'low':60.0,'high':80.0}).draw()))
            npv,npQ,baseQALY = calculateNPVValues(startAge,  annualWage_, targetYear_,  params_)
            npvCases +=Cost( nCasesRem*npv.value,npv.year)
            npQCases += Cost(nCasesRem*npQ.value,npQ.year)
            baseQALYCases += Cost(nCasesRem*baseQALY.value,baseQALY.year)
            npvDead += Cost(nDeadRem*npv.value,npv.year)
            npQDead += Cost(nDeadRem*npQ.value,npQ.year)
            baseQALYDead +=Cost( nDeadRem*baseQALY.value,baseQALY.year)
            
            
                
            #print drugCostsForTherapy
            #print "{0} Dead: {1}".format(k,nDead)
            ### Now ... hospital costs
            
            ### Gather Distributions Needed
            cICUBedDayDict = params_['costs']['icuBedDay']
            cICUBedDay = Cost(Distribution(cICUBedDayDict['distribution']['type'],
                                           cICUBedDayDict['distribution']['args']).draw(),
                              cICUBedDayDict['year'])
            
            cGenWardBedDayDict = params_['costs']['generalWardBedDay']
            cGenWardBedDay = Cost(Distribution(cGenWardBedDayDict['distribution']['type'],
                                               cGenWardBedDayDict['distribution']['args']).draw(),
                                  cGenWardBedDayDict['year'])
            
            nContactsICUDict = params_['contacts']['perDayInICU']
            nContactsICU = Distribution(nContactsICUDict['distribution']['type'],
                                        nContactsICUDict['distribution']['args']).draw()
            
            nContactsGenWardDict = params_['contacts']['perDayInGeneralWards']
            nContactsGenWard = Distribution(nContactsGenWardDict['distribution']['type'],
                                            nContactsGenWardDict['distribution']['args']).draw()
                                            
            attrLOS = Distribution(v['attributableLOS']['distribution']['type'],
                                   v['attributableLOS']['distribution']['args']).draw()
            
            probICU = Distribution(params_['probabilities']['patientInICU']['distribution']['type'],
                                   params_['probabilities']['patientInICU']['distribution']['args']).draw()
            
            hospOutcomesDict = v['hospitalizationOutcomes']['cost']
                                             
            hospOutcomes = Cost(Distribution(hospOutcomesDict['distribution']['type'],
                                             hospOutcomesDict['distribution']['args']).draw(),
                                hospOutcomesDict['year'])
                
            treatDurDays = Distribution(v['treatmentDuration']['distribution']['type'],
                                        v['treatmentDuration']['distribution']['args']).draw()
                                        
#             cICUBedDay = Cost(4750.316316,2015)
#             cGenWardBedDay = Cost(2632.7543,2013)
#             probICU = 0.442002305
#             nContactsICU = 116
#             nContactsGenWard = 68.5
#             attrLOS = 9.5
#             cGloves_ = Cost(0.102123536, 2015)
#             cGowns_ = Cost(0.75,2004)
#             rNurseWage_ = Cost(33.71405953,2014)
#             
            nICU = (nCases*probICU) 
            nGenW = (1.0 - probICU) * nCases
#             print nICU
#             print "nGenW = {0}".format(nGenW)
#             print "LOS = {0}".format(attrLOS)
            ggr = cGloves_.discountedValue(targetYear_).value \
                + cGowns_.discountedValue(targetYear_).value \
                + (rNurseWage_.discountedValue(targetYear_).value / Constants.minutesPerHour)
#             print ggr  
            cggrLOSICU = (cICUBedDay.discountedValue(targetYear_).value + (nContactsICU*ggr))*attrLOS
            cggrLOSGW =  (cGenWardBedDay.discountedValue(targetYear_).value + (nContactsGenWard*ggr))*attrLOS
#             print "CICU = {0}".format(nContactsICU*ggr)
#             print "c2ICU = {0}".format((cICUBedDay.discountedValue(targetYear_).value + (nContactsICU*ggr)))
#             print "c3ICU = {0}".format(nICU*cggrLOSICU)
#             print "cG = {0}".format((nContactsGenWard*ggr))
#             print "cG22={0}".format((cGenWardBedDay.discountedValue(targetYear_).value + (nContactsGenWard*ggr)))
#             print "c2G = {0}".format((cGenWardBedDay.discountedValue(targetYear_).value + (nContactsGenWard*ggr))*attrLOS)
#             print "c3G = {0}".format(nGenW*cggrLOSGW)
#             print cggrLOSICU
            hospCost= Cost(nICU*cggrLOSICU + nGenW*cggrLOSGW,targetYear_)       
            
            outcomeCosts[k]['Hospitalization Cost'] = hospCost
            #print "{0} Hosp = {1}".format(k,hospCost)
            
            utilityWeight = Distribution(v['utilityWeight']['distribution']['type'],
                                         v['utilityWeight']['distribution']['args']).draw()
            QALYLost = Cost(((utilityWeight*baseQALYCases.discountedValue(targetYear_).value)/Constants.daysPerYear)*attrLOS \
                            +npQDead.discountedValue(targetYear_).value,
                            targetYear_)
            
            outcomeCosts[k]['QALYs Lost'] = QALYLost                
            #print "QLost = {0}".format(QALYLost)
            
            ### prod Losses
            prodLosses = Cost(nCases * attrLOS * hourlyWage_.discountedValue(targetYear_).value * Constants.workDayHours,
                              targetYear_)
            
            outcomeCosts[k]['Productivity Lost'] = prodLosses
            #print "ProdLost {0}: {1}".format(k,prodLosses)
            
            #### Losses Mortality
            
            lossesMortality = Cost(npvDead.value,
                                   targetYear_)
            
            outcomeCosts[k]['Losses Due to Mortality'] = lossesMortality
            
            #print "MortLoss {0}: {1}".format(k,lossesMortality)
            
            ### LETS Build the first term of 3rd Party
            thirdPartyCosts = 0.0
            if k == "bacteremia":
                ### doing this to match the complex Excel formula, can do something more sane later
                cHospBacteremia = hospOutcomes.discountedValue(targetYear_)
                cPICC = testsCosts['PICCLine'].discountedValue(targetYear_)
                cBloodCulture = testsCosts['bloodCulture'].discountedValue(targetYear_)
                pICU = probICU
                nContactsDayICU = nContactsICU
                nContactsDay = nContactsGenWard
                LOS_Bacteremia = attrLOS
                pComboWithCarb = therapyProbs['carbapenemContaining']
                pComboNoCarb = therapyProbs['carbapenemSparing']
                pMonotherapy = therapyProbs['monotherapy']
                ## Already Discounted
                cComboCarbContainingDaily = Cost(drugCostsForTherapy['carbapenemContaining'],targetYear_)
                cComboCarbSparingDaily = Cost(drugCostsForTherapy['carbapenemSparing'],targetYear_)
                cMonotherapyDaily = Cost(drugCostsForTherapy['monotherapy'],targetYear_)
                ### From Spreadsheet cells L-M 62
                cCarbContainInitialDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                cCarbSparingInitialDose = Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                cMonotherapyInitalDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                
                tBacteremiaTreatment = treatDurDays
                
                ctCC = pComboWithCarb*((cComboCarbContainingDaily.value*tBacteremiaTreatment) + cCarbContainInitialDose.value)
                ctCS = pComboNoCarb*((cComboCarbSparingDaily.value*tBacteremiaTreatment) + cCarbSparingInitialDose.value)
                ctM  = pMonotherapy*((cMonotherapyDaily.value*tBacteremiaTreatment) + cMonotherapyInitalDose.value)
                cBC = cBloodCulture.value * testsNums['bloodCulture']['bacteremia']
                ggr = (cGloves_.discountedValue(targetYear_).value \
                    + cGowns_.discountedValue(targetYear_).value \
                    + (rNurseWage_.discountedValue(targetYear_).value / Constants.minutesPerHour))
            
                thirdPartyCosts = Cost(nCases*(cHospBacteremia.value + cPICC.value + ctCC + ctCS + ctM + cBC) \
                                        + ((pICU*nCases*nContactsDayICU \
                                        + (1.0-pICU)*nCases*nContactsDay) \
                                        * ggr*LOS_Bacteremia),
                                       targetYear_)
                
            elif k == 'UTI':
                #nCases = 29.673531
#                 cHospUTI = Cost(6900.72812,2012).discountedValue(targetYear_)#hospOutcomes.discountedValue(targetYear_)
#                 cPICC = testsCosts['PICCLine'].discountedValue(targetYear_)
#                 cUA = Cost(2.933778,2015).discountedValue(targetYear_) #testsCosts['UA'].discountedValue(targetYear_)
#                 cUrineCulture = Cost(10.6956948,2015).discountedValue(targetYear_)
#                 pICU = probICU
#                 nContactsDayICU = nContactsICU
#                 nContactsDay = nContactsGenWard
#                 LOS_UTI = 6.0 #attrLOS
#                 pComboWithCarb = therapyProbs['carbapenemContaining']
#                 pComboNoCarb = therapyProbs['carbapenemSparing']
#                 pMonotherapy = therapyProbs['monotherapy']
#                 ## Already Discounted
#                 cComboCarbContainingDaily = Cost(273.82678,targetYear_)#Cost(drugCostsForTherapy['carbapenemContaining'],targetYear_)
#                 cComboCarbSparingDaily = Cost(393.5630219,targetYear_)#Cost(drugCostsForTherapy['carbapenemSparing'],targetYear_)
#                 cMonotherapyDaily = Cost(136.9133886,targetYear_)#Cost(drugCostsForTherapy['monotherapy'],targetYear_)
#                 ### From Spreadsheet cells L-M 62
#                 cCarbContainInitialDose = Cost(0,targetYear_)#Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
#                                               #targetYear_)
#                 cCarbSparingInitialDose = Cost(303.4237654,targetYear_)#Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
#                                           #    targetYear_)
#                 cMonotherapyInitalDose = Cost(0,targetYear_)#Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
#                                          #     targetYear_)
#                 
#                 tUTITreatment = 12.0 #treatDurDays
                cHospUTI = hospOutcomes.discountedValue(targetYear_)
                cPICC = testsCosts['PICCLine'].discountedValue(targetYear_)
                cUA = testsCosts['UA'].discountedValue(targetYear_)
                cUrineCulture = testsCosts['urineCulture'].discountedValue(targetYear_)
                pICU = probICU
                nContactsDayICU = nContactsICU
                nContactsDay = nContactsGenWard
                LOS_UTI = 17.5#attrLOS
                pComboWithCarb = therapyProbs['carbapenemContaining']
                pComboNoCarb = therapyProbs['carbapenemSparing']
                pMonotherapy = therapyProbs['monotherapy']
                ## Already Discounted
                cComboCarbContainingDaily = Cost(drugCostsForTherapy['carbapenemContaining'],targetYear_)
                cComboCarbSparingDaily = Cost(drugCostsForTherapy['carbapenemSparing'],targetYear_)
                cMonotherapyDaily = Cost(drugCostsForTherapy['monotherapy'],targetYear_)
                ### From Spreadsheet cells L-M 62
                cCarbContainInitialDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                cCarbSparingInitialDose = Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                cMonotherapyInitalDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                 
                tUTITreatment = treatDurDays
                
                ctCC = pComboWithCarb*((cComboCarbContainingDaily.value*tUTITreatment) + cCarbContainInitialDose.value)
                ctCS = pComboNoCarb*((cComboCarbSparingDaily.value*tUTITreatment) + cCarbSparingInitialDose.value)
                ctM  = pMonotherapy*((cMonotherapyDaily.value*tUTITreatment) + cMonotherapyInitalDose.value)
                cUA = cUA.value * testsNums['UA']['UTI'] + cUrineCulture.value * testsNums['urineCulture']['UTI']
                ggr = (cGloves_.discountedValue(targetYear_).value \
                    + cGowns_.discountedValue(targetYear_).value \
                    + (rNurseWage_.discountedValue(targetYear_).value / Constants.minutesPerHour))
            
                thirdPartyCosts = Cost(nCases*(cHospUTI.value + cPICC.value + ctCC + ctCS + ctM + cUA) \
                                        + ((pICU*nCases*nContactsDayICU \
                                        + (1.0-pICU)*nCases*nContactsDay) \
                                        * ggr*LOS_UTI),
                                       targetYear_)
                
            elif k == 'intraAbdominalInfection':
#                 nCases = 3.6529
#                 cHospAbdominal = Cost(12949,2012).discountedValue(targetYear_)#hospOutcomes.discountedValue(targetYear_)
#                 cPICC = testsCosts['PICCLine'].discountedValue(targetYear_)
#                 cWoundCulture = Cost(11.5583088,2015).discountedValue(targetYear_)#testsCosts['woundCulture'].discountedValue(targetYear_)
#                 cBloodCulture = Cost(13.59432625,2015).discountedValue(targetYear_)#testsCosts['bloodCulture'].discountedValue(targetYear_)
#                 cCTScan = Cost(253.0916333,2015).discountedValue(targetYear_)#testsCosts['CTAbdominalScan'].discountedValue(targetYear_)
#                 pICU = probICU
#                 nContactsDayICU = nContactsICU
#                 nContactsDay = nContactsGenWard
#                 LOS_IAI = 17.5 #attrLOS
#                 pComboWithCarb = therapyProbs['carbapenemContaining']
#                 pComboNoCarb = therapyProbs['carbapenemSparing']
#                 pMonotherapy = therapyProbs['monotherapy']
#                 ## Already Discounted
#                 cComboCarbContainingDaily = Cost(273.82678,targetYear_)#Cost(drugCostsForTherapy['carbapenemContaining'],targetYear_)
#                 cComboCarbSparingDaily = Cost(393.5630219,targetYear_)#Cost(drugCostsForTherapy['carbapenemSparing'],targetYear_)
#                 cMonotherapyDaily = Cost(136.9133886,targetYear_)#Cost(drugCostsForTherapy['monotherapy'],targetYear_)
#                 ### From Spreadsheet cells L-M 62
#                 cCarbContainInitialDose = Cost(0,targetYear_)#Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
#                                               #targetYear_)
#                 cCarbSparingInitialDose = Cost(303.4237654,targetYear_)#Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
#                                           #    targetYear_)
#                 cMonotherapyInitalDose = Cost(0,targetYear_)#Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
#                                          #     targetYear_)
                 
#                tIAITreatment = 17.5 #treatDurDays
                cHospAbdominal = hospOutcomes.discountedValue(targetYear_)
                cPICC = testsCosts['PICCLine'].discountedValue(targetYear_)
                cWoundCulture = testsCosts['woundCulture'].discountedValue(targetYear_)
                cBloodCulture = testsCosts['bloodCulture'].discountedValue(targetYear_)
                cCTScan = testsCosts['CTAbdominalScan'].discountedValue(targetYear_)
                pICU = probICU
                nContactsDayICU = nContactsICU
                nContactsDay = nContactsGenWard
                LOS_IAI = attrLOS
                pComboWithCarb = therapyProbs['carbapenemContaining']
                pComboNoCarb = therapyProbs['carbapenemSparing']
                pMonotherapy = therapyProbs['monotherapy']
                ## Already Discounted
                cComboCarbContainingDaily = Cost(drugCostsForTherapy['carbapenemContaining'],targetYear_)
                cComboCarbSparingDaily = Cost(drugCostsForTherapy['carbapenemSparing'],targetYear_)
                cMonotherapyDaily = Cost(drugCostsForTherapy['monotherapy'],targetYear_)
                ### From Spreadsheet cells L-M 62
                cCarbContainInitialDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                cCarbSparingInitialDose = Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                cMonotherapyInitalDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                  
                tIAITreatment = treatDurDays
                
                ctCC = pComboWithCarb*((cComboCarbContainingDaily.value*tIAITreatment) + cCarbContainInitialDose.value)
                ctCS = pComboNoCarb*((cComboCarbSparingDaily.value*tIAITreatment) + cCarbSparingInitialDose.value)
                ctM  = pMonotherapy*((cMonotherapyDaily.value*tIAITreatment) + cMonotherapyInitalDose.value)
                #cTests = cWoundCulture.value * 2.0 + cBloodCulture.value * 2.0 + cCTScan.value * 2.0
                cTests = cWoundCulture.value * testsNums['woundCulture']['intraAbdominalInfection'] \
                        + cBloodCulture.value * testsNums['bloodCulture']['intraAbdominalInfection'] \
                        + cCTScan.value * testsNums['CTAbdominalScan']['intraAbdominalInfection']
                #cUA.value * testsNums['UA']['UTI'] + cUrineCulture.value * testsNums['urineCulture']['UTI']
                ggr = (cGloves_.discountedValue(targetYear_).value \
                    + cGowns_.discountedValue(targetYear_).value \
                    + (rNurseWage_.discountedValue(targetYear_).value / Constants.minutesPerHour))
            
                thirdPartyCosts = Cost(nCases*(cHospAbdominal.value + cPICC.value + ctCC + ctCS + ctM + cTests) \
                                        + ((pICU*nCases*nContactsDayICU \
                                        + (1.0-pICU)*nCases*nContactsDay) \
                                        * ggr*LOS_IAI),
                                       targetYear_)
            
            elif k == "ventilatorAssociatedPneumonia":
                cHospVAP = hospOutcomes.discountedValue(targetYear_)
                cPICC = testsCosts['PICCLine'].discountedValue(targetYear_)
                cChestXRay = testsCosts['chestXRay'].discountedValue(targetYear_)
                cBronc = testsCosts['bronchoscopy'].discountedValue(targetYear_)
                cSputumCulture = testsCosts['sputumCultures'].discountedValue(targetYear_)
                cBloodCulture = testsCosts['bloodCulture'].discountedValue(targetYear_)
                pICU = probICU
                nContactsDayICU = nContactsICU
                nContactsDay = nContactsGenWard
                LOS_VAP = attrLOS
                pComboWithCarb = therapyProbs['carbapenemContaining']
                pComboNoCarb = therapyProbs['carbapenemSparing']
                pMonotherapy = therapyProbs['monotherapy']
                ## Already Discounted
                cComboCarbContainingDaily = Cost(drugCostsForTherapy['carbapenemContaining'],targetYear_)
                cComboCarbSparingDaily = Cost(drugCostsForTherapy['carbapenemSparing'],targetYear_)
                cMonotherapyDaily = Cost(drugCostsForTherapy['monotherapy'],targetYear_)
                ### From Spreadsheet cells L-M 62
                cCarbContainInitialDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                cCarbSparingInitialDose = Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                cMonotherapyInitalDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                  
                tVAPTreatment = treatDurDays
                
                ctCC = pComboWithCarb*((cComboCarbContainingDaily.value*tVAPTreatment) + cCarbContainInitialDose.value)
                ctCS = pComboNoCarb*((cComboCarbSparingDaily.value*tVAPTreatment) + cCarbSparingInitialDose.value)
                ctM  = pMonotherapy*((cMonotherapyDaily.value*tVAPTreatment) + cMonotherapyInitalDose.value)
                #cTests = cWoundCulture.value * 2.0 + cBloodCulture.value * 2.0 + cCTScan.value * 2.0
                cTests = cChestXRay.value * testsNums['chestXRay']['ventilationAssociatedPneumonia'] \
                        + cBronc.value * testsNums['bronchoscopy']['ventilationAssociatedPneumonia'] \
                        + cSputumCulture.value * testsNums['sputumCultures']['ventilationAssociatedPneumonia']\
                        + cBloodCulture.value * testsNums['bloodCulture']['ventilationAssociatedPneumonia']
                #cUA.value * testsNums['UA']['UTI'] + cUrineCulture.value * testsNums['urineCulture']['UTI']
                ggr = (cGloves_.discountedValue(targetYear_).value \
                    + cGowns_.discountedValue(targetYear_).value \
                    + (rNurseWage_.discountedValue(targetYear_).value / Constants.minutesPerHour))
            
                thirdPartyCosts = Cost(nCases*(cHospVAP.value + cPICC.value + ctCC + ctCS + ctM + cTests) \
                                        + ((pICU*nCases*nContactsDayICU \
                                        + (1.0-pICU)*nCases*nContactsDay) \
                                        * ggr*LOS_VAP),
                                       targetYear_)
            elif k == "nonVentilatorAssociatedPneumonia":
                cHospNVAP = hospOutcomes.discountedValue(targetYear_)
                cPICC = testsCosts['PICCLine'].discountedValue(targetYear_)
                cChestXRay = testsCosts['chestXRay'].discountedValue(targetYear_)
                cSputumCulture = testsCosts['sputumCultures'].discountedValue(targetYear_)
                cBloodCulture = testsCosts['bloodCulture'].discountedValue(targetYear_)
                pICU = probICU
                nContactsDayICU = nContactsICU
                nContactsDay = nContactsGenWard
                LOS_NVAP = attrLOS
                pComboWithCarb = therapyProbs['carbapenemContaining']
                pComboNoCarb = therapyProbs['carbapenemSparing']
                pMonotherapy = therapyProbs['monotherapy']
                ## Already Discounted
                cComboCarbContainingDaily = Cost(drugCostsForTherapy['carbapenemContaining'],targetYear_)
                cComboCarbSparingDaily = Cost(drugCostsForTherapy['carbapenemSparing'],targetYear_)
                cMonotherapyDaily = Cost(drugCostsForTherapy['monotherapy'],targetYear_)
                ### From Spreadsheet cells L-M 62
                cCarbContainInitialDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                cCarbSparingInitialDose = Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                cMonotherapyInitalDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear_).value,
                                              targetYear_)
                  
                tNVAPTreatment = treatDurDays
                
                ctCC = pComboWithCarb*((cComboCarbContainingDaily.value*tNVAPTreatment) + cCarbContainInitialDose.value)
                ctCS = pComboNoCarb*((cComboCarbSparingDaily.value*tNVAPTreatment) + cCarbSparingInitialDose.value)
                ctM  = pMonotherapy*((cMonotherapyDaily.value*tNVAPTreatment) + cMonotherapyInitalDose.value)
                cTests = cChestXRay.value * testsNums['chestXRay']['nonventilationAssociatedPneumonia'] \
                        + cSputumCulture.value * testsNums['sputumCultures']['nonventilationAssociatedPneumonia']\
                        + cBloodCulture.value * testsNums['bloodCulture']['nonventilationAssociatedPneumonia']
                #cTests = cChestXRay.value * 2.0 + cSputumCulture.value  + cBloodCulture.value
                ggr = (cGloves_.discountedValue(targetYear_).value \
                    + cGowns_.discountedValue(targetYear_).value \
                    + (rNurseWage_.discountedValue(targetYear_).value / Constants.minutesPerHour))
            
                #print cTests
                thirdPartyCosts = Cost(nCases*(cHospNVAP.value + cPICC.value + ctCC + ctCS + ctM + cTests) \
                                        + ((pICU*nCases*nContactsDayICU \
                                        + (1.0-pICU)*nCases*nContactsDay) \
                                        * ggr*LOS_NVAP),
                                       targetYear_)
            
            outcomeCosts[k]['thirdParty Costs'] = thirdPartyCosts
                        
    return outcomeCosts
    
    
def extractNewColonized(abbrev,specialDict):
    ncList = getTimeSeriesList(abbrev,specialDict,'localtiernewcolonized')
    sums = {}
    for dayVec,curves in ncList:
        for tpl, curve in curves.items():
            sums[CareTier.names[tpl]] = sum(curve)
    
    return sums

def extractContactPrecautionDays(abbrev,specialDict):
    cpList = getTimeSeriesList(abbrev,specialDict,'localtierCP')
    sums = {}
    for dayVec,curves in cpList:
        for tpl,curve in curves.items():
            sums[CareTier.names[tpl]] = sum(curve)
    return sums

def extractCREBundlesGiven(abbrev,specialDict):
    cpList = getTimeSeriesList(abbrev,specialDict,'localtierCREBundle')
    sums = {}
    for dayVec,curves in cpList:
        for tpl,curve in curves.items():
            sums[CareTier.names[tpl]] = sum(curve)
    return sums   

def extractXDROCounts(abbrev,specialDict):
    cpList = print_xdro_counts.getTimeSeriesList(abbrev, specialDict, 'localtierarrivals')
    sumT = 0.0
    for dayVec,curves in cpList:
        for tpl,curve in curves.items():
            sumT += sum(curve)
    return sumT

def getNewCols(abbrevs,note,i,facDict,xdroabbrevs,newColsReturn):
    print note
    specialDict = mergeNotesFiles([note], False)
    specialXDRODict = print_xdro_counts.mergeNotesFiles([note], False)
    newColL = {}
    newColAll = 0.0
    for abbrev in abbrevs:
        
        newColDict = extractNewColonized(abbrev,specialDict)
        for k,v in newColDict.items():
            newColAll += v
    
    cpDaysByType = {}
    for abbrev in abbrevs:
        cpForThisPlace = extractContactPrecautionDays(abbrev,specialDict)
        facCat = facDict[abbrev]['category']
        if facCat not in cpDaysByType.keys():
            cpDaysByType[facCat] = 0.0
        for k,v in cpForThisPlace.items():
            cpDaysByType[facCat] += v 
 
   
    creBundlesAll = 0.0
    
    for abbrev in abbrevs:
        creBundlesDict = extractCREBundlesGiven(abbrev, specialDict)
        for k,v in creBundlesDict.items():
            creBundlesAll += v

    xdroAll = 0.0
    for abbrev in abbrevs:
        if abbrev in xdroabbrevs:
            xdroAll += extractXDROCounts(abbrev,specialXDRODict)
    newColsReturn[i] = (newColAll,cpDaysByType,creBundlesAll,xdroAll)	 

def determin_costs(newColsArray, cpDaysArray, creBundlesArray, xdroArray, params, opts, fracAttribMort, notes, i, costsOfRealizationTmp): 
    costsOfReal = {}
    randomInt = random.randint(0,len(notes)-1)
    newColAll = newColsArray[randomInt]
    cpDaysByType = cpDaysArray[randomInt]
    creBundlesAll = creBundlesArray[randomInt]
    xdroAll = xdroArray[randomInt] 
    personAge = int(round(Distribution('uniform',{'low':60.0,'high':80.0}).draw()))
    annualWage = Cost(Distribution(params['costs']['annualWage']['distribution']['type'],
                                   params['costs']['annualWage']['distribution']['args']).draw(),
                      params['costs']['annualWage']['year'])


    hourlyWage = Cost(Distribution(params['costs']['hourlyWage']['distribution']['type'],
                                   params['costs']['hourlyWage']['distribution']['args']).draw(),
                      params['costs']['hourlyWage']['year'])
    
    
    
    ## for now, compute some common costs:
    cGlovesDict = params['interventionParameters']['glovesPair']['cost']
    cGloves = Cost(Distribution(cGlovesDict['distribution']['type'],
                                cGlovesDict['distribution']['args']).draw(),
                   cGlovesDict['year'])
            
    cGownsDict = params['interventionParameters']['gowns']['cost']
    cGowns = Cost(Distribution(cGownsDict['distribution']['type'],
                               cGownsDict['distribution']['args']).draw(),
                  cGownsDict['year'])
    
    rNurseWageDict = params['costs']['regNurseHourlyWage']
    rNurseWage = Cost(Distribution(rNurseWageDict['distribution']['type'],
                                   rNurseWageDict['distribution']['args']).draw(),
                      rNurseWageDict['year'])
    if opts.xdroyamlfile: 
        xdroCosts = computeCostOfXDROReg(xdroAll,opts.targetyear, params) # replace with RHEA
    else:
        xdroCosts = Cost(0.0,opts.targetyear)
    contPrecCosts = Cost(0.0,opts.targetyear)
    for fType,cpDays in cpDaysByType.items():
        cpst = computeCostsOfContactPrecautions(cpDays, fType, rNurseWage, cGloves, cGowns, opts.targetyear, params) 
        contPrecCosts += cpst 
    bundleCosts = computeCostsOfBundles(creBundlesAll, 0, opts.targetyear, params) # replace with RHEA
    outcomesDict = determineOutcomes(newColAll,0,infectionGivenCol,fracAttribMort,cGloves,cGowns,rNurseWage, personAge, annualWage, hourlyWage, opts.targetyear, params) # Replace with RHEA
    #costsOfRealization[i]['outcomesDict'] = outcomesDict
    
    totalSoc = 0.0
    totalQALY = 0.0
    totalHosp = 0.0
    total3rd = 0.0
    totalCases = 0.0
    totalDead = 0.0
    
    for oC,oCD in outcomesDict.items():
        if oC not in ['pneumoniaAll']:
            #print oC
            totalCases += oCD['nCases']
            totalDead += oCD['nDead']
            totalQALY += oCD['QALYs Lost'].value
            totalHosp += oCD['Hospitalization Cost'].value
            total3rd += oCD['thirdParty Costs'].value
            totalSoc += oCD['thirdParty Costs'].value + oCD['Productivity Lost'].value + oCD['Losses Due to Mortality'].value
    
    
    costsOfReal['NewColonizations'] = newColAll       
    costsOfReal['Intervention Costs'] = Cost(xdroCosts.value + contPrecCosts.value + bundleCosts.value,opts.targetyear)
    costsOfReal['QALYs Lost'] = Cost(totalQALY,opts.targetyear)
    costsOfReal['Hospitalization Costs'] = Cost(totalHosp,opts.targetyear)
    costsOfReal['thirdParty Costs'] = Cost(total3rd,opts.targetyear)
    costsOfReal['Societal Costs'] = Cost(totalSoc,opts.targetyear)
    costsOfReal['nCases'] = totalCases
    costsOfReal['nDead'] = totalDead
    costsOfReal['cpDays'] = sum([x for k,x in cpDaysByType.items()])
    costsOfReal['creBundles'] = creBundlesAll
    costsOfRealizationTmp[i] = costsOfReal
    
def main():
    
    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl [--out outname.yaml] run_descr.yaml
    """)
    
    parser.add_option('-n', '--notes', action='append', type='string',
                     help="Notes filename - may be repeated")
    parser.add_option('-o', '--out', action='store', type='string',
                     help="Output file (defaults to %s)" % DEFAULT_OUT_FILE)
    parser.add_option('-g','--glob', action='store_true',
                      help=("Apply filename globbing for notes files."
                            "  (Remember to protect the filename string from the shell!)"))
    parser.add_option('-r','--nreals',type='int',default=1000,
                     help='number of realizations of the costing model to run through')
    parser.add_option('-y','--targetyear',type='int',default=2017,
                      help='year that you would like the costing to be computed in')
    parser.add_option('-x','--xdroyamlfile',type='string', default=None,
                      help="specify the xdro scenario yaml file if there is one, if not, XDRO won't be costed")
    parser.add_option('-m','--nprocs',type='int',default=1,
                      help='number of cpus to run the costing model over')
    
    opts, args = parser.parse_args()
    
    inputDict = checkInputFileSchema(args[0], os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    if opts.out:
        outFileName = opts.out
    else:
        outFileName = DEFAULT_OUT_FILE
        
    if 'trackedFacilities' not in inputDict or not len(inputDict['trackedFacilities']):
        raise RuntimeError('Run description file does not list any tracked facilities')

    modelDir = inputDict['modelDir']
    pyrheautils.PATH_STRING_MAP['MODELDIR'] = os.path.abspath(modelDir)
    implDir = pyrheautils.pathTranslate(inputDict['facilityImplementationDir'])
    pyrheautils.PATH_STRING_MAP['IMPLDIR'] = implDir

    facDirList = [pyrheautils.pathTranslate(fPth) for fPth in inputDict['facilityDirs']]
    facDict = readFacFiles(facDirList)
    
    notes = []
    if opts.glob:
        notes = glob.glob('{0}'.format(opts.notes[0]))
    else:
        notes = opts.notes
    xdroabbrevs = []
    if opts.xdroyamlfile:
        with open(opts.xdroyamlfile,'rb') as f:
            xdroparams = yaml.load(f)
    xdroabbrevs = xdroparams['locationsImplementingScenario']['locAbbrevList']

    newColsArray = [0.0 for x in range(0,len(notes))]
    cpDaysArray = [0.0 for x in range(0,len(notes))]
    creBundlesArray = [0.0 for x in range(0,len(notes))]
    xdroArray = [0.0 for x in range(0,len(notes))] 
    abbrevs = []
    if 'trackedFacilities' in inputDict:
        for abbrev in inputDict['trackedFacilities']:
            if abbrev in facDict:
                abbrevs.append(abbrev)
                
    nprocs = opts.nprocs
    
    ### This part is the parallel reading of input data into the array
    for i in range(0,len(notes),nprocs):
        manager = Manager()
        newColsReturn = manager.dict()
        jobs = []
        end = i + nprocs
        if end > len(notes):
            end = len(notes)
        
        for j in range(i,end):
            p = Process(target=getNewCols,args = (abbrevs,notes[j],j,facDict,xdroabbrevs,newColsReturn))
            jobs.append(p)
            p.start()
        
        for proc in jobs:
            proc.join()
        for k,v in newColsReturn.items():
            newColsArray[k] = v[0]
            cpDaysArray[k] = v[1]
            creBundlesArray[k] = v[2] 
            xdroArray[k] = v[3]
        #print cpDaysArray   
        #print xdroArray
        
        
    with open("costModel_ChicagoLand.yaml","rb") as f:
        params = yaml.load(f)
    
    nReals = opts.nreals
    
    costs = {'Societal Costs': {x:{} for x in fractionAttribMorts}, 
             'Intervention Costs': {x:{} for x in fractionAttribMorts},
             'thirdParty Costs':{x:{} for x in fractionAttribMorts},
             'Hospitalization Costs':{x:{} for x in fractionAttribMorts},
             'Colonizations':{x:{} for x in fractionAttribMorts}, 
             'Deaths': {x:{} for x in fractionAttribMorts},
             'CPDays': {x:{} for x in fractionAttribMorts}, 
             'QALYs Lost': {x:{} for x in fractionAttribMorts},
             'CRE Bundles': {x:{} for x in fractionAttribMorts}
             }
    
    costsOfRealization = [{} for x in range(0,nReals)]
    for fracAttribMort in fractionAttribMorts:
        for i in range(0,nReals,nprocs):
            manager = Manager()
            costsOfRealizationTmp = manager.dict()
            jobs = []
            end = i + nprocs
            if end > nReals:
                end = nReals
                
            for j in range(i,end):
                p = Process(target=determin_costs, args=(newColsArray, cpDaysArray, creBundlesArray, xdroArray, params, opts, fracAttribMort, notes, j, costsOfRealizationTmp ))
                jobs.append(p)
                p.start()
            for proc in jobs:
                proc.join()
        
            for k,v in costsOfRealizationTmp.items():
                costsOfRealization[k] = {}
            
                for h,w in v.items():
                    costsOfRealization[k][h] = w
          
        a = np.array([costsOfRealization[i]['Societal Costs'].value for i in range(0,nReals)])
        costs['Societal Costs'][fracAttribMort]['mean'] = float(np.mean(a))
        costs['Societal Costs'][fracAttribMort]['median'] = float(np.median(a))
        costs['Societal Costs'][fracAttribMort]['stdev'] =  float(np.std(a,ddof=1))
        costs['Societal Costs'][fracAttribMort]['5%CI'] =  float(st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[0])
        costs['Societal Costs'][fracAttribMort]['95%CI'] =  float(st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[1])
        
        a = np.array([costsOfRealization[i]['Intervention Costs'].value for i in range(0,nReals)])
        costs['Intervention Costs'][fracAttribMort]['mean'] =  float( np.mean(a))
        costs['Intervention Costs'][fracAttribMort]['median'] =  float(np.median(a))
        costs['Intervention Costs'][fracAttribMort]['stdev'] =  float(np.std(a,ddof=1))
        costs['Intervention Costs'][fracAttribMort]['5%CI'] =  float(st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[0])
        costs['Intervention Costs'][fracAttribMort]['95%CI'] =  float(st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[1])
    
        a = np.array([costsOfRealization[i]['Hospitalization Costs'].value for i in range(0,nReals)])
        costs['Hospitalization Costs'][fracAttribMort]['mean'] =  float(np.mean(a))
        costs['Hospitalization Costs'][fracAttribMort]['median'] =  float(np.median(a))
        costs['Hospitalization Costs'][fracAttribMort]['stdev'] =  float(np.std(a,ddof=1))
        costs['Hospitalization Costs'][fracAttribMort]['5%CI'] =  float(st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[0])
        costs['Hospitalization Costs'][fracAttribMort]['95%CI'] =  float(st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[1])
        
        a = np.array([costsOfRealization[i]['thirdParty Costs'].value for i in range(0,nReals)])
        costs['thirdParty Costs'][fracAttribMort]['mean'] =  float(np.mean(a))
        costs['thirdParty Costs'][fracAttribMort]['median'] =  float(np.median(a))
        costs['thirdParty Costs'][fracAttribMort]['stdev'] =  float(np.std(a,ddof=1))
        costs['thirdParty Costs'][fracAttribMort]['5%CI'] =  float(st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[0])
        costs['thirdParty Costs'][fracAttribMort]['95%CI'] =  float(st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[1])
        
        a = np.array([costsOfRealization[i]['QALYs Lost'].value for i in range(0,nReals)])
        costs['QALYs Lost'][fracAttribMort]['mean'] =  float(np.mean(a))
        costs['QALYs Lost'][fracAttribMort]['median'] =  float(np.median(a))
        costs['QALYs Lost'][fracAttribMort]['stdev'] =  float(np.std(a,ddof=1))
        costs['QALYs Lost'][fracAttribMort]['5%CI'] =  float(st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[0])
        costs['QALYs Lost'][fracAttribMort]['95%CI'] =  float(st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[1])
        
        a = np.array([costsOfRealization[i]['nCases'] for i in range(0,nReals)])
        costs['Colonizations'][fracAttribMort]['mean'] =  float(np.mean(a))
        costs['Colonizations'][fracAttribMort]['median'] = np.median(a)
        costs['Colonizations'][fracAttribMort]['stdev'] = np.std(a,ddof=1)
        costs['Colonizations'][fracAttribMort]['5%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[0]
        costs['Colonizations'][fracAttribMort]['95%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[1]
        a = np.array([costsOfRealization[i]['cpDays'] for i in range(0,nReals)])
        costs['CPDays'][fracAttribMort]['mean'] =  float(np.mean(a))
        costs['CPDays'][fracAttribMort]['median'] = np.median(a)
        costs['CPDays'][fracAttribMort]['stdev'] = np.std(a,ddof=1)
        costs['CPDays'][fracAttribMort]['5%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[0]
        costs['CPDays'][fracAttribMort]['95%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[1]
        
        a = np.array([costsOfRealization[i]['creBundles'] for i in range(0,nReals)])
        costs['CRE Bundles'][fracAttribMort]['mean'] =  float(np.mean(a))
        costs['CRE Bundles'][fracAttribMort]['median'] = np.median(a)
        costs['CRE Bundles'][fracAttribMort]['stdev'] = np.std(a,ddof=1)
        costs['CRE Bundles'][fracAttribMort]['5%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[0]
        costs['CRE Bundles'][fracAttribMort]['95%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[1]
        
        a = np.array([costsOfRealization[i]['nDead'] for i in range(0,nReals)])
        costs['Deaths'][fracAttribMort]['mean'] =  float(np.mean(a))
        costs['Deaths'][fracAttribMort]['median'] = np.median(a)
        costs['Deaths'][fracAttribMort]['stdev'] = np.std(a,ddof=1)
        costs['Deaths'][fracAttribMort]['5%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[0]
        costs['Deaths'][fracAttribMort]['95%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))[1]
         
    for f,g in costs.items():
        for k,d in g.items():
            print "{6},{0},{1},{2},{3},{4},{5}".format(k,d['mean'],d['median'],d['stdev'],d['5%CI'],d['95%CI'],f)
            
    with open(outFileName,'wb') as f:
        yaml.dump(costs,f,indent=4,encoding='utf-8',width=130,explicit_start=True)
        
if __name__ == "__main__":
    main()
