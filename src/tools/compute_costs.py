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
import yaml

import numpy as np
import scipy.stats as st

targetYear = 2017
incidentCarriers = 393.41
numXDROReg = 393.41
numCREBund = 9682.83
totCarriersGettingCRE = 17.95
carrierCPdays = 11.21
nonCarrierCPDays = 3350.67
fractionAttribMort = 0.26
infectionGivenCol = 0.15

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
        

def calculateNPVValues(startAge_,annualWage_,params_,discountRate_=0.03):
    
    annualWageDiscounted = annualWage_.discountedValue(targetYear)

    
    lifeExpAtAge = int(round(params_['lifetimeProdLossParams']['expectancyAtAge'][startAge_]))
    
    #print lifeExpAtAge-startAge_
    npvProdLoss = Cost(sum([annualWageDiscounted.valueInYear(targetYear + x).value for x in range(0,lifeExpAtAge)]),targetYear)
    
    #print npvProdLoss
    
    qWeight = 1.0
    for i in params_['lifetimeProdLossParams']['npvQALYWeight']['startAges']:
        #print "{0}: {1}".format(startAge_,i)
        if i > startAge_:
            break
        qWeight = params_['lifetimeProdLossParams']['npvQALYWeight']['weights'][params_['lifetimeProdLossParams']['npvQALYWeight']['startAges'].index(i)]
    
    
    npvQ = 0.0
    for i in range(0,lifeExpAtAge):
        npvQ += Cost(qWeight,targetYear).valueInYear(targetYear+i).value
    
    npvQCost = Cost(npvQ,targetYear)
    
    return npvProdLoss, npvQCost, Cost(qWeight,targetYear)
 
def computeCostOfXDROReg(numPatients_,params_):
    
    tRegistryLogin = params_['interventionParameters']['personnelTimeXDROLoginMinutes']
    tRegistryPerPatient = params_['interventionParameters']['personnelTimeXDROSearchMinutesPerPatient']
    
    cInfectControlWageDict = params_['interventionParameters']['personnelXDROWage']
    cInfectControlWage = Cost(Distribution(cInfectControlWageDict['cost']['distribution']['type'],
                                           cInfectControlWageDict['cost']['distribution']['args']).draw(),
                              cInfectControlWageDict['cost']['year'])
    
    value = ((tRegistryLogin*Constants.daysPerYear) + (numPatients_*tRegistryPerPatient))\
          * (cInfectControlWage.discountedValue(targetYear).value/Constants.minutesPerHour)

    return Cost(value,targetYear)

def getNContactsPerDayByType(locationType_,params_):
    contactRates = params_['contacts']
    if locationType_ == 'LTACH':
        locationType_ = 'GeneralWards'
    ### Update these with RHEA monikers, may need a translate
    return Distribution(contactRates['perDayIn{0}'.format(locationType_)]['distribution']['type'],
                        contactRates['perDayIn{0}'.format(locationType_)]['distribution']['args']).draw()
                        
def computeCostsOfContactPrecautions(contactPrecautionDays_, locationType_, rNurseWage_, cGloves_, cGowns_, params_):
    
    contactRate = getNContactsPerDayByType(locationType_,params_)
    
    
    
    tContactPrecautionDict = params_['interventionParameters']['timeDonDoffMinutes']
    tCPMinutes = Distribution(tContactPrecautionDict['distribution']['type'],
                              tContactPrecautionDict['distribution']['args']).draw()
    #print "tCPMinutes = {0}".format(tCPMinutes)
    
#     cGlovesDict = params_['interventionParameters']['glovesPair']['cost']
#     cGloves = Cost(Distribution(cGlovesDict['distribution']['type'],
#                                 cGlovesDict['distribution']['args']).draw(),
#                    cGlovesDict['year'])
#     #print 'cGloves = {0}'.format(cGloves)
#     
#     cGownsDict = params_['interventionParameters']['gowns']['cost']
#     cGowns = Cost(Distribution(cGownsDict['distribution']['type'],
#                                cGownsDict['distribution']['args']).draw(),
#                    cGownsDict['year'])
    
    
    #print (contactPrecautionDays_ * contactRate)
    value = (contactPrecautionDays_ \
          * contactRate) \
          * ((rNurseWage_.discountedValue(targetYear).value/60.0) \
          * tCPMinutes \
          + cGloves_.discountedValue(targetYear).value \
          + cGowns_.discountedValue(targetYear).value)
    
    return Cost(value,targetYear)

def computeCostsOfBundles(numBundles_,numBaths_,params_):
    chgWipeDict = params_['interventionParameters']['chgWipesPerBath']['cost']
    cChgWipesPerBath = Cost(Distribution(chgWipeDict['distribution']['type'],
                                        chgWipeDict['distribution']['args']).draw(),
                           chgWipeDict['year'])
    
    screenDict = params_['interventionParameters']['screeningTotal']['cost']
    cScreening = Cost(Distribution(screenDict['distribution']['type'],
                                   screenDict['distribution']['args']).draw(),
                           screenDict['year'])
    
    #print cChgWipesPerBath.discountedValue(targetYear)
    #print cScreening.discountedValue(targetYear)
    value = numBundles_ * cScreening.discountedValue(targetYear).value \
          + numBaths_ * cChgWipesPerBath.discountedValue(targetYear).value
    
    return Cost(value, targetYear)      
        
                                        
def determineOutcomes(nIncidence_,nCarriersCRE_, probInfect_, attribMort_, cGloves_, cGowns_, rNurseWage_, personAge_, annualWage_, hourlyWage_, params_):
    outcomesDict = params_['outcomes']
    outcomeCosts = {x:{} for x in outcomesDict.keys()}
    
    #nOutsDict = {x:Distribution(outcomesDict[x]['probability']['distribution']['type'],
    #                            outcomesDict[x]['probability']['distribution']['args']).draw() \
    #             * (nIncidence_+nCarriersCRE_) \
    #             * probInfect_ \
    #             for x in outcomesDict.keys()}
    therapyDict = params_['therapy']
    
    npv,npQ,baseQALY = calculateNPVValues(personAge_, annualWage_, params_)
    
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
                
            
                         
                nDead += Distribution(mD['distribution']['type'],
                                      mD['distribution']['args']).draw() \
                       * nCases * attribMort_ * therapyProbs[mk]
            
                
                drugCostsForTherapy = {}
                for th,thD in therapyDrugReg.items():
                    drugCostsForTherapy[th] = 0.0
                    for d,dD in thD.items():
                        drugCost = drugCosts[d]
                        drugCostsForTherapy[th] += drugCost.discountedValue(targetYear).value * dD
                
            outcomeCosts[k]['nDead'] = nDead    
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
            ggr = cGloves_.discountedValue(targetYear).value \
                + cGowns_.discountedValue(targetYear).value \
                + (rNurseWage_.discountedValue(targetYear).value / Constants.minutesPerHour)
#             print ggr  
            cggrLOSICU = (cICUBedDay.discountedValue(targetYear).value + (nContactsICU*ggr))*attrLOS
            cggrLOSGW =  (cGenWardBedDay.discountedValue(targetYear).value + (nContactsGenWard*ggr))*attrLOS
#             print "CICU = {0}".format(nContactsICU*ggr)
#             print "c2ICU = {0}".format((cICUBedDay.discountedValue(targetYear).value + (nContactsICU*ggr)))
#             print "c3ICU = {0}".format(nICU*cggrLOSICU)
#             print "cG = {0}".format((nContactsGenWard*ggr))
#             print "cG22={0}".format((cGenWardBedDay.discountedValue(targetYear).value + (nContactsGenWard*ggr)))
#             print "c2G = {0}".format((cGenWardBedDay.discountedValue(targetYear).value + (nContactsGenWard*ggr))*attrLOS)
#             print "c3G = {0}".format(nGenW*cggrLOSGW)
#             print cggrLOSICU
            hospCost= Cost(nICU*cggrLOSICU + nGenW*cggrLOSGW,targetYear)       
            
            outcomeCosts[k]['Hospitalization Cost'] = hospCost
            #print "{0} Hosp = {1}".format(k,hospCost)
            
            utilityWeight = Distribution(v['utilityWeight']['distribution']['type'],
                                         v['utilityWeight']['distribution']['args']).draw()
            QALYLost = Cost(nCases*((utilityWeight*baseQALY.discountedValue(targetYear).value)/Constants.daysPerYear)*attrLOS \
                            + nDead * npQ.discountedValue(targetYear).value,
                            targetYear)
            
            outcomeCosts[k]['QALYs Lost'] = QALYLost                
            #print "QLost = {0}".format(QALYLost)
            
            ### prod Losses
            prodLosses = Cost(nCases * attrLOS * hourlyWage_.discountedValue(targetYear).value * Constants.workDayHours,
                              targetYear)
            
            outcomeCosts[k]['Productivity Lost'] = prodLosses
            #print "ProdLost {0}: {1}".format(k,prodLosses)
            
            #### Losses Mortality
            
            lossesMortality = Cost(nDead * npv.value,
                                   targetYear)
            
            outcomeCosts[k]['Losses Due to Mortality'] = lossesMortality
            
            #print "MortLoss {0}: {1}".format(k,lossesMortality)
            
            ### LETS Build the first term of 3rd Party
            thirdPartyCosts = 0.0
            if k == "bacteremia":
                ### doing this to match the complex Excel formula, can do something more sane later
                cHospBacteremia = hospOutcomes.discountedValue(targetYear)
                cPICC = testsCosts['PICCLine'].discountedValue(targetYear)
                cBloodCulture = testsCosts['bloodCulture'].discountedValue(targetYear)
                pICU = probICU
                nContactsDayICU = nContactsICU
                nContactsDay = nContactsGenWard
                LOS_Bacteremia = attrLOS
                pComboWithCarb = therapyProbs['carbapenemContaining']
                pComboNoCarb = therapyProbs['carbapenemSparing']
                pMonotherapy = therapyProbs['monotherapy']
                ## Already Discounted
                cComboCarbContainingDaily = Cost(drugCostsForTherapy['carbapenemContaining'],targetYear)
                cComboCarbSparingDaily = Cost(drugCostsForTherapy['carbapenemSparing'],targetYear)
                cMonotherapyDaily = Cost(drugCostsForTherapy['monotherapy'],targetYear)
                ### From Spreadsheet cells L-M 62
                cCarbContainInitialDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                cCarbSparingInitialDose = Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                cMonotherapyInitalDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                
                tBacteremiaTreatment = treatDurDays
                
                ctCC = pComboWithCarb*((cComboCarbContainingDaily.value*tBacteremiaTreatment) + cCarbContainInitialDose.value)
                ctCS = pComboNoCarb*((cComboCarbSparingDaily.value*tBacteremiaTreatment) + cCarbSparingInitialDose.value)
                ctM  = pMonotherapy*((cMonotherapyDaily.value*tBacteremiaTreatment) + cMonotherapyInitalDose.value)
                cBC = cBloodCulture.value * testsNums['bloodCulture']['bacteremia']
                ggr = (cGloves_.discountedValue(targetYear).value \
                    + cGowns_.discountedValue(targetYear).value \
                    + (rNurseWage_.discountedValue(targetYear).value / Constants.minutesPerHour))
            
                thirdPartyCosts = Cost(nCases*(cHospBacteremia.value + cPICC.value + ctCC + ctCS + ctM + cBC) \
                                        + ((pICU*nCases*nContactsDayICU \
                                        + (1.0-pICU)*nCases*nContactsDay) \
                                        * ggr*LOS_Bacteremia),
                                       targetYear)
                
            elif k == 'UTI':
                #nCases = 29.673531
#                 cHospUTI = Cost(6900.72812,2012).discountedValue(targetYear)#hospOutcomes.discountedValue(targetYear)
#                 cPICC = testsCosts['PICCLine'].discountedValue(targetYear)
#                 cUA = Cost(2.933778,2015).discountedValue(targetYear) #testsCosts['UA'].discountedValue(targetYear)
#                 cUrineCulture = Cost(10.6956948,2015).discountedValue(targetYear)
#                 pICU = probICU
#                 nContactsDayICU = nContactsICU
#                 nContactsDay = nContactsGenWard
#                 LOS_UTI = 6.0 #attrLOS
#                 pComboWithCarb = therapyProbs['carbapenemContaining']
#                 pComboNoCarb = therapyProbs['carbapenemSparing']
#                 pMonotherapy = therapyProbs['monotherapy']
#                 ## Already Discounted
#                 cComboCarbContainingDaily = Cost(273.82678,targetYear)#Cost(drugCostsForTherapy['carbapenemContaining'],targetYear)
#                 cComboCarbSparingDaily = Cost(393.5630219,targetYear)#Cost(drugCostsForTherapy['carbapenemSparing'],targetYear)
#                 cMonotherapyDaily = Cost(136.9133886,targetYear)#Cost(drugCostsForTherapy['monotherapy'],targetYear)
#                 ### From Spreadsheet cells L-M 62
#                 cCarbContainInitialDose = Cost(0,targetYear)#Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
#                                               #targetYear)
#                 cCarbSparingInitialDose = Cost(303.4237654,targetYear)#Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
#                                           #    targetYear)
#                 cMonotherapyInitalDose = Cost(0,targetYear)#Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
#                                          #     targetYear)
#                 
#                 tUTITreatment = 12.0 #treatDurDays
                cHospUTI = hospOutcomes.discountedValue(targetYear)
                cPICC = testsCosts['PICCLine'].discountedValue(targetYear)
                cUA = testsCosts['UA'].discountedValue(targetYear)
                cUrineCulture = testsCosts['urineCulture'].discountedValue(targetYear)
                pICU = probICU
                nContactsDayICU = nContactsICU
                nContactsDay = nContactsGenWard
                LOS_UTI = 17.5#attrLOS
                pComboWithCarb = therapyProbs['carbapenemContaining']
                pComboNoCarb = therapyProbs['carbapenemSparing']
                pMonotherapy = therapyProbs['monotherapy']
                ## Already Discounted
                cComboCarbContainingDaily = Cost(drugCostsForTherapy['carbapenemContaining'],targetYear)
                cComboCarbSparingDaily = Cost(drugCostsForTherapy['carbapenemSparing'],targetYear)
                cMonotherapyDaily = Cost(drugCostsForTherapy['monotherapy'],targetYear)
                ### From Spreadsheet cells L-M 62
                cCarbContainInitialDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                cCarbSparingInitialDose = Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                cMonotherapyInitalDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                 
                tUTITreatment = treatDurDays
                
                ctCC = pComboWithCarb*((cComboCarbContainingDaily.value*tUTITreatment) + cCarbContainInitialDose.value)
                ctCS = pComboNoCarb*((cComboCarbSparingDaily.value*tUTITreatment) + cCarbSparingInitialDose.value)
                ctM  = pMonotherapy*((cMonotherapyDaily.value*tUTITreatment) + cMonotherapyInitalDose.value)
                cUA = cUA.value * testsNums['UA']['UTI'] + cUrineCulture.value * testsNums['urineCulture']['UTI']
                ggr = (cGloves_.discountedValue(targetYear).value \
                    + cGowns_.discountedValue(targetYear).value \
                    + (rNurseWage_.discountedValue(targetYear).value / Constants.minutesPerHour))
            
                thirdPartyCosts = Cost(nCases*(cHospUTI.value + cPICC.value + ctCC + ctCS + ctM + cUA) \
                                        + ((pICU*nCases*nContactsDayICU \
                                        + (1.0-pICU)*nCases*nContactsDay) \
                                        * ggr*LOS_UTI),
                                       targetYear)
                
            elif k == 'intraAbdominalInfection':
#                 nCases = 3.6529
#                 cHospAbdominal = Cost(12949,2012).discountedValue(targetYear)#hospOutcomes.discountedValue(targetYear)
#                 cPICC = testsCosts['PICCLine'].discountedValue(targetYear)
#                 cWoundCulture = Cost(11.5583088,2015).discountedValue(targetYear)#testsCosts['woundCulture'].discountedValue(targetYear)
#                 cBloodCulture = Cost(13.59432625,2015).discountedValue(targetYear)#testsCosts['bloodCulture'].discountedValue(targetYear)
#                 cCTScan = Cost(253.0916333,2015).discountedValue(targetYear)#testsCosts['CTAbdominalScan'].discountedValue(targetYear)
#                 pICU = probICU
#                 nContactsDayICU = nContactsICU
#                 nContactsDay = nContactsGenWard
#                 LOS_IAI = 17.5 #attrLOS
#                 pComboWithCarb = therapyProbs['carbapenemContaining']
#                 pComboNoCarb = therapyProbs['carbapenemSparing']
#                 pMonotherapy = therapyProbs['monotherapy']
#                 ## Already Discounted
#                 cComboCarbContainingDaily = Cost(273.82678,targetYear)#Cost(drugCostsForTherapy['carbapenemContaining'],targetYear)
#                 cComboCarbSparingDaily = Cost(393.5630219,targetYear)#Cost(drugCostsForTherapy['carbapenemSparing'],targetYear)
#                 cMonotherapyDaily = Cost(136.9133886,targetYear)#Cost(drugCostsForTherapy['monotherapy'],targetYear)
#                 ### From Spreadsheet cells L-M 62
#                 cCarbContainInitialDose = Cost(0,targetYear)#Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
#                                               #targetYear)
#                 cCarbSparingInitialDose = Cost(303.4237654,targetYear)#Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
#                                           #    targetYear)
#                 cMonotherapyInitalDose = Cost(0,targetYear)#Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
#                                          #     targetYear)
                 
#                tIAITreatment = 17.5 #treatDurDays
                cHospAbdominal = hospOutcomes.discountedValue(targetYear)
                cPICC = testsCosts['PICCLine'].discountedValue(targetYear)
                cWoundCulture = testsCosts['woundCulture'].discountedValue(targetYear)
                cBloodCulture = testsCosts['bloodCulture'].discountedValue(targetYear)
                cCTScan = testsCosts['CTAbdominalScan'].discountedValue(targetYear)
                pICU = probICU
                nContactsDayICU = nContactsICU
                nContactsDay = nContactsGenWard
                LOS_IAI = attrLOS
                pComboWithCarb = therapyProbs['carbapenemContaining']
                pComboNoCarb = therapyProbs['carbapenemSparing']
                pMonotherapy = therapyProbs['monotherapy']
                ## Already Discounted
                cComboCarbContainingDaily = Cost(drugCostsForTherapy['carbapenemContaining'],targetYear)
                cComboCarbSparingDaily = Cost(drugCostsForTherapy['carbapenemSparing'],targetYear)
                cMonotherapyDaily = Cost(drugCostsForTherapy['monotherapy'],targetYear)
                ### From Spreadsheet cells L-M 62
                cCarbContainInitialDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                cCarbSparingInitialDose = Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                cMonotherapyInitalDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                  
                tIAITreatment = treatDurDays
                
                ctCC = pComboWithCarb*((cComboCarbContainingDaily.value*tIAITreatment) + cCarbContainInitialDose.value)
                ctCS = pComboNoCarb*((cComboCarbSparingDaily.value*tIAITreatment) + cCarbSparingInitialDose.value)
                ctM  = pMonotherapy*((cMonotherapyDaily.value*tIAITreatment) + cMonotherapyInitalDose.value)
                #cTests = cWoundCulture.value * 2.0 + cBloodCulture.value * 2.0 + cCTScan.value * 2.0
                cTests = cWoundCulture.value * testsNums['woundCulture']['intraAbdominalInfection'] \
                        + cBloodCulture.value * testsNums['bloodCulture']['intraAbdominalInfection'] \
                        + cCTScan.value * testsNums['CTAbdominalScan']['intraAbdominalInfection']
                #cUA.value * testsNums['UA']['UTI'] + cUrineCulture.value * testsNums['urineCulture']['UTI']
                ggr = (cGloves_.discountedValue(targetYear).value \
                    + cGowns_.discountedValue(targetYear).value \
                    + (rNurseWage_.discountedValue(targetYear).value / Constants.minutesPerHour))
            
                thirdPartyCosts = Cost(nCases*(cHospAbdominal.value + cPICC.value + ctCC + ctCS + ctM + cTests) \
                                        + ((pICU*nCases*nContactsDayICU \
                                        + (1.0-pICU)*nCases*nContactsDay) \
                                        * ggr*LOS_IAI),
                                       targetYear)
            
            elif k == "ventilatorAssociatedPneumonia":
                cHospVAP = hospOutcomes.discountedValue(targetYear)
                cPICC = testsCosts['PICCLine'].discountedValue(targetYear)
                cChestXRay = testsCosts['chestXRay'].discountedValue(targetYear)
                cBronc = testsCosts['bronchoscopy'].discountedValue(targetYear)
                cSputumCulture = testsCosts['sputumCultures'].discountedValue(targetYear)
                cBloodCulture = testsCosts['bloodCulture'].discountedValue(targetYear)
                pICU = probICU
                nContactsDayICU = nContactsICU
                nContactsDay = nContactsGenWard
                LOS_VAP = attrLOS
                pComboWithCarb = therapyProbs['carbapenemContaining']
                pComboNoCarb = therapyProbs['carbapenemSparing']
                pMonotherapy = therapyProbs['monotherapy']
                ## Already Discounted
                cComboCarbContainingDaily = Cost(drugCostsForTherapy['carbapenemContaining'],targetYear)
                cComboCarbSparingDaily = Cost(drugCostsForTherapy['carbapenemSparing'],targetYear)
                cMonotherapyDaily = Cost(drugCostsForTherapy['monotherapy'],targetYear)
                ### From Spreadsheet cells L-M 62
                cCarbContainInitialDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                cCarbSparingInitialDose = Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                cMonotherapyInitalDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                  
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
                ggr = (cGloves_.discountedValue(targetYear).value \
                    + cGowns_.discountedValue(targetYear).value \
                    + (rNurseWage_.discountedValue(targetYear).value / Constants.minutesPerHour))
            
                thirdPartyCosts = Cost(nCases*(cHospVAP.value + cPICC.value + ctCC + ctCS + ctM + cTests) \
                                        + ((pICU*nCases*nContactsDayICU \
                                        + (1.0-pICU)*nCases*nContactsDay) \
                                        * ggr*LOS_VAP),
                                       targetYear)
            elif k == "nonVentilatorAssociatedPneumonia":
                cHospNVAP = hospOutcomes.discountedValue(targetYear)
                cPICC = testsCosts['PICCLine'].discountedValue(targetYear)
                cChestXRay = testsCosts['chestXRay'].discountedValue(targetYear)
                cSputumCulture = testsCosts['sputumCultures'].discountedValue(targetYear)
                cBloodCulture = testsCosts['bloodCulture'].discountedValue(targetYear)
                pICU = probICU
                nContactsDayICU = nContactsICU
                nContactsDay = nContactsGenWard
                LOS_NVAP = attrLOS
                pComboWithCarb = therapyProbs['carbapenemContaining']
                pComboNoCarb = therapyProbs['carbapenemSparing']
                pMonotherapy = therapyProbs['monotherapy']
                ## Already Discounted
                cComboCarbContainingDaily = Cost(drugCostsForTherapy['carbapenemContaining'],targetYear)
                cComboCarbSparingDaily = Cost(drugCostsForTherapy['carbapenemSparing'],targetYear)
                cMonotherapyDaily = Cost(drugCostsForTherapy['monotherapy'],targetYear)
                ### From Spreadsheet cells L-M 62
                cCarbContainInitialDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                cCarbSparingInitialDose = Cost(0.60*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                cMonotherapyInitalDose = Cost(0.25*drugCosts['tigercyclineCAP50mg12hrs'].discountedValue(targetYear).value,
                                              targetYear)
                  
                tNVAPTreatment = treatDurDays
                
                ctCC = pComboWithCarb*((cComboCarbContainingDaily.value*tNVAPTreatment) + cCarbContainInitialDose.value)
                ctCS = pComboNoCarb*((cComboCarbSparingDaily.value*tNVAPTreatment) + cCarbSparingInitialDose.value)
                ctM  = pMonotherapy*((cMonotherapyDaily.value*tNVAPTreatment) + cMonotherapyInitalDose.value)
                cTests = cChestXRay.value * testsNums['chestXRay']['nonventilationAssociatedPneumonia'] \
                        + cSputumCulture.value * testsNums['sputumCultures']['nonventilationAssociatedPneumonia']\
                        + cBloodCulture.value * testsNums['bloodCulture']['nonventilationAssociatedPneumonia']
                #cTests = cChestXRay.value * 2.0 + cSputumCulture.value  + cBloodCulture.value
                ggr = (cGloves_.discountedValue(targetYear).value \
                    + cGowns_.discountedValue(targetYear).value \
                    + (rNurseWage_.discountedValue(targetYear).value / Constants.minutesPerHour))
            
                #print cTests
                thirdPartyCosts = Cost(nCases*(cHospNVAP.value + cPICC.value + ctCC + ctCS + ctM + cTests) \
                                        + ((pICU*nCases*nContactsDayICU \
                                        + (1.0-pICU)*nCases*nContactsDay) \
                                        * ggr*LOS_NVAP),
                                       targetYear)
            
            outcomeCosts[k]['thirdParty Costs'] = thirdPartyCosts
            
            #print "{0} 3rd Party: {1}".format(k,thirdPartyCosts)
                        
            
    return outcomeCosts
    
        
#     ### two special cases 
#     nOutsDict['ventilatorAssociatedPneumonia'] = Distribution(outcomesDict['ventilatorAssociatedPneumonia']['probability']['distribution']['type'],
#                                                               outcomesDict['ventilatorAssociatedPneumonia']['probability']['distribution']['args']).draw() \
#                                                * nOutsDict['pneumoniaAll']
#     nOutsDict['nonVentilatorAssociatedPneumonia'] = nOutsDict['pneumoniaAll'] - nOutsDict['ventilatorAssociatedPneumonia']
#     
#     
#     
#     print nOutsDict
                            
def main():
    with open("costModel_ChicagoLand.yaml","rb") as f:
        params = yaml.load(f)
    
    nReals = 1000
    
    costsOfRealization = [{} for x in range(nReals)]
   
    for i in range(0,nReals):
        print "starting realization {0}".format(i)
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
        
        costsOfRealization[i]['xdroCosts'] = computeCostOfXDROReg(138,params) # replace with RHEA
        costsOfRealization[i]['contPrecCosts'] = computeCostsOfContactPrecautions(3362,'LTACH',cGloves,cGowns,rNurseWage, params) # Replace with place specific members
        costsOfRealization[i]['bundleCosts'] = computeCostsOfBundles(9682.8270318252, 0, params) # replace with RHEA
        outcomesDict = determineOutcomes(393.41,17.95,infectionGivenCol,fractionAttribMort,cGloves,cGowns,rNurseWage, personAge, annualWage, hourlyWage, params) # Replace with RHEA
        costsOfRealization[i]['outcomesDict'] = outcomesDict
        
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
                
                
        costsOfRealization[i]['Intervention Costs'] = Cost(costsOfRealization[i]['xdroCosts'].value + costsOfRealization[i]['contPrecCosts'].value + costsOfRealization[i]['bundleCosts'].value,targetYear)
        costsOfRealization[i]['QALYs Lost'] = Cost(totalQALY,targetYear)
        costsOfRealization[i]['Hospitalization Costs'] = Cost(totalHosp,targetYear)
        costsOfRealization[i]['thirdParty Costs'] = Cost(total3rd,targetYear)
        costsOfRealization[i]['Societal Costs'] = Cost(totalSoc,targetYear)
        costsOfRealization[i]['nCases'] = totalCases
        costsOfRealization[i]['nDead'] = totalDead
        
        
    
    costs = {'Societal Costs': {}, 'Intervention Costs': {},'thirdParty Costs':{},'Hospitalization Costs':{}, 'Infections':{}, 'Deaths': {}, 'QALYs Lost': {}}
    
    a = np.array([costsOfRealization[i]['Societal Costs'].value for i in range(0,nReals)])
    costs['Societal Costs']['mean'] = np.mean(a)
    costs['Societal Costs']['median'] = np.median(a)
    costs['Societal Costs']['stdev'] = np.std(a,ddof=1)
    costs['Societal Costs']['95%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))
    
    a = np.array([costsOfRealization[i]['Intervention Costs'].value for i in range(0,nReals)])
    costs['Intervention Costs']['mean'] = np.mean(a)
    costs['Intervention Costs']['median'] = np.median(a)
    costs['Intervention Costs']['stdev'] = np.std(a,ddof=1)
    costs['Intervention Costs']['95%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))

    a = np.array([costsOfRealization[i]['Hospitalization Costs'].value for i in range(0,nReals)])
    costs['Hospitalization Costs']['mean'] = np.mean(a)
    costs['Hospitalization Costs']['median'] = np.median(a)
    costs['Hospitalization Costs']['stdev'] = np.std(a,ddof=1)
    costs['Hospitalization Costs']['95%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))
    
    a = np.array([costsOfRealization[i]['thirdParty Costs'].value for i in range(0,nReals)])
    costs['thirdParty Costs']['mean'] = np.mean(a)
    costs['thirdParty Costs']['median'] = np.median(a)
    costs['thirdParty Costs']['stdev'] = np.std(a,ddof=1)
    costs['thirdParty Costs']['95%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))
    
    a = np.array([costsOfRealization[i]['QALYs Lost'].value for i in range(0,nReals)])
    costs['QALYs Lost']['mean'] = np.mean(a)
    costs['QALYs Lost']['median'] = np.median(a)
    costs['QALYs Lost']['stdev'] = np.std(a,ddof=1)
    costs['QALYs Lost']['95%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))
    
    a = np.array([costsOfRealization[i]['nCases'] for i in range(0,nReals)])
    costs['Infections']['mean'] = np.mean(a)
    costs['Infections']['median'] = np.median(a)
    costs['Infections']['stdev'] = np.std(a,ddof=1)
    costs['Infections']['95%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))
    
    a = np.array([costsOfRealization[i]['nDead'] for i in range(0,nReals)])
    costs['Deaths']['mean'] = np.mean(a)
    costs['Deaths']['median'] = np.median(a)
    costs['Deaths']['stdev'] = np.std(a,ddof=1)
    costs['Deaths']['95%CI'] = st.t.interval(0.95,len(a)-1, loc=np.mean(a),scale=st.sem(a))
    
    
    for k,d in costs.items():
        print "{0},{1},{2},{3},{4},{5}".format(k,d['mean'],d['median'],d['stdev'],d['95%CI'][0],d['95%CI'][1])
        
if __name__ == "__main__":
    main()
