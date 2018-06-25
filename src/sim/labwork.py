'''
Created on Jun 25, 2018

@author: welling
'''

from random import random
import logging

from quilt.peopleplaces import FutureMsg

class LabWorkMsg(FutureMsg):
    idCounter = 0

    @classmethod
    def nextId(cls):
        rslt = cls.idCounter
        cls.idCounter += 1
        return rslt

    def __init__(self, baseName, patch, payload, destAddr, arrivalTime, debug=False):
        super(LabWorkMsg, self).__init__(baseName + '_lab_work_%d' % self.nextId(),
                                         patch, payload, destAddr, arrivalTime,
                                         debug=debug)

class LabWork(object):
    '''
    This class provides support for blood tests and such.

    The name doesn't include either blood or test to avoid horror and/or confusion.
    '''

    def __init__(self, sensitivity, specificity, delayDays, debug=False):
        '''
        sensitivity is the sensitivity (effectiveness) of the test.
        specificity is the specificity (1.0 - false positive rate) of the test
        delayDays is the number of days in the future for results arrival
        '''
        self.sensitivity = sensitivity
        self.falsePosRate = 1.0 - specificity
        self.delayDays = delayDays
        self.debug = debug

    def trueTestFun(self, patientStatus):
        """
        Returns a boolean giving the correct test result
        """
        raise RuntimeError('Derived class must implement this method')

    @classmethod
    def posAction(cls, patientRecord):
        """
        Implements what happens if the test returns True (perhaps erroneously).
        Returns the updated patientRecord.
        """
        raise RuntimeError('Derived class must implement this method')

    @classmethod
    def negAction(cls, patientRecord):
        """
        Implements what happens if the test returns False (perhaps erroneously).
        Returns the updated patientRecord.
        """
        raise RuntimeError('Derived class must implement this method')

    def performLab(self, patientId, patientStatus, ward, timeNow):
        """
        Perform the test.  The results will be added to the patient's record after a delay.

        patientId is a patient's unique ID.
        patientStatus is a PatientStatus instance.
        ward is the ward in which the test is to be performed.
        timeNow is the current time
        """
        if timeNow is None:
            timeNow = 0  # Deal with possible messages during initialization
        trueRslt = self.trueTestFun(patientStatus)
        if trueRslt:
            rslt = (random() <= self.sensitivity)
        else:
            rslt = (random() <= self.falsePosRate)
        labWorkMsg = LabWorkMsg(ward._name, ward.patch, (patientId, rslt), ward.getReqQueueAddr(),
                                timeNow + self.delayDays, self.debug)
        ward.patch.launch(labWorkMsg, timeNow)

    @classmethod
    def handleLabMsg(cls, fac, msgType, payload, timeNow):
        print 'handler for %s' % cls.__name__
        assert cls == msgType, 'lab message of type %s passed to %s' % (msgType.__name__,
                                                                        cls.__name__)
        patientId, rslt = payload
        with fac.getPatientRecord(patientId) as pRec:
            if rslt:
                pRec = cls.posAction(pRec)
            else:
                pRec = cls.negAction(pRec)
                
                
class SwabTest(LabWork):
    def __init__(self, debug=False):
        super(self, SwabTest).__init__(_constants['swabDiagnosticSensitivity']['value'],
                                       _constants['swabDiagnosticSpecificity']['value'],
                                       _constants['swabDelayDays']['value'],
                                       debug=debug)
    
    @classmethod
    def posAction(cls, patientRecord):
        patientRecord.carriesPth = True
        return patientRecord

    @classmethod
    def negAction(cls, patientRecord):
        patientRecord.carriesPth = False
        return patientRecord

                
                
