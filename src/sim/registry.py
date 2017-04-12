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

import logging
from random import choice
from collections import defaultdict

from phacsl.utils.collections.phacollections import SingletonMetaClass
import quilt.peopleplaces as peopleplaces
from pathogenbase import PthStatus

logger = logging.getLogger(__name__)

class RegistryUpdateMsg(peopleplaces.SimpleMsg):
    pass

class RegistryGroupUpdateMsg(peopleplaces.SimpleMsg):
    pass

class RegistryUpdateQueue(peopleplaces.RequestQueue):
    pass

class RegistryManager(peopleplaces.Manager):
    """Manager agent for the per-patch registry"""

    def __init__(self, name, patch, managementBase):
        super(RegistryManager, self).__init__(name, patch, managementBase)
        self.newRegList = []  # Format is [(patientId, condition, status), ...]

    def perTickActions(self, timeNow):
        """
        Update the registries on other patches about new registrations.
        """
        if self.newRegList:
            facAddrL = [tpl[1] for tpl in self.patch.serviceLookup('RegistryUpdateQueue')
                        if not self.patch.isLocal(tpl[1])]
            payload = self.newRegList[:]
            for fA in facAddrL:
                msg = RegistryGroupUpdateMsg(self.name + '_groupUpdateMsg',
                                             self.patch, payload, fA, debug=True)
                self.patch.launch(msg, timeNow)
            self.newRegList = []

class RegistryCore(object):
    """
    There is one instance of this class per quilt.patches.PatchGroup.
    Any agent may read its state during that agent's time slice, but
    the state can only be set by the manager during the manager's
    time slice.
    """
    __metaclass__ = SingletonMetaClass

    def __init__(self):
        self.setDict = defaultdict(set)

    def getPatientStatus(self, condition, patientId):
        """Return patient status from local cache"""
        return patientId in self.setDict[condition]

    def setLclPatientStatus(self, condition, patientId, status=True):
        """Set patient status, but only in local cache"""
        logger.debug('setLclPatientStatus %s %s %s', condition, patientId, status)
        if status:
            self.setDict[condition].add(patientId)
        elif patientId in self.setDict[condition]:
            self.setDict[condition].remove(patientId)

class Registry(peopleplaces.ManagementBase):
    """
    Implements a registry for patient status.  There is one registry per patch,
    but registry info must be shared across all patches.  It is the job of this
    class and its manager agent to keep the global representation consistent.
    """

    def __init__(self, name, patch, managerClass=None, reqQueueClasses=None):
        reqQueueClasses = (reqQueueClasses if reqQueueClasses else [RegistryUpdateQueue])
        managerClass = (managerClass if managerClass else RegistryManager)
        super(Registry, self).__init__(name, patch,
                                       managerClass=managerClass,
                                       reqQueueClasses=reqQueueClasses)
        self.core = RegistryCore()
        self.noteHolder = None

    def setNoteHolder(self, noteHolder):
        self.noteHolder = noteHolder
        self.noteHolder.addNote({'name': self.name})

    def getNoteHolder(self):
        return self.noteHolder

    def getMsgPayload(self, msgType, patientAgent):
        if issubclass(msgType, RegistryUpdateMsg):
            return (patientAgent.id,
                    str(patientAgent.ward.iA),
                    patientAgent.getPthDiagnosis() != PthStatus.CLEAR)
        else:
            return super(Registry, self).getMsgPayload(msgType, patientAgent)

    @classmethod
    def buildMsgPayload(cls, msgType, patientAgent, condition):
        """So that registry clients can build a payload without a Registry instance"""
        if issubclass(msgType, RegistryUpdateMsg):
            return (patientAgent.id, condition,
                    patientAgent.getPthDiagnosis() != PthStatus.CLEAR)
        else:
            raise RuntimeError('Cannot build payload for a %s' % msgType.__name__)

    def handleIncomingMsg(self, msgType, payload, timeNow):
        if issubclass(msgType, RegistryUpdateMsg):
            patientId, condition, status = payload
            self.core.setLclPatientStatus(condition, patientId, status)
            self.manager.newRegList.append(((patientId, condition, status)))
        elif issubclass(msgType, RegistryGroupUpdateMsg):
            # Payload format is that of registry.newRegList
            for patientId, condition, status in payload:
                self.core.setLclPatientStatus(condition, patientId, status)
        return timeNow

    @classmethod
    def getCore(cls):
        """
        Return the singleton core for this Registry
        """
        if RegistryCore in SingletonMetaClass._instances:
            return SingletonMetaClass._instances[RegistryCore]
        else:
            dummyCore = RegistryCore() # Which will persist because it is a singleton
            return dummyCore

    @classmethod
    def getPatientStatus(cls, condition, patientId):
        return cls.getCore().getPatientStatus(condition, patientId)

    @classmethod
    def registerPatientStatus(cls, patientAgent, condition, timeNow):
        facAddr = choice([tpl[1] for tpl in patientAgent.patch.serviceLookup('RegistryUpdateQueue')
                          if patientAgent.patch.group.isLocal(tpl[1])])
        payload = cls.buildMsgPayload(RegistryUpdateMsg, patientAgent, condition)
        patientAgent.patch.launch(RegistryUpdateMsg(patientAgent.name + '_registryUpdateMsg',
                                                    patientAgent.patch, payload, facAddr),
                                  timeNow)
