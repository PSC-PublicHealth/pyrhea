#! /usr/bin/env python

_rhea_svn_id_ = "$Id$"

from mpi4py import MPI
import numpy as np


class MsgTypes():
    GATE = 0
    ENDOFDAY = 1
    ANTI_ENDOFDAY = 2
    SINGLE_ENDOFDAY = 3
    SINGLE_ANTI_ENDOFDAY = 4


def getCommWorld():
    """Provide easy access to the world to packages that don't want to know about MPI"""
    return MPI.COMM_WORLD


class VectorClock(object):
    def __init__(self, commSize, rank):
        self.rank = rank
        self.vec = np.zeros(commSize, dtype=np.int32)

    def incr(self):
        self.vec[self.rank] += 1

    def merge(self, foreignVec):
        """ This operation does not include incrememting the local time """
        self.vec = np.maximum(self.vec, foreignVec)

    def max(self):
        return np.amax(self.vec)

    def min(self):
        return np.amin(self.vec)


class GblAddr(object):
    def __init__(self, rank, lclId):
        self.rank = rank
        self.lclId = lclId

    def __str__(self):
        return '%d_%d' % (self.rank, self.lclId)


class NetworkInterface(object):
    def __init__(self, comm, sync=True):
        self.comm = comm
        self.vclock = VectorClock(self.comm.size, self.comm.rank)
        self.outgoingDict = {}
        self.outstandingSendReqs = []
        self.outstandingRecvReqs = []
        self.expectFrom = set()
        self.clientIncomingCallbacks = {}
        self.sync = sync
        self.endOfDayMsg = []
        self.nEndOfDayNeighbors = 0
        self.incomingLclMessages = []

    def getGblAddr(self, lclId):
        return GblAddr(self.comm.rank, lclId)

    def getLclAddr(self, gblAddr):
        assert gblAddr.rank == self.comm.rank, 'Network object %s is not local!' % gblAddr
        return gblAddr.lclId

    def barrier(self):
        self.comm.Barrier()

    def enqueue(self, msgType, thing, srcAddr, gblAddr):
        toRank = gblAddr.rank
        if toRank not in self.outgoingDict:
            self.outgoingDict[toRank] = []
        self.outgoingDict[toRank].append((srcAddr, gblAddr, msgType, thing))

    def sendEOD(self, date):
        """Send an end-of-day message through all outgoing gates"""
        self.endOfDayMsg.append((MsgTypes.ENDOFDAY, date))

    def sendAntiEOD(self, date):
        """Send an revoke-end-of-day message through all outgoing gates"""
        self.endOfDayMsg.append((MsgTypes.ANTI_ENDOFDAY, date))

    def expect(self, srcAddr, destAddr, handleIncoming):
        """
        the handleIncoming is a callback with the signature

            handleIncoming(msgType, incomingTuple)

        There is no way to drop a rank from the expected source set because one can
        never be sure there is no straggler message from that rank
        """
        self.expectFrom.add(srcAddr.rank)
        assert destAddr.rank == self.comm.rank, "Cannot deliver to foreign object %s" % destAddr
        self.clientIncomingCallbacks[(srcAddr.rank, srcAddr.lclId,
                                      destAddr.lclId)] = handleIncoming

    def startRecv(self):
        for srcRank in self.expectFrom:
            self.outstandingRecvReqs.append(self.comm.irecv(None, srcRank, MPI.ANY_TAG))

    def _innerRecv(self, tpl):
        msgType, srcTag, destTag, partTpl = tpl
        # print ('######## %s: got %d agents for %s, time=%s, vtime=%s' %
        #        (self.name, len(agentList), destTag, tm, vtm))
        self.clientIncomingCallbacks[(srcTag.rank, srcTag.lclId, destTag.lclId)](msgType, partTpl)

    def finishRecv(self):
        self.vclock.incr()  # must happen before incoming messages arrive
        for tpl in self.incomingLclMessages:
            self._innerRecv(tpl)
        self.incomingLclMessages = []
        while True:
            if not self.outstandingRecvReqs:
                # print '######## %s nothing to recv' % self.name
                break
            # print ('######## %s entering recv testany on %d outstanding reqs' %
            #        (self.name, len(self.outstandingRecvReqs)))
            idx, flag, msg = MPI.Request.testany(self.outstandingRecvReqs)
            if idx >= 0:
                self.outstandingRecvReqs.pop(idx)
            if flag:
                vtm = msg[0]
                #
                # Handle vtime order issues here
                #
                self.vclock.merge(vtm)
                for tpl in msg[1:]:
                    self._innerRecv(tpl)
            else:
                # print '######## %s empty recv queue' % self.name
                break
        self.outstandingRecvReqs = []
        if self.sync:
            self.comm.Barrier()

    def startSend(self):
        vTimeNow = self.vclock.vec
        for destRank, msgList in self.outgoingDict.items():
            if destRank == self.comm.rank:
                # local message
                for srcTag, destTag, msgType, cargo in msgList:
                    self.incomingLclMessages.append((msgType, srcTag, destTag, cargo))
            else:
                bigCargo = [vTimeNow]
                for srcTag, destTag, msgType, cargo in msgList:
                    bigCargo.append((msgType, srcTag, destTag, cargo))
                # print '######### %s sent %s to %s' % (self.name, len(bigCargo), destRank)
                bigCargo.extend(self.endOfDayMsg)
                self.outstandingSendReqs.append(self.comm.isend(bigCargo, destRank))
        self.outgoingDict.clear()
        self.endOfDayMsg = []

    def finishSend(self):
        # print ('######## %s entering send waitall on %d requests' %
        #        (self.name, len(self.outstandingSendReqs)))
        result = MPI.Request.waitall(self.outstandingSendReqs)  # @UnusedVariable
        # print '######## %s finished send waitall; result was %s' % (self.name, result)
        self.outstandingSendReqs = []
