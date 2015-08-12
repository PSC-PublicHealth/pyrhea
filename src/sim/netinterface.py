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

_rhea_svn_id_ = "$Id$"

from mpi4py import MPI
import numpy as np
import types
from collections import namedtuple, deque


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


_InnerGblAddr = namedtuple('_innerGblAddr', ['rank', 'lclId'])


class GblAddr(_InnerGblAddr):

    def getLclAddr(self):
        return self.lclId

    def getPatchAddr(self):
        if isinstance(self.lclId, types.TupleType):
            return GblAddr(self.rank, self.lclId[0])
        else:
            return GblAddr(self.rank, self.lclId)

    def __str__(self):
        if isinstance(self.lclId, types.TupleType):
            return '%d_%d_%d' % (self.rank, self.lclId[0], self.lclId[1])
        else:
            return '%d_%d' % (self.rank, self.lclId)

    def __lt__(self, other):
        return (self.rank < other.rank
                or (self.rank == other.rank and self.lclId < other.lclId))

    def __le__(self, other):
        return self < other or self == other

    def __eq__(self, other):
        return self.rank == other.rank and self.lclId == other.lclId

    def __ne__(self, other):
        return self.rank != other.rank or self.lclId != other.lclId

    def __gt__(self, other):
        return (self.rank > other.rank
                or (self.rank == other.rank and self.lclId > other.lclId))

    def __ge__(self, other):
        return self > other or self == other


class NetworkInterface(object):
    MPI_TAG_MORE = 1
    MPI_TAG_END = 2
#     maxChunksPerMsg = 256000
#     irecvBufferSize = 1024*1024
    maxChunksPerMsg = 32
    irecvBufferSize = 1024*1024

    def __init__(self, comm, sync=True, deterministic=False):
        self.comm = comm
        self.vclock = VectorClock(self.comm.size, self.comm.rank)
        self.outgoingDict = {}
#         self.outstandingSendReqs = deque()
#         self.outstandingRecvReqs = deque()
        self.outstandingSendReqs = []
        self.outstandingRecvReqs = []
        self.expectFrom = set()
        self.clientIncomingCallbacks = {}
        self.sync = sync
        self.deterministic = deterministic
        self.endOfDayMsg = []
        self.nEndOfDayNeighbors = 0
        self.incomingLclMessages = []

    def getGblAddr(self, lclId):
        return GblAddr(self.comm.rank, lclId)

    def isLocal(self, gblAddr):
        return gblAddr.rank == self.comm.rank

    def barrier(self):
        self.comm.Barrier()

    def enqueue(self, msgType, thing, srcAddr, gblAddr):
        toRank = gblAddr.rank
        if toRank not in self.outgoingDict:
            self.outgoingDict[toRank] = []
        self.outgoingDict[toRank].append((srcAddr, gblAddr, msgType, thing))

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
        if self.deterministic:
            l = [a for a in self.expectFrom]
            l.sort()
            for srcRank in l:
                if srcRank != self.comm.rank:
                    buf = bytearray(NetworkInterface.irecvBufferSize)
                    self.outstandingRecvReqs.append(self.comm.irecv(buf, srcRank, MPI.ANY_TAG))
        else:
            for srcRank in self.expectFrom:
                if srcRank != self.comm.rank:
                    buf = bytearray(NetworkInterface.irecvBufferSize)
                    self.outstandingRecvReqs.append(self.comm.irecv(buf, srcRank, MPI.ANY_TAG))

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
                break
            if self.deterministic:
                s = MPI.Status()
                msg = MPI.Request.wait(self.outstandingRecvReqs[0], s)
                self.outstandingRecvReqs.popleft()
                tag = s.Get_tag()
                if tag == NetworkInterface.MPI_TAG_MORE:
                    print 'MORE from %s' % s.Get_source()
                    buf = bytearray(NetworkInterface.irecvBufferSize)
                    self.outstandingRecvReqs.append(self.comm.irecv(buf, s.Get_source(),
                                                                    MPI.ANY_TAG))
                vtm = msg[0]
                #
                # Handle vtime order issues here
                #
                self.vclock.merge(vtm)
                for tpl in msg[1:]:
                    self._innerRecv(tpl)
            elif self.sync:
                s = MPI.Status()
                idx, msg = MPI.Request.waitany(self.outstandingRecvReqs, s)
                self.outstandingRecvReqs.pop(idx)
                tag = s.Get_tag()
                if tag == NetworkInterface.MPI_TAG_MORE:
                    print 'MORE from %s' % s.Get_source()
                    buf = bytearray(NetworkInterface.irecvBufferSize)
                    self.outstandingRecvReqs.append(self.comm.irecv(buf, s.Get_source(),
                                                                    MPI.ANY_TAG))
                vtm = msg[0]
                #
                # Handle vtime order issues here
                #
                self.vclock.merge(vtm)
                for tpl in msg[1:]:
                    self._innerRecv(tpl)
            else:
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
#         self.outstandingRecvReqs = deque()
        self.outstandingRecvReqs = []

    def startSend(self):
        vTimeNow = self.vclock.vec
        if self.deterministic:
            l = self.outgoingDict.keys()
            l.sort()
            for destRank in l:
                msgList = self.outgoingDict[destRank][:]
                msgList.sort()
                if destRank == self.comm.rank:
                    # local message
                    for srcTag, destTag, msgType, cargo in msgList:
                        self.incomingLclMessages.append((msgType, srcTag, destTag, cargo))
                else:
                    bigCargo = [vTimeNow]
                    for srcTag, destTag, msgType, cargo \
                            in msgList[0:NetworkInterface.maxChunksPerMsg]:
                        bigCargo.append((msgType, srcTag, destTag, cargo))
                    msgList = msgList[NetworkInterface.maxChunksPerMsg:]
                    if msgList:
                        req = self.comm.isend(bigCargo, destRank,
                                              tag=NetworkInterface.MPI_TAG_MORE)
                    else:
                        bigCargo.extend(self.endOfDayMsg)
                        req = self.comm.isend(bigCargo, destRank,
                                              tag=NetworkInterface.MPI_TAG_END)
                    self.outstandingSendReqs.append(req)

        else:
            for destRank, msgList in self.outgoingDict.items():
                if destRank == self.comm.rank:
                    # local message
                    for srcTag, destTag, msgType, cargo in msgList:
                        self.incomingLclMessages.append((msgType, srcTag, destTag, cargo))
                else:
                    while msgList:
                        bigCargo = [vTimeNow]
                        for srcTag, destTag, msgType, cargo \
                                in msgList[0:NetworkInterface.maxChunksPerMsg]:
                            bigCargo.append((msgType, srcTag, destTag, cargo))
                        msgList = msgList[NetworkInterface.maxChunksPerMsg:]
                        if msgList:
                            req = self.comm.isend(bigCargo, destRank,
                                                  tag=NetworkInterface.MPI_TAG_MORE)
                        else:
                            bigCargo.extend(self.endOfDayMsg)
                            req = self.comm.isend(bigCargo, destRank,
                                                  tag=NetworkInterface.MPI_TAG_END)
                        self.outstandingSendReqs.append(req)
                    # print '######### %s sent %s to %s' % (self.name, len(bigCargo), destRank)
        self.outgoingDict.clear()
        self.endOfDayMsg = []

    def finishSend(self):
        # print ('######## %s entering send waitall on %d requests' %
        #        (self.name, len(self.outstandingSendReqs)))
        result = MPI.Request.waitall(self.outstandingSendReqs)  # @UnusedVariable
        # print '######## %s finished send waitall; result was %s' % (self.name, result)
#         self.outstandingSendReqs = deque()
        self.outstandingSendReqs = []
