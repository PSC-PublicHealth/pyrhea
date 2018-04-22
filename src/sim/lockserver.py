#! /usr/bin/env python

###################################################################################
# Copyright   2018, Pittsburgh Supercomputing Center (PSC).  All Rights Reserved. #
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

import select
import socket
import signal
import random
import sys
import os
import logging
import optparse

#############################################################
#
# A couple of helper functions from stack overflow
#
#############################################################

LOGGER = logging.getLogger(__name__)

# from https://stackoverflow.com/questions/6229073/how-to-make-a-python-dictionary-that-returns-key-for-keys-missing-from-the-dicti
_DEFAULT_DICT_ID = 0
'''
Unique arbitrary identifier with which to uniquify the classname of the next
:func:`DefaultDict`-derived type.
'''

def DefaultDict(keygen):
    '''
    Sane **default dictionary** (i.e., dictionary implicitly mapping a missing
    key to the value returned by a caller-defined callable passed both this
    dictionary and that key).

    The standard :class:`collections.defaultdict` class is sadly insane,
    requiring the caller-defined callable accept *no* arguments. This
    non-standard alternative requires this callable accept two arguments:

    #. The current instance of this dictionary.
    #. The current missing key to generate a default value for.

    Parameters
    ----------
    keygen : CallableTypes
        Callable (e.g., function, lambda, method) called to generate the default
        value for a "missing" (i.e., undefined) key on the first attempt to
        access that key, passed first this dictionary and then this key and
        returning this value. This callable should have a signature resembling:
        ``def keygen(self: DefaultDict, missing_key: object) -> object``.
        Equivalently, this callable should have the exact same signature as that
        of the optional :meth:`dict.__missing__` method.

    Returns
    ----------
    MappingType
        Empty default dictionary creating missing keys via this callable.
    '''

    # Global variable modified below.
    global _DEFAULT_DICT_ID

    # Unique classname suffixed by this identifier.
    default_dict_class_name = 'DefaultDict' + str(_DEFAULT_DICT_ID)

    # Increment this identifier to preserve uniqueness.
    _DEFAULT_DICT_ID += 1

    # Dynamically generated default dictionary class specific to this callable.
    default_dict_class = type(
        default_dict_class_name, (dict,), {'__missing__': keygen,})

    # Instantiate and return the first and only instance of this class.
    return default_dict_class()


def getIp():
    """
    gets the ip addr of the default route (or at least the route to the IP addr in
    the connect statement). Doesn't send any packets, just "connects" a datagram socket

    stolen from a stackoverflow response
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        # doesn't have to be reachable
        if 'SLURM_LAUNCH_NODE_IPADDR' in os.environ:
            target = os.environ['SLURM_LAUNCH_NODE_IPADDR']
        else:
            target = '10.255.255.255'
        s.connect((target, 1))
        IP = s.getsockname()[0]
    except:
        IP = '127.0.0.1'
    finally:
        s.close()
        LOGGER.info('Using IP addr %s', IP)
    return IP


################################################################
#
# This section is the server code
#
# LockServer() is the main class to run the server
# At the bottom of this file, this is called from the main
# hook if this file is run standalone rather than imported.
#
################################################################

class ServerLock(object):
    def __init__(self, lName):
        LockDict[lName] = self
        self.lName = lName
        self.clientsWaiting = []
        self.accessCount = 0

    def tryRequest(self, shared):
        """
        tries to request this lock either shared or exclusive

        returns True if successful or False otherwise
        """
        if self.accessCount == 0:
            self.shared = shared
            self.accessCount += 1
            return True

        if self.shared and shared:
            self.accessCount += 1
            return True

        return False

    def request(self, shared, client, wait):
        """
        requests lock asynchronously,
        returns True if lock is available immediately otherwise False
        if wait is set, adds client to the wait list and the client will be notified when the lock
        is given to them.
        """
        if self.tryRequest(shared):
            return True

        if wait:
            self.clientsWaiting.append((client,shared))
        return False

    def release(self):
        self.accessCount -= 1

        if self.accessCount > 0:
            return

        if not self.clientsWaiting:
            del LockDict[self.lName]
            return

        client, shared = self.clientsWaiting.pop(0)
        self.tryRequest(shared)  # this should be guaranteed to work
        client.notify()

        if not shared:  # if it's an exclusive lock, there's nothing more to do
            return

        # if it's a shared lock, check the remaining clients and see which others want access
        clients = self.clientsWaiting
        self.clientsWaiting = []
        for client,shared in clients:
            if shared:
                self.tryRequest(shared)  # again this should always work
                client.notify()
            else:
                self.clientsWaiting.append((client, shared))

    def clearWaiting(self, client, shared):
        self.clientsWaiting.remove((client, shared))


class LockConnection(object):
    def __init__(self, clientSocket, address):
        self.fileno = clientSocket.fileno()
        ClientDict[self.fileno] = self
        self.clientSocket = clientSocket
        self.address = address
        LOGGER.debug("new client from %s on fd %s", address, self.fileno)
        self.readBuf = ""
        self.locks = {}  # True if shared, False if exclusive
        self.waiting = False  # True if waiting for a lock

        self.clientSocket.setblocking(0)

    def read(self):
        try:
            data = self.clientSocket.recv(select.PIPE_BUF)
            if len(data) == 0:
                self.killClient()
                return
            self.readBuf += data
        except:
            LOGGER.warning("got exception on read")
            self.killClient()
            return
        self.processReadBuf()

    def processReadBuf(self):

        # process read buffer for as long as we have complete commands and aren't waiting for a lock
        while self.waiting is False:
            line, nl, rest = self.readBuf.partition('\n')
            if nl != '\n':
                return
            self.readBuf = rest
            line = line.strip()
            cmd, sp, lName = line.partition(' ')
            if sp != ' ':
                self.killClient()
                return

            lName = lName.strip()
            if cmd == "release":
                if lName not in self.locks:
                    self.send("ERROR %s not already locked\n"%lName)
                    continue
                LockDict[lName].release()
                del self.locks[lName]
                self.send("RELEASED %s\n"%lName)
                continue

            # all other commands require an lName that isn't already in locks
            if lName in self.locks:
                self.send("ERROR %s already locked\n"%lName)
                continue

            if cmd == "xlock":
                self.request(lName, shared=False, wait=False)
            elif cmd == "slock":
                self.request(lName, shared=True, wait=False)
            elif cmd == "xlockwait":
                self.request(lName, shared=False, wait=True)
            elif cmd == "slockwait":
                self.request(lName, shared=True, wait=True)
            else:
                self.killClient()
                return

    def request(self, lName, shared, wait):
        if LockDict[lName].request(shared=shared, client=self, wait=wait):
            self.locks[lName] = shared
            self.send("ACQUIRED %s\n"%lName)
            return

        if wait:
            self.waitName = lName
            self.waitShared = shared
            self.waiting = True
            return

        self.send("FAILED %s\n"%lName)

    def notify(self):
        self.locks[self.waitName] = self.waitShared
        self.send("ACQUIRED %s\n"%self.waitName)
        self.waiting = False

    def send(self, data):
        """
        send data back to the client.  For now I'm going to say that pipelining requests is
        not legal and if a send would block then the client is behaving badly and should die

        this will fail on unicode strings
        """
        if len(data) != self.clientSocket.send(data):
            self.killClient()

    def killClient(self):
        for lName,share in self.locks.items():
            LockDict[lName].release()

        if self.waiting:
            LockDict[self.waitName].clearWaiting(self, self.waitShared)

        try:
            self.clientSocket.close()
        except:
            pass

        del ClientDict[self.fileno]
        LOGGER.debug('killClient on fd %s', self.fileno)

    def socketError(self):
        self.killClient()



class LockServer(object):
    def __init__(self, host='', port=29292):
        self.host = host
        self.port = port
        self.run = True

    def serve(self):
        self.listenSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.listenSocket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.listenSocket.bind((self.host, self.port))
        self.listenSocket.listen(5)

        ClientDict[self.listenSocket.fileno()] = self

        while self.run:
            try:
                readList, writeList, exceptList = select.select(ClientDict.keys(), # @UnusedVariable
                                                                [],
                                                                ClientDict.keys())
                for e in exceptList:
                    # we're dealing with some error on the socket - just deal with it and stay alive
                    try:
                        ClientDict[e].socketError()
                    except Exception, ee:  # @UnusedVariable
                        pass

                for r in readList:
                    if r in ClientDict:
    #                try:  # it's extremely possible that r has been expunged from ClientDict by now
                        ClientDict[r].read()
    #                except:
    #                    pass
            except select.error, se:  # @UnusedVariable
                pass  # select throws an error when we catch a signal

    def shutdown(self):
        """
        Trigger exit of the main loop
        """
        self.run = False
        # the signal that got us here will have broken us out of select

    def cleanup(self):
        """
        Perform cleanup actions
        """
        pass

    def read(self):
        newSock, addr = self.listenSocket.accept()
        LockConnection(newSock, addr)

    def socketError(self):
        # something has gone to hell with the server.  Just shut down
        sys.exit()


LockDict = DefaultDict(lambda dd,key: ServerLock(key)) # key is name of lock, val is ServerLock
ClientDict = {} # key is fileno of clientsocket, val is LockConnection

################################################################
#
# This section is to set up a few defaults and make the server
# easily findable in many environments.
#
################################################################

DefaultLockHost = 'localhost'
DefaultLockPort = 29292
DefaultLockClient = None
DefaultLockFile = 'lockserver.info'

def setDefaultLockServer(host, port):
    global DefaultLockHost
    global DefaultLockPort

    DefaultLockHost = host
    DefaultLockPort = port

def _tryPort(port):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = False
    try:
        sock.bind(("0.0.0.0", port))
        result = True
    except:
        LOGGER.info("Port %s is in use", port)
    sock.close()
    return result

def findFreePort():
    port = DefaultLockPort
    portRange = (DefaultLockPort, DefaultLockPort + 100)
    maxTries = portRange[1] - portRange[0]
    for i in range(maxTries):
        if _tryPort(port):
            break
        port = random.randint(*portRange)
    else:
        sys.exit('No ports available between %s and %s' % portRange)
    return port

def discoverLockServer(fileName = None):
    host = DefaultLockHost
    port = DefaultLockPort
    pid = None

    if fileName is None:
        fileName = DefaultLockFile

    with open(fileName) as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        if line.startswith('#'):
            continue
        key,eq,val = line.partition('=')
        if eq != '=':
            raise RuntimeError("invalid line %s in lockfile pointer %s"%(line, fileName))
        key = key.strip()
        val = val.strip()
        if key == 'host':
            host = val
        elif key == 'port':
            port = int(val)
        elif key == 'pid':
            pid = int(val)

    return host, port, pid


def writeLockServerHostFile(fileName=None, host=None, port=None):
    if fileName is None:
        fileName = DefaultLockFile
    if host is None:
        host = getIp()
    if port is None:
        port = DefaultLockPort

    with open(fileName, "w") as f:
        f.write("host = %s\n"%host)
        f.write("port = %s\n"%port)
        f.write("pid = %s\n"%os.getpid())

    return fileName


def rmLockServerHostFile(fileName=None):
    if fileName is None:
        fileName = DefaultLockFile


################################################################
#
# This is the client code
#
# use Lock() and SharedLock() as the main classes and it's
# easiest to use them in the form of a context manager ie:
#
# with Lock(someNameStr):
#     some code that needs exclusive access
#
#   or
#
# with SharedLock(someNameStr):
#     some code that needs shared lock access
#
################################################################


class LockClient(object):
    "handles communications with lock server.  Not intended to be directly used by user"
    def __init__(self, host=None, port=None):
        if port is None:
            port = DefaultLockPort
        if host is None:
            host = DefaultLockHost
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((host, port))


    def getLock(self, cmd, lName):
        "issue any of the lock commands, responses are same"
        req = cmd + " " + lName + "\n"
        self.sock.send(req)
        s = ""
        while "\n" not in s:
            s += self.sock.recv(select.PIPE_BUF)

        assert s.endswith('\n'), "invalid strings coming from lock server"
        s = s.rstrip()
        result, sp, msg = s.partition(' ')
        if result == "ACQUIRED":
            return True
        if result == "FAILED":
            return False
        # anything else is an error
        raise RuntimeError("invalid respone from lockserver")

    def releaseLock(self, lName):
        req = "release %s\n"%lName
        self.sock.send(req)
        s = ""
        while "\n" not in s:
            s += self.sock.recv(select.PIPE_BUF)

        assert s.endswith('\n'), "invalid strings coming from lock server"
        s = s.rstrip()
        result, sp, msg = s.partition(' ')

        assert result=="RELEASED", "invalid response to release from lock server"
        return True

    def close(self):
        try:
            self.sock.close()
        except:
            pass

class Lock(object):
    """
    This is our basic client lock class.  If used with defaults it will read the
    """
    cmdDict = {(False, False): "xlock",
               (False, True): "xlockwait",
               (True, False): "slock",
               (True, True): "slockwait"}

    def __init__(self, lockName, shared=False, wait=True, discoverServer=True,
                 host=None, port=None, lockClient=None):
        """
        wait can't be set to False if using this as a context manager
        """
        self.lName = lockName
        self.shared = shared
        self.wait = wait
        self.host = host
        self.port = port
        self.discoverServer=discoverServer
        self.lockClient = lockClient

    def __enter__(self):
        if self.wait == False:
            raise RuntimeError("Can't use Lock() in context manager with wait set to False")

        if not self.lock():
            raise RuntimeError("Failed to acquire lock %s"%self.lName)

        LOGGER.debug('enter Lock %s %s', self.lName, ('shared' if self.shared else ''))

    def __exit__(self, *args):
        self.release()
        LOGGER.debug('exit Lock %s %s', self.lName, ('shared' if self.shared else ''))

    def lock(self):
        global DefaultLockClient

        if self.lockClient is None:
            if DefaultLockClient is None:
                if self.discoverServer is False:
                    DefaultLockClient = LockClient(self.host, self.port)
                else:
                    if self.discoverServer is True:
                        filename = None
                    else:
                        filename = self.discoverServer
                    self.host, self.port, _pid =  discoverLockServer(filename)
                    DefaultLockClient = LockClient(self.host, self.port)

            self.lockClient = DefaultLockClient
        cmd = Lock.cmdDict[(self.shared, self.wait)]
        return self.lockClient.getLock(cmd, self.lName)

    def release(self):
        return self.lockClient.releaseLock(self.lName)


class SharedLock(Lock):
    def __init__(self, *args, **kwargs):
        kwargs['shared']=True
        super(SharedLock, self).__init__(*args, **kwargs)


################################################
#
#  main hook
#
################################################

def main():
    parser = optparse.OptionParser(usage="""
    %prog [--debug]
    """)
    parser.add_option("-d", "--debug", action="store_true",
                      help="debugging output")
    opts, args = parser.parse_args()
    if args:
        parser.error('Invalid arguments %s' % args)
    if opts.debug:
        logLvl = 'debug'
    else:
        logLvl = 'info'
    logging.basicConfig(level=getattr(logging, logLvl.upper(), None))
    try:
        host, port, pid = discoverLockServer()
        if (host is not None and port is not None and pid is not None
                and pid != os.getpid()):
            # Valid file; now check for running lockserver
            os.kill(pid, 0) # Throws an error if that pid is not live, no effect if it is
            # OK, so there is already a live server pointed to by this file
            sys.exit('A running lock server is already present')
    except Exception, e:
        pass # No old file or server
    port = findFreePort()
    LOGGER.info('Using port %s', port)

    hostFileName = writeLockServerHostFile(port=port)
    lockServer = LockServer(port=port)

    def handle_signal(sig, frame):
        lockServer.shutdown()

    if os.name != "nt":
        signal.signal(signal.SIGTERM, handle_signal)

    try:
        lockServer.serve()  # should never exit
    finally:
        lockServer.cleanup()
        os.remove(hostFileName)
        LOGGER.info('exiting')

if __name__ == "__main__":
    main()

