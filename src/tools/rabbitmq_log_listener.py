#!/usr/bin/env python
import pika
import pickle
import logging
from optparse import OptionParser

logger = logging.getLogger(__name__)


def callback(ch, method, properties, pickledBody):
    bodyDict = pickle.loads(pickledBody)
    # print bodyDict
    numLevel = getattr(logging, bodyDict['level'])
    extras = {newNm: bodyDict[oldNm] for oldNm, newNm in [('level', 'srclevel'),
                                                          ('pathname', 'srcpathname'),
                                                          ('lineno', 'srclineno')]}
    logger.log(numLevel, bodyDict['message'], extra=extras)


def main():
    parser = OptionParser(usage="""
    %prog [-v][-d]
    """)
    parser.add_option("-v", "--verbose", action="store_true",
                      help="show INFO messages and above")
    parser.add_option("-d", "--debug", action="store_true",
                      help="show DEBUG messages and above")
    parser.add_option("-l", "--line", action="store_true",
                      help="show originating line number")

    opts, args = parser.parse_args()
    if args:
        parser.error('No arguments expected')

    if opts.verbose:
        logLevel = 'INFO'
    elif opts.debug:
        logLevel = 'DEBUG'
    else:
        logLevel = 'WARNING'
    if opts.line:
        fmt = "%(levelname)s: %(srcpathname)s:%(srclineno)s: %(message)s"
    else:
        fmt = "%(levelname)s %(message)s"
    parser.destroy()

    connection = pika.BlockingConnection(pika.ConnectionParameters(host='localhost'))
    channel = connection.channel()

    channel.queue_declare(queue='pyrhea_logger')
    logging.basicConfig(level=logLevel, format=fmt)

    print 'Waiting for messages. To exit press CTRL+C'
    channel.basic_consume(callback, queue='pyrhea_logger', no_ack=True)
    channel.start_consuming()

    logging.shutdown()

############
# Main hook
############

if __name__ == "__main__":
    main()
