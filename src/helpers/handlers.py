import os
from helpers import parameters as params
import logging.handlers
#from pathos.multiprocessing import ProcessingPool
import multiprocessing
import time
from threading import Thread
import multiprocessing, threading, logging, sys, traceback, StringIO, Queue


def createDirectory(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def GetProjectPaths(results_path):
    cancer_type = params.GetCancerType()
    if (cancer_type):
        cancer_dir_path = "/".join([results_path, cancer_type])
    else:
        cancer_dir_path = results_path

    haplotype_path = "/".join([cancer_dir_path, "haplotypedir"])
    log_path = "/".join([cancer_dir_path, "logs"])
    tmpbams_path = "/".join([cancer_dir_path, "tmpbams"])
    finalbams_path = "/".join([cancer_dir_path, "finalbams"])
    logfile = "/".join([log_path, "debug.log"])

    createDirectory(results_path)
    createDirectory(cancer_dir_path)
    createDirectory(haplotype_path)
    createDirectory(log_path)
    createDirectory(tmpbams_path)
    createDirectory(finalbams_path)
    return (haplotype_path, cancer_dir_path, tmpbams_path, finalbams_path, log_path, logfile)


def GetLoggings(logfile):
    terminating = multiprocessing.Event()
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)
    logQueue = multiprocessing.Queue(16)
    filehandler = MultiProcessingLogHandler(logging.FileHandler(logfile), logQueue)
    logger.addHandler(filehandler)
    filehandler.setLevel(logging.DEBUG)
    return (terminating, logger, logQueue)


def GetHetBamPaths(tmpbams_path, chr, event):
    NHET = "/".join([tmpbams_path, str(chr + event).upper() + '_NH.bam'])
    HET = "/".join([tmpbams_path, str(chr + event).upper() + '_H.bam'])
    return (HET, NHET)


class MultiProcessingLogHandler(logging.Handler):
    def __init__(self, handler, queue, child=False):
        logging.Handler.__init__(self)

        self._handler = handler
        self.queue = queue

        if child == False:
            self.shutdown = False
            self.polltime = 1
            t = threading.Thread(target=self.receive)
            t.daemon = True
            t.start()

    def setFormatter(self, fmt):
        logging.Handler.setFormatter(self, fmt)
        self._handler.setFormatter(fmt)

    def receive(self):
        while (self.shutdown == False) or (self.queue.empty() == False):
            try:
                record = self.queue.get(True, self.polltime)
                self._handler.emit(record)
            except Queue.Empty, e:
                pass

    def send(self, s):
        self.queue.put(s)

    def _format_record(self, record):
        ei = record.exc_info
        if ei:
            dummy = self.format(record)  # just to get traceback text into record.exc_text
            record.exc_info = None  # to avoid Unpickleable error

        return record

    def emit(self, record):
        try:
            s = self._format_record(record)
            self.send(s)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

    def close(self):
        time.sleep(self.polltime + 1)  # give some time for messages to enter the queue.
        self.shutdown = True
        time.sleep(self.polltime + 1)  # give some time for the server to time out and see the shutdown

    def __del__(self):
        self.close()  # hopefully this aids in orderly shutdown when things are going poorly.


class bidir_itr(object):
    def __init__(self, collection):
        self.collection = collection
        self.index = 0

    def next(self):
        try:
            result = self.collection[self.index]
            self.index += 1
        except IndexError:
            raise StopIteration
        return result

    def prev(self):
        self.index -= 1
        if self.index < 0:
            raise StopIteration
        return self.collection[self.index]

    def __iter__(self):
        return self
