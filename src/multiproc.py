#!/usr/bin/env python

import pysam
import sys
import pybedtools
import os
import subprocess
from random import random
import random
import argparse
from uuid import uuid4
from re import sub
from itertools import izip
import commands

import vcf
import gzip
import shutil
import traceback
import time
import multiprocessing
from multiprocessing import Pool
from contextlib import closing
from pathos.multiprocessing import ProcessingPool 
import signal
import itertools

#handling concurrency
import logging.handlers
from functools import partial
import multiprocessing, threading, logging, sys, traceback,  StringIO, Queue
from functools import partial
from itertools import chain
from threading import Thread
import fnmatch
import settings


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
            dummy = self.format(record) # just to get traceback text into record.exc_text
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
        time.sleep(self.polltime+1) # give some time for messages to enter the queue.
        self.shutdown = True
        time.sleep(self.polltime+1) # give some time for the server to time out and see the shutdown

    def __del__(self):
        self.close() # hopefully this aids in orderly shutdown when things are going poorly.

