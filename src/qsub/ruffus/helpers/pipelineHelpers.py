import re
import os
import getpass
import subprocess
from ruffus.proxy_logger import *
import logging
import time
import random
from time import gmtime, strftime
from helpers import runIDHelpers as rid
from helpers import parameters as params
import helpers
from re import sub
from ruffus import *

configReader = params.GetConfigReader()

def GetLogFile(logger_name):
    """returns the log file"""
    
    results_path = configReader.get('RESULTS', 'results_path')    
    cancer_type = params.GetCancerType() 
    cancer_dir_path = "/".join([results_path, cancer_type])
    
    log = []
    
    
    log_file_path = cancer_dir_path + '/' + \
        params.GetProjectName(
        ) + '_' + rid.GetRunID() + '/batch_logs/'
    log_file = log_file_path + 'pipeline_batch_log_' + \
        params.GetProjectName() + '_' + rid.GetRunID() + '.log'

    try:
        os.makedirs(log_file_path)
    except:
        pass

    logger_args = {}
    logger_args["file_name"] = log_file
    logger_args["level"] = logging.DEBUG
    logger_args["rotating"] = True
    logger_args["maxBytes"] = 10000000
    logger_args["backupCount"] = 10
    logger_args[
        "formatter"] = "[%(asctime)s] [%(name)s] [%(levelname)s] - %(message)s"

    logger_proxy, logger_mutex = make_shared_logger_and_proxy(
        setup_std_shared_logger, logger_name, logger_args)
    log.append(logger_proxy)
    log.append(logger_mutex)

    return log


def Logging(log_type, log, message):
    """write to specified log"""
    logger_proxy = log[0]
    logger_mutex = log[1]

    if log_type == 'INFO':
        with logger_mutex:
            logger_proxy.info(message)
    elif log_type == 'WARNING':
        with logger_mutex:
            logger_proxy.warning(message)
    elif log_type == 'ERROR':
        with logger_mutex:
            logger_proxy.error(message)


def LogParameters(logger_name):
    """Logs the parameters used in the run"""
    log = GetLogFile(logger_name)
    Logging('INFO', log, '-------------------------' +
            logger_name + ' CONFIG FILE ARGUMENTS-------------------------')

    for option in configReader.options(logger_name):
        # logger_name is used as the section name in the log file
        Logging('INFO', log, '[{0}] {1}: {2}'.format(
            logger_name, option, configReader.get(logger_name, option)))


def LogInfile(infile):
    """Logs the contents of the infile"""
    log = GetLogFile('MAIN')
    Logging('INFO', log, '-------------------------INPUT FILE-------------------------')
    for line in infile:
        Logging('INFO', log, line)


def LogArguments(argdict, version):
    """Logs the command line arguments, program version, date"""

    log = GetLogFile('MAIN')
    Logging('INFO', log, '\n-----------------------------------------------\nPipeline Start:' \
        '\n-----------------------------------------------\nPROGRAM: ' + version)
    Logging('INFO', log, 'RUNDATE: ' +
            strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    Logging('INFO', log, 'USER: ' + getpass.getuser())
    Logging('INFO', log, '-------------------------COMMAND LINE ARGUMENTS-------------------------')
    print '--------------------------------------------------------'
    print 'Pipeline invoked with the following arguments:'
    print '--------------------------------------------------------'
    for a, v in argdict:
        line = a + ' = ' + str(v)
        print line
        Logging('INFO', log, '\t' + line)


def LogInputCheck(msg):
    """Logs the results of InputChecker.py module"""

    # to screen
    hdr = '\n-----------------------------------------------\nChecking Input:' \
          '\n-----------------------------------------------\n'
    footer = '\n-----------------------------------------------\n'
    print hdr + msg + footer

    # to log
    log = GetLogFile('MAIN')
    Logging('INFO', log, hdr + msg)


def runCommand(cmd):
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
    except OSError as e:
        print("Execution failed:", e)


def RunTask(command, num_cpu, mem_usage, sample_id, software):
    """assignes a task to the cluster"""     
    results_path = configReader.get('RESULTS', 'results_path')    
    cancer_type = params.GetCancerType() 
    cancer_dir_path = "/".join([results_path, cancer_type])
    
    time.sleep(random.uniform(params.num_samples, params.num_samples + 20))
    log_path = cancer_dir_path + '/' + params.GetProjectName() + '_' + rid.GetRunID()  +  '/logs/' 
    try:
        os.makedirs(log_path)
    except:
        pass
    cluster_cmd = params.GetQsubStatement().format(
        log_path, str(num_cpu), mem_usage)
    process = []
    subprocess.call('chmod +x ' + command, shell=True)
    cluster_cmd = cluster_cmd + ' ' + command

    max_jobs = int(configReader.get('CLUSTER', 'max_jobs'))
    current = 0

    jobcheck = subprocess.Popen('qstat -u {0} | wc -l'.format(
        getpass.getuser()), stdout=subprocess.PIPE, shell=True)
    for line in jobcheck.stdout:
        current = int(re.sub(r'[\r\n]', '', line)) - 2

    while current >= max_jobs:
        time.sleep(10)
        jobcheck2 = subprocess.Popen('qstat -u {0} | wc -l'.format(
            getpass.getuser()), stdout=subprocess.PIPE, shell=True)
        for line in jobcheck2.stdout:
            current = int(re.sub(r'[\r\n]', '', line)) - 2
    
    
    task = subprocess.Popen(cluster_cmd, shell=True)
    process.append(task)
    process.append(cluster_cmd)
    process.append(0)
    process.append(log_path + os.path.basename(command))

    return process


def GetScriptPath(sample_id, software):
    """returns path to store scripts"""
    
    results_path = configReader.get('RESULTS', 'results_path')
    cancer_type = params.GetCancerType() 
    cancer_dir_path = "/".join([results_path, cancer_type])
    script_path = cancer_dir_path + '/' + params.GetProjectName() + '_' + rid.GetRunID(
    )  + '/scripts/'

    try:
        os.makedirs(script_path)
    except:
        pass
    return script_path


def CheckSentinel(sentinel_list, log, log_msg):
    """Checks if sentinel file exists"""
    if len(sentinel_list) != 0:
        for sentinel_file in sentinel_list:
            if sentinel_file[0] != '':
                if os.path.isfile(sentinel_file[0]) != True:
                    Logging('ERROR', log, log_msg + 'Could not find previous file ' +
                            sentinel_file[0] + '. Ignoring current task...')
                    return False
    return True


def CheckTaskStatus(taskList, output_sentinel, log, log_msg):
    """Checks which tasks on the queue are done"""

    failed = False
    while len(taskList) > 0:
        for task in taskList:
            if task[0].poll() != None:
                   # Error code 192 is the exit code for one of the programs in
                   # the pipeline
                if task[0].returncode > 0 and task[0].returncode != 192:
                    if task[2] == 0:
                        task[0] = subprocess.Popen(task[1], shell=True)
                        task[2] += 1
                        Logging('WARNING', log, log_msg + 'Attempt#' + str(
                            task[2]) + ' failed, retrying: ' + task[1] + ' Error Code: ' \
                                + str(task[0].returncode))
                        Logging('WARNING', log, log_msg + 'See: ' + task[3])
                    else:
                        Logging(
                            'ERROR', log, log_msg + 'Failed to run: ' + task[1])
                        Logging('ERROR', log, log_msg + 'See: ' + task[3])
                        failed = True
                        taskList.remove(task)
                else:
                    taskList.remove(task)
            time.sleep(10)

    if failed:
        Logging('ERROR', log, log_msg + 'Completed with errors')
    else:
        Logging('INFO', log, log_msg + 'Completed Sucessfully!')
        subprocess.call('touch ' + output_sentinel, shell=True)
    return not failed 

