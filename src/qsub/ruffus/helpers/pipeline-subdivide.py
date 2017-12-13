import os
import sys
from re import sub
from ruffus import *
from helpers import depthOfCoverageTasks, depthOfCoverageHelpers, pipelineHelpers
from helpers import runIDHelpers as rid
from helpers import parameters as params
import subprocess 

current_path = params.GetProgramPath()
log = pipelineHelpers.GetLogFile('DepthOfCoverage')

        
def runCommand(cmd):
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
    except OSError as e:
        print("Execution failed:", e)