import argparse
import os
import time
from ruffus import *
from helpers import runIDHelpers as rid
from helpers import parameters as params
from helpers import inputChecker as ic
from helpers import pipelineHelpers
from helpers import taskHelpers
from helpers import bamgineerHelpers as bamhelp

__version__ = 'Bamgineer v0.1.0'

events= [['gain'],['loss']]
@parallel(events)
def parallel_events(evt):
    if(evt == 'gain'):
        pipeline_run([pipeline.subsample_gain], multiprocess=4)
       
    elif(evt == 'loss'):
        pipeline_run([pipeline.subsample_loss], multiprocess=2)
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='adds CN spikes to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-outbam', dest='outBamFile', required=True,
                        help='.bam file name for output')  
    parser.add_argument('-cnv_amp', dest='cnvAmpFile', required=False,
                        help='CNV amplification .bed file name')
    parser.add_argument('-cnv_del', dest='cnvDelFile', required=False,
                        help='CNV deletion .bed file name')
    parser.add_argument('-i', '--infile', action='store', dest='infile', required=True,
        help='2 column (tab delimited) input file: ID  tumour_path')
    parser.add_argument('-cancertype', dest='cancerType', required=False,
                        help='acronyms for cancer type')
    parser.add_argument('-splitbamdir', dest='splitbams', required=False,
                        help='input bam split by chromosomes')
    parser.add_argument('-c', '--configFile', action='store', required=True, dest='configfile',
                        help='path to config_file.cfg')
    parser.add_argument('-phase',dest= 'phase', action="store_true")
    parser.add_argument('-ctDNA', dest='ctDNA', action="store_true")
    
    parser.add_argument(
        '--project-name', action='store', required=True, dest='project_name',
        help='name of the project')
    parser.add_argument(
        '-r', '--run_id', action='store', type=int, dest='custom_run_id',
        help='manually enter run_id (used to re-launch a run)', default=False)
    
    args = parser.parse_args()
    t0 = time.time()
   
    args = parser.parse_args()
    current_path = os.path.dirname(os.path.realpath(__file__))
   
    outbamfn = args.outBamFile
    configReader = params.GetConfigReader()
   
    params.InitConfigReader(args.configfile)
    params.SetProjectName(args.project_name)
    params.SetProgramPath(current_path)
    params.SetInfile(args.infile)
    
    params.SetGainCNV(args.cnvAmpFile)
    params.SetLossCNV(args.cnvDelFile)
    params.SetCancerType(args.cancerType)
    params.SetOutputFileName(args.outBamFile)
    
 
    # set the run id used for this run
    if args.custom_run_id == False:
        # increment and set runID
        rid.IncRunID(params.GetProjectName(), current_path)
    else:
        # use runid specified on the command line -r
        rid.SetRunID(args.custom_run_id, current_path)

        # check that the user specified run direcotry exists
        if not os.path.exists(taskHelpers.GetProjectNamePathRunID()[0]):
            print (str(taskHelpers.GetProjectNamePathRunID())[0])
            raise Exception('The project folder for this run does not exist')
        
    from helpers import pipeline

    def LaunchPipeline():
        """initiates pipeline found in pipeline.py"""
        print 'Bamgineer COMPLETE'

    # log command line arguments
    pipelineHelpers.LogArguments(vars(args).items(), __version__)

    # log config file arguments
    pipelineHelpers.LogParameters('CLUSTER')

    ## parse the input file
    infile_list = taskHelpers.ParseInfile(args.infile)
    line_count = len(infile_list)
    params.num_samples = len(infile_list)
    
    ## log input bam file
    pipelineHelpers.LogInfile(infile_list)

    # verify inputs, print to screen and log
    (infile_pass, infile_msg) = ic.CheckInfile(args.infile)
    (config_pass, config_msg) = ic.CheckConfig(args.configfile)
    verify_msg = infile_msg + '\n' + config_msg
    pipelineHelpers.LogInputCheck(verify_msg)


    if not (infile_pass and config_pass):
        raise Exception('Invalid input.')

    pipeline_msg = '\n---------------------------------\n' \
        'Running pipeline\n' \
        '---------------------------------'
    print pipeline_msg
    log = pipelineHelpers.GetLogFile('MAIN')
    pipelineHelpers.Logging('INFO', log, pipeline_msg)

    num_procs = line_count * 2
    
    if( args.phase):    
       
        if(not args.splitbams):
            pipeline_run([pipeline.find_roi_bam], multiprocess=num_procs, verbose=1)
            pipeline_run([pipeline.complete_pipeline])
            
        else:
            
            params.SetSplitBamsPath(args.splitbams)
            ##to do
    else:
        print('user must provide phase VCFs') #user provides phased VCF
    t1 = time.time()