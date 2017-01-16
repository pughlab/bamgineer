import argparse
import os
from ruffus import *
from helpers import runIDHelpers as rid
from helpers import parameters as params
from helpers import inputChecker as ic
from helpers import pipelineHelpers
from helpers import taskHelpers

__version__ = 'DepthOfCoverage v0.1.0'

# main scope
if __name__ == '__main__':

    # command line arguments
    parser = argparse.ArgumentParser(description=__version__)
    
    ####      SOROUSH        #####
    ##############################
   
    
    
    
    parser.add_argument('-outbam', dest='outBamFile', required=True,
                        help='.bam file name for output')
    parser.add_argument('-vcf', dest='vcfFile', required=True,
                        help='input vcf file ')
    #Phasing 
    parser.add_argument('-cnv_amp_pat', dest='cnvAmpPatFile', required=False,
                        help='CNV amplification .bed file name')
    parser.add_argument('-cnv_del_pat', dest='cnvDelPatFile', required=False,
                        help='CNV deletion .bed file name')
    parser.add_argument('-cnv_amp_mat', dest='cnvAmpMatFile', required=False,
                        help='CNV amplification .bed file name')
    parser.add_argument('-cnv_del_mat', dest='cnvDelMatFile', required=False,
                        help='CNV deletion .bed file name')   
    parser.add_argument('-inbam', dest='inbamFile', required=False,
                        help='sam/bam file from which to obtain reads')
    parser.add_argument('-exons', dest='exonsFile', required=True,
                        help='Exon .bed file name')
    parser.add_argument('-cancertype', dest='cancerType', required=False,
                        help='acronyms for cancer type')
    parser.add_argument('-phase',dest= 'phase', action="store_true")
    ##############################
    parser.add_argument(
        '--project-name', action='store', required=True, dest='project_name',
        help='name of the project')
    parser.add_argument(
        '-i', '--infile', action='store', dest='infile', required=True,
        help='2 column (tab delimited) input file: ID  tumour_path')
    parser.add_argument(
        '--project-path', action='store', required=True, dest='project_path',
        help='specify results directory')
    parser.add_argument(
        '-c', '--configFile', action='store', required=True, dest='configfile',
        help='/path/to/config_file.cfg')
    parser.add_argument(
        '-r', '--run_id', action='store', type=int, dest='custom_run_id',
        help='manually enter run_id (used to re-launch a run)', default=False)
    #parser.add_argument('--version', action='version', version=__version__)
        
    args = parser.parse_args()
    current_path = os.path.dirname(os.path.realpath(__file__))
    
    # write command line arguments to parameters module
    params.InitConfigReader(args.configfile)
    params.SetProjectName(args.project_name)
    params.SetProjectPath(args.project_path)
    params.SetReferenceVersion('hg19')
    params.SetReference('hg19')
    params.SetProgramPath(current_path)
    params.SetInfile(args.infile)

    #SOROUSH
    params.SetPatGainCNV(args.cnvAmpPatFile)
    params.SetPatLossCNV(args.cnvDelPatFile)
    params.SetMatGainCNV(args.cnvAmpMatFile)
    params.SetMatLossCNV(args.cnvDelMatFile)
    params.SetCancerType(args.cancerType)
    params.SetVCF(args.vcfFile)
    params.SetExonPath(args.exonsFile)
    params.SetOutputFileName(args.outBamFile)


    # set the run id used for this run
    if args.custom_run_id == False:
        # increment and set runID
        rid.IncRunID(params.GetProjectName(), current_path)
    else:
        # use runid specified on the command line -r
        rid.SetRunID(args.custom_run_id, current_path)

        # check that the user specified run direcotry exists
        if not os.path.exists(taskHelpers.GetProjectNamePathRunID()):
            print (str(taskHelpers.GetProjectNamePathRunID()))
            raise Exception('The project folder for this run does not exist')

    
    # ----------------------------------------------------------------------------#
    #  Launch Depth Of Coverage Pipeline
    # ----------------------------------------------------------------------------#
    from helpers import pipeline

  
    def LaunchPipeline():
        """initiates pipeline found in pipeline.py"""
        print 'DepthOfCoverage COMPLETE'

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

    # terminate run if the inputChecker returns any failures
    if not (infile_pass and config_pass):
        raise Exception('Invalid input.')

    pipeline_msg = '\n---------------------------------\n' \
        'Running pipeline\n' \
        '---------------------------------'
    print pipeline_msg
    log = pipelineHelpers.GetLogFile('MAIN')
    pipelineHelpers.Logging('INFO', log, pipeline_msg)

    # LaunchPipeline()
    #num_procs = line_count * 3 number of chrom
    num_procs = line_count * 16
    #pipeline_run([LaunchPipeline], multiprocess=num_procs, verbose=1)
    
    #pipeline_run([pipeline.SplitBams, pipeline.SortByName, pipeline.FindRoiBam], multiprocess=num_procs, verbose=1)
    pipeline_run([pipeline.CompletePipeline], multiprocess=num_procs, verbose=1)
   

