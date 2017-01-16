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


@subdivide(["1.input_file"],
            regex(r"(.+).bam"),      # match file prefix
           [r"\1.chr1",
            r"\1.chr2",
            r"\1.chr3"])

def SplitBam():
    





##global split_finished
##split_finished = False
#@files(depthOfCoverageTasks.SplitBamsTaskList)
#def SplitBams(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
#    """splits bam according to chromosome"""
#
#    task_list = []
#    log_msg = ' [SplitBams] ' + '[' + sample_id + '] '
#    #fin = False #SOROUSH
#
#    
#    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
#    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
#        samtools = '/mnt/work1/software/samtools/1.2/bin/samtools'
#        python = sys.executable
#        script_path = pipelineHelpers.GetScriptPath(
#            sample_id, depthOfCoverageHelpers.name)
#        depthOfCoverage_mem = depthOfCoverageHelpers.GetDepthOfCoverageMem('med')
#     
#        exons = params.GetBedFile()
#        command=  "".join(["""awk '($1 ~ "chr"){print $0 >> $1".bed" }' """, exons]) 
#           
#        os.chdir(script_path)
#        runCommand(command)
#        
#        chr_count = 1
#        
#        for  outp in outputs[0]:
#             
#            num = str(chr_count)         
#            script = open(
#                '{0}splitbam_chr{1}.sh'.format(script_path, 
#                                                     num), 'w')
#            script.write('#!/bin/bash\n\n')
#            exons = params.GetBedFile()
#            
#            script.write('{rc} view -bh {bam} chr{num} > '
#                         '{out}\n'.format(rc=samtools, bam=inputs[0][0],
#                                          num=num, path=current_path, out=outp))
#            script.write('{rc} index {out}\n'.format(rc=samtools, path=current_path, out=outp))
#            script.close()
#
#            process = pipelineHelpers.RunTask(
#                os.path.abspath(script.name), 4, depthOfCoverage_mem,
#                sample_id,  depthOfCoverageHelpers.name)
#           
#            task_list.append(process)
#            chr_count = chr_count + 1
#    
#            # Checks which tasks are complete, and reruns tasks that have failed
#            
#        pipelineHelpers.CheckTaskStatus(
#                task_list, output_sentinel, log, log_msg)
#    #if(fin):
#        pipelineHelpers.Logging('INFO', log, log_msg + 'Finished XXX')
#    
#    #with io.FileIO("split.log", "w") as file:
#    #    file.write("split done!")
#    #split_finished = True
#
##depth_finished = False
##@active_if(split_finished)
#@follows(SplitBams)
##@check_if_uptodate(SplitBams)
#@files(depthOfCoverageTasks.DepthOfCovTaskList)
#def DepthOfCov(inputs, output_sentinel, outputs, sample_id, prev_sentinel): #SOROUSH added task_done
#
#    """depth of coverage for each chromosome"""
#
#    task_list = []
#    log_msg = ' [DepthOfCov] ' + '[' + sample_id + '] '
#
#    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
#    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
#        #gatk = params.GetGATK()
#        gatk ='java -Xmx8g -jar /mnt/work1/software/gatk/3.0-0/GenomeAnalysisTK.jar'
#        reference = params.GetReference()
#        exons = params.GetBedFile()
#       
#        
#        python = sys.executable
#        script_path = pipelineHelpers.GetScriptPath(
#            sample_id, depthOfCoverageHelpers.name)
#        depthOfCoverage_mem = depthOfCoverageHelpers.GetDepthOfCoverageMem('med')
#        
#       
#        #for chr in chr_list:
#        chr_count = 1
#        for outp, inp in zip( outputs[0], inputs[0]):
#         
#            #outp = sub('.bam$', '', outp)
#            
#            num = str(chr_count)
#            chr = str("chr"+num)
#            
#
#            script = open(
#                '{0}depthcov_chr{1}.sh'.format(script_path,
#                                                     num), 'w')
#            # Write GATK command to the script depending on the chromosome
#            
#            script.write('#!/bin/bash\n')
#            script.write('module load  gatk/3.0-0\n')
#            script.write('{gatk} -R {ref} -T DepthOfCoverage -o {out} '
#                         '-I {inp} -L {bedfn} --filter_mismatching_base_and_quals ' #SOROUSH remove the flag later
#                         '\n'.format(gatk =gatk, inp=inp,
#                                           path=current_path, out=outp, bedfn= "/".join([script_path, chr+'.bed']), ref=reference))
#                        
#            script.close()
#
#            process = pipelineHelpers.RunTask(
#                os.path.abspath(script.name), 4, depthOfCoverage_mem,
#                sample_id,  depthOfCoverageHelpers.name)
#            # Add task to tasklist
#            task_list.append(process)
#            chr_count = chr_count + 1
#
#        # Checks which tasks are complete, and reruns tasks that have failed
#        pipelineHelpers.CheckTaskStatus(
#            task_list, output_sentinel, log, log_msg)
#        
#    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished')
#    
#    #with io.FileIO("depth.log", "w") as file:
#    #    file.write("depth done!")
#    depth_finished = True
#
#
##@active_if(depth_finished)
#@follows(DepthOfCov)
##@check_if_uptodate(DepthOfCov)
#@files(depthOfCoverageTasks.MergeCovFilesTaskList)
#
#def MergeCovFiles(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
#    
#    """merge the results and create final files"""
#
#    task_list = []
#    log_msg = ' [MergeCoverageFiles] ' + '[' + sample_id + '] '
#
#    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
#    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
#
#        python = sys.executable
#        script_path = pipelineHelpers.GetScriptPath(
#            sample_id, depthOfCoverageHelpers.name)
#        depthOfCoverage_mem = depthOfCoverageHelpers.GetDepthOfCoverageMem('med')
#
#    file_types =['cov', 'cov.sample_cumulative_coverage_counts',
#                 'cov.sample_cumulative_coverage_proportions','cov.sample_interval_statistics',
#                 'cov.sample_statistics','cov.sample_interval_summary','cov.sample_summary']
#        
#    for type in file_types:
#        
#        for inp in inputs:
#           
#            script = open(
#                    '{0}merge_{1}.sh'.format(script_path, type), 'w')
#            script.write('#!/bin/bash\n')   
#
#            fname,ext = str(inp[0]).split(".", 1)
#            
#            if (ext == type):
#                type = sub('cov', '', type,1)
#                
#                script.write('head -1 {inpt} > {head}\n'.format(inpt = str(inp[0]), head = "/".join([script_path, 'header.txt'+type])))
#                
#                
#                for count in range (1 , 23):
#                    str2 = "'{}'".format('1d')
#                    script.write('sed -i {str} {infile}\n'.format(str=str2, infile=str(inp[count - 1])))
#                 
#                script.write('cat {head} '.format(head = "/".join([script_path, 'header.txt'+type])))    
#                for count in range (1 , 23):
#                    
#                    script.write(str(inp[count - 1]) + ' ')
#                
#            else:
#                continue
#            script.write('> ' + str(outputs[0][0]) + str(type)) #outputs[0][0]
#                    
#            script.close()
#         
#            process = pipelineHelpers.RunTask(
#            os.path.abspath(script.name), 4, depthOfCoverage_mem,
#                sample_id,  depthOfCoverageHelpers.name)
#            
#            task_list.append(process)
#    
#    pipelineHelpers.CheckTaskStatus(
#            task_list, output_sentinel, log, log_msg)
#    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished')