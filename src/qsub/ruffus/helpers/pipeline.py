import os
import sys
import re
from re import sub
from ruffus import *
from helpers import bamgineerTasks, bamgineerHelpers, pipelineHelpers
from helpers import runIDHelpers as rid
from helpers import parameters as params
from Methods import Functions as fns
import taskHelpers
import subprocess

import pysam
import pybedtools

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
#from pathos.multiprocessing import ProcessingPool
import signal
import itertools

#handling concurrency
import logging.handlers
from functools import partial
#import multiprocessing, threading, logging, sys, traceback,  StringIO, Queue
from functools import partial
from itertools import chain
from threading import Thread
import fnmatch
import inspect
import importlib


current_path = params.GetProgramPath()
pat_gain_event = params.GetPatGainCNV()
pat_loss_event = params.GetPatLossCNV()
mat_gain_event = params.GetMatGainCNV()
mat_loss_event = params.GetMatLossCNV()
cancer_type = params.GetCancerType()
vcf_path = params.GetVCF()
exon_path = params.GetExonPath()

log = pipelineHelpers.GetLogFile('Bamgineer')


#SOROUSH
(results_path, intermediate_path,
 sentinel_path,cancerDir,tmpdir,tmpbams,splittmpbams, finalbams,patfinalbams,matfinalbams) = taskHelpers.GetProjectPaths(bamgineerHelpers.name)


import vcf
import gzip
import shutil
chr_list = [12, 13,14]
haplotype_list=['pat','mat']
event_list=['gain','loss']


def runCommand(cmd):
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
    except OSError as e:
        print("Execution failed:", e)

@files(bamgineerTasks.SplitBamsTaskList)
def SplitBams(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """splits bam according to chromosome"""

    task_list = []
    log_msg = ' [SplitBams] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        samtools = '/mnt/work1/software/samtools/1.2/bin/samtools'
        python = sys.executable
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamgineerHelpers.name)
        bamgineer_mem = bamgineerHelpers.GetBamgineerMem('med')

        exons = params.GetBedFile()
        command=  "".join(["""awk '($1 ~ "chr"){print $0 >> $1".bed" }' """, exons])

        os.chdir(script_path)
        runCommand(command)

        initialize(pat_gain_event,pat_loss_event,mat_gain_event,mat_loss_event,cancer_type  ,vcf_path, exon_path)
        #chr_count = 21

        for  outp, num in zip(outputs[0], chr_list):

            #num = str(chr_count)
            script = open(
                '{0}splitbam_chr{1}.sh'.format(script_path,
                                                     num), 'w')

            script.write('#!/bin/bash\n\n')
            exons = params.GetBedFile()

            script.write('{rc} view -bh {bam} chr{num} > '
                         '{out}\n'.format(rc=samtools, bam=inputs[0][0],
                                          num=num, out=outp))
            script.close()

            process = pipelineHelpers.RunTask(
                os.path.abspath(script.name), 4, bamgineer_mem,
                sample_id,  bamgineerHelpers.name)

            task_list.append(process)
            #chr_count = chr_count + 1

        # Checks which tasks are complete, and reruns tasks that have failed

        pipelineHelpers.CheckTaskStatus(
                task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished 1')


@follows(SplitBams)
@files(bamgineerTasks.SortByNameTaskList)
def SortByName(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """sorting bam file by name"""

    task_list = []
    log_msg = ' [SortByName] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        sambamba = '/mnt/work1/software/sambamba/0.5.4/sambamba'
        python = sys.executable
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamgineerHelpers.name)
        bamgineer_mem = bamgineerHelpers.GetBamgineerMem('med')

        #chr_count = 21
        for outp, inp, num  in zip( outputs[0], inputs[0], chr_list):

                #num = str(chr_count)
                script = open(
                    '{0}sortbyname_chr{1}.sh'.format(script_path,
                                                         num), 'w')

                script.write('#!/bin/bash\n\n')


                script.write('{sb} sort -n {inp} '
                             '-o {outp} -t 4 \n'.format(sb=sambamba, inp = inp,
                                               path=current_path, outp=outp))

                script.close()

                process = pipelineHelpers.RunTask(
                    os.path.abspath(script.name), 16, bamgineer_mem,
                    sample_id,  bamgineerHelpers.name)

                task_list.append(process)
                #chr_count = chr_count + 1

        # Checks which tasks are complete, and reruns tasks that have failed
        pipelineHelpers.CheckTaskStatus(
                    task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished 2')


@follows(SortByName)
@files(bamgineerTasks.FindRoiBamTaskList)
def FindRoiBam(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """finding ROI bam for each haplotype/event/chr"""

    task_list = []
    log_msg = ' [FindRoiBam] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        sambamba = '/mnt/work1/software/sambamba/0.5.4/sambamba'
        bedtools = '/mnt/work1/software/bedtools/2.23.0/bin/bedtools'
        python = sys.executable
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamgineerHelpers.name)
        bamgineer_mem = bamgineerHelpers.GetBamgineerMem('med')


        for hap,  event in itertools.product(haplotype_list,  event_list):
            splittmpbams_hap = "/".join([splittmpbams, hap])
            tmpdir_hap = "/".join([tmpdir, hap])
            hapev = eval(hap +'_' + event +'_event')
            if(not hapev is None):#input event (specific haplotype)
              
                for inp,chr   in zip(  inputs[0], chr_list):
    
                        #name= os.path.basename(outp)
                        #het & non het
                        op = "/".join([splittmpbams_hap, 'chr'+str(chr)+'_'+event+'_het_roi.bam'])
                        op2 = "/".join([splittmpbams_hap, 'chr'+str(chr) + '_'+event +  "_non_het_roi.bam"])
                        
                        bedfn = "/".join([tmpdir_hap, event+'_het_chr' + str(chr) + '.bed'])
                        bedfn2= "/".join([tmpdir_hap,event + '_non_het_chr'+ str(chr) + '.bed'])
    
                        script = open(
                            '{0}find_roi_chr{1}_{2}_{3}_h.sh'.format(script_path,
                                                                 chr, event, hap), 'w')

                        script.write('#!/bin/bash\n\n')
                        script.write('{bt} pairtobed -abam {inp} '
                                     '-b {bf} -type either > {outp} \n'.format(bt=bedtools, inp = inp,
                                                       bf=bedfn, outp=op))
                        script.write('{sb} sort {outp} '
                                      '{outpsorted} \n'.format(sb=sambamba,bt=bedtools, outp=op, outpsorted= sub('.bam$','.sorted',op) ))
                        
                        script.write('rm {outp} \n'.format( outp=op))
                        script.close()
    
                        process = pipelineHelpers.RunTask(
                            os.path.abspath(script.name), 4, bamgineer_mem,
                            sample_id,  bamgineerHelpers.name)
                        
                        
                        script2 = open(
                            '{0}find_roi_chr{1}_{2}_{3}_nh.sh'.format(script_path,
                                                                 chr, event, hap), 'w')
                        
                        
                        script2.write('#!/bin/bash\n\n')
                        script2.write('{bt} pairtobed -abam {inp} '
                                     '-b {bf} -type either > {outp} \n'.format(bt=bedtools, inp = inp, bf=bedfn2, outp=op2))
                        script2.write('{sb} sort {outp} '
                                      '{outpsorted2} \n'.format(sb=sambamba,bt=bedtools, outp=op2, outpsorted2= sub('.bam$','.sorted',op2) ))
                        
                        script2.write('rm {outp} \n'.format( outp=op2))
                        script2.close()
    
                        process2 = pipelineHelpers.RunTask(
                            os.path.abspath(script2.name), 16, bamgineer_mem,
                            sample_id,  bamgineerHelpers.name)
                        
                        task_list.append(process)
                        task_list.append(process2)
                        
        pipelineHelpers.CheckTaskStatus(
                    task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished 3')


@follows(FindRoiBam)
@files(bamgineerTasks.MutateReadsTaskList)
def MutateReads(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """mutating reads and finding reads not matching hg19 at germline SNP locations"""
    task_list = []
    log_msg = ' [FindRoiBam] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):

        python = sys.executable
        script_path = pipelineHelpers.GetScriptPath(
                sample_id, bamgineerHelpers.name)
        bamgineer_mem = bamgineerHelpers.GetBamgineerMem('med')


        for  inp  in inputs[0]:
                    
            path,fn= os.path.split("/".join([os.getcwd(),inp]))
            chr, event, ht, bn= re.split(r'_',fn)
            
            
            for hap  in haplotype_list:
                splittmpbams_hap = "/".join([splittmpbams, hap])
                tmpdir_hap = "/".join([tmpdir, hap])
                hapev = eval(hap +'_' + event +'_event')
                if(not hapev is None):#input event (specific haplotype)
                   
                    command = " ".join(["samtools view", inp, "| less | head -1 | wc -l" ])
                    nline= subprocess.check_output(command, shell = True) 
                    
                    
                    if ((int(nline) == 0)):
                        os.remove(inp)
                        print('removing 1 '+ inp )
                        
                    
                    else:
                        
                        op = "/".join([splittmpbams_hap, str(chr) +'_'+event+'_het_alt_roi.bam'])
                        op2 = "/".join([splittmpbams_hap, str(chr) + '_'+event +  "_non_het_roi.bam"])
                        
                        
                        bedfn = "/".join([tmpdir_hap,  event + '_het_snp_' + str(chr) +  '.bed'])
                        bedfn2= "/".join([tmpdir_hap, event + '_non_het_'+ str(chr) + '.bed'])
                        
                        script = open(
                            '{0}mutate_{1}_{2}_{3}.sh'.format(script_path,
                                                                 chr, event, hap), 'w')
                        script.write('#!/bin/bash\n\n')
                        script.write('python {path}/mutateReads.py {bf} {inp} '
                                     ' {outp} {tmp}\n'.format(fns=fns,inp = inp,  outp=op,bf=bedfn, path=current_path, tmp=tmpdir))
                        
                        script.close()
                        
                        
                        process = pipelineHelpers.RunTask( 
                            os.path.abspath(script.name), 16, bamgineer_mem,
                            sample_id, bamgineerHelpers.name)
                        task_list.append(process)
                  
        
        # Checks which tasks are complete, and reruns tasks that have failed
        pipelineHelpers.CheckTaskStatus(
                    task_list, output_sentinel, log, log_msg)
       
         
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished 4')


@follows(MutateReads)
@files(bamgineerTasks.RemoveOverlappingReadsTaskList)
def RemoveOverlapping(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """mutating reads and finding reads not matching hg19 at germline SNP locations"""
    task_list = []
    log_msg = ' [FindRoiBam] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):

        python = sys.executable
        script_path = pipelineHelpers.GetScriptPath(
                sample_id, bamgineerHelpers.name)
        bamgineer_mem = bamgineerHelpers.GetBamgineerMem('med')


        for  in1, in2  in zip(inputs[0], inputs[1]):
            
            path,fn= os.path.split("/".join([os.getcwd(),in1]))
            chr, event, ht, bn= re.split(r'_',fn)
            for hap  in haplotype_list:
                splittmpbams_hap = "/".join([splittmpbams, hap])
                overlapsbam  = "/".join([splittmpbams_hap, chr + event + "_overlaps.bam"])
                hapev = eval(hap +'_' + event +'_event')
                NONHET = "/".join([splittmpbams_hap,  (chr + '_'+ event).upper()  +'_NONHET.bam'])
                HET_ALT =  "/".join([splittmpbams_hap, (chr + '_'+event).upper() +'_HET_ALT.bam'])
                if(not hapev is None):
                    tmpdir_hap = "/".join([tmpdir, hap])
                    bedfn = "/".join([tmpdir_hap, event+ '_het_'+ chr+'.bed'])
    
                    script = open(
                                '{0}removeoverlap_{1}_{2}_{3}.sh'.format(script_path,
                                                                     chr, event, hap), 'w')
                    script.write('#!/bin/bash\n\n')
                    script.write('python {path}/removeOverlaps.py {nhbam} {altbam} '
                                 ' {nhuniq} {huniq} {tmp}\n'.format(nhbam=in2, altbam=in1 , nhuniq=NONHET, huniq=HET_ALT,path=current_path, tmp=splittmpbams_hap))
                    
                    script.close()
                    
                    
                    process = pipelineHelpers.RunTask( 
                        os.path.abspath(script.name), 8, bamgineer_mem,
                        sample_id, bamgineerHelpers.name)
                    task_list.append(process)
                            
                 
        # Checks which tasks are complete, and reruns tasks that have failed
        pipelineHelpers.CheckTaskStatus(
                    task_list, output_sentinel, log, log_msg)
       
         
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished 5')


@follows(RemoveOverlapping)
@files(bamgineerTasks.ImplementGainLossTaskList)
def ImplementGainLoss (inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """mutating reads and finding reads not matching hg19 at germline SNP locations"""
    task_list = []
    log_msg = ' [FindRoiBam] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):

        python = sys.executable
        script_path = pipelineHelpers.GetScriptPath(
                sample_id, bamgineerHelpers.name)
        bamgineer_mem = bamgineerHelpers.GetBamgineerMem('med')
        
        for  in1, in2  in zip(inputs[0], inputs[1]):
            
            path,fn= os.path.split("/".join([os.getcwd(),in1]))
            
            chr, event, ht, bm = re.split(r'_',fn)
            
           
            for hap  in haplotype_list:
                splittmpbams_hap = "/".join([splittmpbams, hap])
                overlapsbam  = "/".join([splittmpbams_hap, chr + event + "_overlaps.bam"])
                hapev = eval(hap +'_' + event.lower() +'_event')
                 
                if(not hapev is None):
                    tmpdir_hap = "/".join([tmpdir, hap])
                    final_hap = "/".join([finalbams, hap])
                   
                    if(event == "GAIN"):
                        script = open(
                                    '{0}gain_{1}_{2}_{3}.sh'.format(script_path,
                                                                         chr.lower(), event.lower(), hap), 'w')
                        script.write('#!/bin/bash\n\n')
                        script.write('python {path}/implementGain.py "{ch}" "{ev}"'
                                     ' {splttmp} {fbams}\n'.format(ch=chr, ev= event, fbams = final_hap,
                                                                               path=current_path, splttmp=splittmpbams_hap))
                        
                        script.close()
                        
                        
                        process = pipelineHelpers.RunTask( 
                            os.path.abspath(script.name), 16, bamgineer_mem,
                            sample_id, bamgineerHelpers.name)
                        task_list.append(process)
                     
                    elif(event == "LOSS"):
                        script = open(
                                    '{0}loss_{1}_{2}_{3}.sh'.format(script_path,
                                                                         chr.lower(), event.lower(), hap), 'w')
                        script.write('#!/bin/bash\n\n')
                        script.write('python {path}/implementLoss.py "{ch}" "{ev}"'
                                     ' {splttmp} {fbams} {splitpath}\n'.format(ch=chr.lower(), ev= event, fbams = final_hap,
                                                                               path=current_path, splttmp=splittmpbams_hap, splitpath=intermediate_path))
                        
                        script.close()
                        
                        
                        process = pipelineHelpers.RunTask( 
                            os.path.abspath(script.name), 8, bamgineer_mem,
                            sample_id, bamgineerHelpers.name)
                        task_list.append(process)
                 
        # Checks which tasks are complete, and reruns tasks that have failed
        pipelineHelpers.CheckTaskStatus(
                    task_list, output_sentinel, log, log_msg)
       
         
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished 6')


@follows(ImplementGainLoss)
@files(bamgineerTasks.SortMergeTaskList)
def SortMerge (inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """merging and sorting the final bams"""
    task_list = []
    log_msg = ' [FindRoiBam] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):

        python = sys.executable
        sambamba = '/mnt/work1/software/sambamba/0.5.4/sambamba'
        script_path = pipelineHelpers.GetScriptPath(
                sample_id, bamgineerHelpers.name)
        bamgineer_mem = bamgineerHelpers.GetBamgineerMem('high')
        mergedbamname = params.GetOutputFileName()
       
        script = open('{0}mergesort.sh'.format(script_path), 'w')
        script.write('#!/bin/bash\n\n')
        script.write('python {path}/mergesort.py '
                                     ' {mergedfinal} {finalbamdir}\n'.format(path=current_path,  mergedfinal=mergedbamname, finalbamdir=finalbams,))

        script.close()
                        
                        
        process = pipelineHelpers.RunTask( os.path.abspath(script.name), 4, bamgineer_mem,
                            sample_id, bamgineerHelpers.name)
        task_list.append(process)
                 
        # Checks which tasks are complete, and reruns tasks that have failed
        pipelineHelpers.CheckTaskStatus(
                    task_list, output_sentinel, log, log_msg)
       
         
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished 7')
 
@follows(SortMerge)   
def CompletePipeline():
    """cleaninng up"""
    #shutil.rmtree(tmpbams)
    print "finished"

###########################################################################################################################
def removeEmptyBams(bamdir):
   for file in os.listdir(bamdir):  
     if file.endswith(".bam"):
        command = " ".join(["samtools view", "/".join([bamdir, file]), "| less | head -1 | wc -l" ])
        nline =  nline= subprocess.check_output(command, shell = True)  
        #print('command: '+ command+ '   ' +str(nline))
        if (os.path.isfile( "/".join([bamdir, file])) and (int(nline) == 0)):
                os.remove("/".join([bamdir, file]))
                print(' removing ' + "/".join([bamdir, file]))

def createDirectory(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def phaseVCF(vcf, phasevcffn):
    print ("___ phasing vcf file ___")
    if(not vcf.endswith('.vcf.gz')):
        gzipFile(vcf)
        vcf = vcf+'.gz'

    path, vcffn = os.path.split(vcf)
    phasevcffn = sub('.vcf.gz$', '_phased', vcffn)
    command = " ".join(["/mnt/work1/software/java/8/jdk1.8.0_45/bin/java -Xmx6g -jar /mnt/work1/users/pughlab/projects/Benchmarking/Beagle/beagle.09Nov15.d2a.jar", "gt="+vcf, "out="+ "/".join([tmpdir, phasevcffn])])
    print(command)
    runCommand (command)
    return phasevcffn

def thinVCF(invcf, outvcf):
   command = " ".join(["vcftools --vcf", invcf, "--thin 25 --out", outvcf,  "--recode"])
   print("thin VCF called with command: "+command )
   runCommand(command)

def convertvcftobed(vcf, bed):

    vcfh = open(vcf, 'r')
    bedh = open(bed, 'w')

    for line in vcfh:
        c = line.strip('\n').split("\t")

        if (not line.startswith('#') and len(c) >= 5 and (len(c[3])+len(c[4]) == 2)):
           start = int(c[1]) - 1
           bedh.write(c[0]+'\t'+str(start) +'\t'+ str(c[1])+'\t' + str(c[3]) + '\t' + str(c[4]) + '\n') #chr start stop ref alt

    bedh.close()

def gzipFile(filename):
    with open(filename, 'rb') as f_in, gzip.open(filename+'.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)



def getVCFPaternalMaternal(phasedvcf, pat, mat):
    out_pat = open(pat, 'w')
    out_mat = open(mat, 'w')


    if(phasedvcf.endswith('.vcf.gz')):
        vcfh = gzip.GzipFile(phasedvcf, 'rb')


        for line in vcfh:
            c = line.strip('\n').split("\t")
            if (len(c) == 10 ):
                if(c[9] == '0|1:1'):
                    out_pat.write(line)
                    continue
                elif(c[9] == '1|0:1'):
                     out_mat.write(line)
                     continue

            elif(line.startswith('#')):
                out_pat.write(line)
                out_mat.write(line)
                continue

    out_pat.close()
    out_mat.close()

def subtractBeds(bedfn1, bedfn2, diffn):
    cwd = os.path.dirname(__file__)
    bed1fullpath = os.path.realpath(bedfn1.format(cwd))
    bed2fullpath = os.path.realpath(bedfn2.format(cwd))

    bed1 = pybedtools.example_bedtool(bed1fullpath)
    bed2 = pybedtools.example_bedtool(bed2fullpath)

    print (bed2fullpath +"\n" +bed2fullpath +"\n"+ diffn)
    f = open(diffn, 'w')
    print >> f, bed1.subtract(bed2, A = True)
    f.close()


def intersectBed(bed1fn, bed2fn, intersectfile, wa=False):
    cwd = os.path.dirname(__file__)

    bed1fncompletepath = os.path.realpath(bed1fn.format(cwd))
    bed2fncompletepath = os.path.realpath(bed2fn.format(cwd))

    bed1 = pybedtools.example_bedtool(bed1fncompletepath)
    bed2 = pybedtools.example_bedtool(bed2fncompletepath)

    f = open(intersectfile, 'w')
    if(wa == False):
        print("intersect bed called with wa=False")
        print >> f, bed1.intersect(bed2,u=True)
    elif(wa==True):
        print("intersect bed called with wa=True")
        print >> f, bed1.intersect(bed2,u=True,wa=True)

    f.close()
    return intersectfile

def extendbed(inbed, outbed):
    inbedh = open(inbed, 'r')
    outbedh = open(outbed, 'w')

    for line in inbedh:
        c = line.strip('\n').split("\t")

        if (len(c) >= 3 ):
           start = int(c[1]) - 25
           stop = int(c[1]) + 25
           outbedh.write(c[0]+'\t'+str(start) +'\t'+ str(stop) +'\n')

    outbedh.close()

def splitBed(bedfn, event):
    path, filename = os.path.split(bedfn)

    command=  "".join(["""awk '($1 ~ "chr"){print $0 >> """ ,'"{}"'.format(event), """$1".bed"}' """, bedfn])
    os.chdir(path)
    runCommand(command)

def initialize(patgaincnv, patlosscnv, matgaincnv, matlosscnv, cancerType, vcfpath, exons):
    global inbam, tmpdir, tmpbams, splitbams, finalbams, vcffiltered, vcffilteredtobed, vcffilteredtobed_padded,  gaincnv, losscnv, ref, splittmpbams ,cancerDir, matfinalbams, patfinalbams
    global matvcffiltered, matvcffilteredtobed, matvcffilteredtobed_padded,  patvcffiltered, patvcffilteredtobed, patvcffilteredtobed_padded, tmpdirpat, tmpdirmat  #Phase params

    tmpdirpat = "/".join([tmpdir, "pat"])
    tmpdirmat = "/".join([tmpdir, "mat"])
    createDirectory(patfinalbams)
    createDirectory(matfinalbams)
    createDirectory(tmpdirpat)
    createDirectory(tmpdirmat)
    createDirectory(splittmpbams)
    createDirectory(finalbams)

    try:
        #logger.debug(' --- Initialization called --- ')
        vpath, vcf = os.path.split(vcfpath)
        phasedvcf = "/".join([tmpdir, sub('.vcf$', '_phased.vcf.gz', vcf)])
        patvcf = "/".join([tmpdir,"pat" ,"het.vcf"])
        matvcf = "/".join([tmpdir, "mat", "het.vcf"])
        patvcffiltered = "/".join([tmpdir, "pat","het_filtered"])
        matvcffiltered = "/".join([tmpdir, "mat","het_filtered"])
        patvcffilteredtobed = "/".join([tmpdir, "pat","het_filtered.bed"])
        matvcffilteredtobed = "/".join([tmpdir, "mat","het_filtered.bed"])
        patvcffilteredtobed_padded = "/".join([tmpdir, "pat","mat_het_padded.bed"]) #NOT USED YET
        matvcffilteredtobed_padded = "/".join([tmpdir, "mat","het_padded.bed"]) #NOT USED YET


        phaseVCF(vcfpath, phasedvcf)
        getVCFPaternalMaternal(phasedvcf, patvcf, matvcf)
        thinVCF(patvcf, patvcffiltered)
        thinVCF(matvcf, matvcffiltered)
        convertvcftobed(patvcffiltered+".recode.vcf", patvcffilteredtobed)
        convertvcftobed(matvcffiltered+".recode.vcf", matvcffilteredtobed)

        for hap in haplotype_list:
            splittmpbams_hap = "/".join([splittmpbams, hap])
            tmpdir_hap = "/".join([tmpdir, hap])
            #splitbams_hap = "/".join([splitbams, hap])
            finalbams_hap = "/".join([finalbams, hap])

            createDirectory(splittmpbams_hap)
            createDirectory(tmpdir_hap)

            for  event in event_list:
                roibed = "/".join([tmpdir,  hap ,event + "_roi.bed"])
                exonsinroibed = "/".join([tmpdir, hap ,  event + "_exons_in_roi.bed"])
                exonsinhetbed = "/".join([tmpdir, hap , event + "_exons_in_het.bed"])
                nonhetbed = "/".join([tmpdir, hap , event + "_non_het.bed"])
                hetbed = "/".join([tmpdir, hap , event + "_het.bed"])
                hetsnpbed = "/".join([tmpdir,  hap , event + "_het_snp.bed"])
                cnvsovelrapexons = "/".join([tmpdir, hap ,  event + "_FINAL_CNVS.bed"])

                if  (eval(str((hap +event + 'cnv')))):

                    hapcnvevent = eval((hap +event + 'cnv'))
                    intersectBed( exons, hapcnvevent, exonsinroibed, wa=True)
                    intersectBed(  hapcnvevent, exons, cnvsovelrapexons, wa=True)
                    intersectBed( exonsinroibed, patvcffilteredtobed, exonsinhetbed, wa=True)
                    intersectBed(patvcffilteredtobed, exonsinroibed, hetsnpbed, wa=True)
                    extendbed(patvcffilteredtobed, patvcffilteredtobed_padded) # add padding around snp location +/-50
                    subtractBeds(exonsinroibed, patvcffilteredtobed_padded , nonhetbed) # so that
                    intersectBed( patvcffilteredtobed, exonsinroibed, hetbed, wa=True)

                    splitBed(hetsnpbed,  event+'_het_snp_')
                    splitBed(hetbed, event+'_het_')

                    nonhetbed = "/".join([tmpdir, hap, event + "_non_het.bed"])
                    splitBed(nonhetbed, event + '_non_het_')

    except:

        raise
    return




