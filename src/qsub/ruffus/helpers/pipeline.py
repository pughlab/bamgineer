import pysam
import ntpath
from helpers import parameters as params
from helpers import bamgineerHelpers as bamhelp
from helpers import bamgineerTasks, pipelineHelpers
import time
import os
import random
from helpers import pipelineHelpers
from ruffus import *
import logging, sys
from re import sub
import taskHelpers
from itertools import izip
import itertools
import re

global bases
bases = ('A','T','C','G')

log = pipelineHelpers.GetLogFile('Bamgineer')

java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path,bamutil_path=params.GetSoftwarePath()

import utils
import vcf
import gzip
import shutil
chr_list = [4, 5]
event_list=['gain','loss']

def initPool(queue, level, terminating_):
    """
    This causes the logging module to be initialized with the necessary info
    in pool threads to work correctly.
    """
    logging.getLogger('').setLevel(level)
    global terminating
    terminating = terminating_


def initialize():
    try:
        utils.createDirectory(results_path)
        utils.createDirectory(cancer_dir_path)
        utils.createDirectory(haplotype_path)
        utils.createDirectory(tmpbams_path)
        utils.createDirectory(finalbams_path)  
        
        event_list=['gain','loss']
        gaincnv = params.GetGainCNV()
        losscnv = params.GetLossCNV()
        
        pipelineHelpers.Logging("INFO", log, " --- Initializing input files  --- ")
        vcf_path = bamhelp.GetVCF()
        exons_path = bamhelp.GetExons()
        reference_path = bamhelp.GetRef()
        vpath, vcf = os.path.split(vcf_path)
        phasedvcf = "/".join([results_path, sub('.vcf$', '_phased.vcf.gz', vcf)])
        vcftobed =  "/".join([results_path, sub('.vcf$', '.bed', vcf)])
        
        hap1vcf = "/".join([results_path,"hap1_het.vcf"])
        hap2vcf = "/".join([results_path, "hap2_het.vcf"])
        hap1vcffiltered = "/".join([results_path, "hap1_het_filtered"])
        hap2vcffiltered = "/".join([results_path, "hap2_het_filtered"])
        hap1vcffilteredtobed = "/".join([results_path, "hap1_het_filtered.bed"])
        hap2vcffilteredtobed = "/".join([results_path, "hap2_het_filtered.bed"])
        phased_bed =  "/".join([results_path, "PHASED.BED"])
        
        
        utils.phaseVCF(vcf_path, phasedvcf)
        utils.getVCFHaplotypes(phasedvcf, hap1vcf, hap2vcf)
        utils.thinVCF(hap1vcf, hap1vcffiltered)
        utils.thinVCF(hap2vcf, hap2vcffiltered)
        utils.convertvcftobed(hap1vcffiltered+".recode.vcf", hap1vcffilteredtobed)
        utils.convertvcftobed(hap2vcffiltered+".recode.vcf", hap2vcffilteredtobed)
       
        cmd1 = """sed -i 's/$/\thap1/' """+ hap1vcffilteredtobed
        cmd2 = """sed -i 's/$/\thap2/' """+ hap2vcffilteredtobed
        cmd3 = "cat " + hap1vcffilteredtobed + " " + hap2vcffilteredtobed + " > " + 'tmp.bed'
        cmd4 = "sort -V -k1,1 -k2,2 tmp.bed > " + phased_bed  
            
        utils.runCommand(cmd1)
        utils.runCommand(cmd2)
        utils.runCommand(cmd3)
        utils.runCommand(cmd4)
        os.remove('tmp.bed')  
        
        for  event in event_list: 
            roibed = "/".join([haplotype_path,  event + "_roi.bed"])
            exonsinroibed = "/".join([haplotype_path,   event + "_exons_in_roi.bed"])
            nonhetbed = "/".join([haplotype_path, event + "_non_het.bed"])
            hetbed = "/".join([haplotype_path, event + "_het.bed"])
            hetsnpbed = "/".join([haplotype_path,  event + "_het_snp.bed"])
            
            utils.intersectBed( exons_path, locals()[event + 'cnv'], exonsinroibed, wa=True)
            utils.intersectBed(phased_bed, exonsinroibed, hetsnpbed, wa=True)
            utils.splitBed(exonsinroibed, event+'_exons_in_roi_')
            utils.splitBed(hetsnpbed, event+'_het_snp_')

    except:  
        pipelineHelpers.Logging("INFO", log, "Initialization error !")
        raise
    
    pipelineHelpers.Logging("ERROR",log, "--- initialization complete ---")    
    return 


@files(bamgineerTasks.split_bam_task_list)
def split_bams(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """splits bam according to chromosome"""
    task_list = []
    log_msg = ' [SplitBams] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        samtools = samtools_path
        python = sys.executable
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('med')

        exons = params.GetExonPath()
        command=  "".join(["""awk '($1 ~ "chr"){print $0 >> $1".bed" }' """, exons])

        initialize()
        
        os.chdir(script_path)
        #utils.runCommand(command)
        for  outp, num in zip(outputs[0], chr_list):
            script = open(
                '{0}splitbam_chr{1}.sh'.format(script_path,
                                                     num), 'w')

            script.write('#!/bin/bash\n\n')
            script.write(command +'\n')
            script.write('{rc} view -bh {bam} chr{num} > '
                         '{out}\n'.format(rc=samtools, bam=inputs[0][0],
                                          num=num, out=outp))
            script.close()

            process = pipelineHelpers.RunTask(
                os.path.abspath(script.name), 4, bamgineer_mem,
                sample_id,  bamhelp.name)
            task_list.append(process)
           
        pipelineHelpers.CheckTaskStatus(
                task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished task1')


@follows(split_bams)
@files(bamgineerTasks.sort_by_name_task_list)
def sort_by_name(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """sorting bam file by name"""

    task_list = []
    log_msg = ' [SortByName] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        sambamba = sambamba_path
        python = sys.executable
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('med')
        for outp, inp, num  in zip( outputs[0], inputs[0], chr_list):

                script = open(
                    '{0}sortbyname_chr{1}.sh'.format(script_path,
                                                         num), 'w')

                script.write('#!/bin/bash\n\n')
                script.write('{sb} sort -n {inp} '
                             '-o {outp} -t 4 \n'.format(sb=sambamba, inp = inp, outp=outp))
                script.close()

                process = pipelineHelpers.RunTask(
                    os.path.abspath(script.name), 4, bamgineer_mem,
                    sample_id,  bamhelp.name)

                task_list.append(process)
        pipelineHelpers.CheckTaskStatus(
                    task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished task2')


@follows(sort_by_name)
@files(bamgineerTasks.find_roi_bam_task_list)
def find_roi_bam(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """finding ROI bam for each haplotype/event/chr"""
    task_list = []
    log_msg = ' [FindRoiBam] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        sambamba = sambamba_path
        bedtools = bedtools_path
        python = sys.executable
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('med')

        for event,chr in itertools.product(event_list,chr_list):
            roi,sortbyname,sortbyCoord, hetsnp = init_file_names('chr'+str(chr), event, tmpbams_path, haplotype_path)
            exonsinroibed = "/".join([haplotype_path,   event + "_exons_in_roi_"+ 'chr'+str(chr) +'.bed'])
            roisort = sub('.bam$', '.sorted.bam', roi)
            if(os.path.isfile(exonsinroibed)):
    
                script = open(
                    '{0}find_roi_chr{1}_{2}.sh'.format(script_path,
                                                         chr, event), 'w')
                script.write('#!/bin/bash\n\n')
                script.write('sort -u {exonbed} -o {exonbed} \n'.format(exonbed=exonsinroibed))
                script.write('{bt} pairtobed -abam {inp} '
                             '-b {bf} -type either > {outp} \n'.format(bt=bedtools, inp = sortbyname,
                                               bf=exonsinroibed, outp=roi))
                script.write('{sb} sort {outp} -o '
                              '{outpsorted} \n'.format(sb=sambamba, outp=roi, outpsorted= roisort))  
                script.write('rm {outp} \n'.format( outp=roi))
                script.close()
                process = pipelineHelpers.RunTask(
                    os.path.abspath(script.name), 1, bamgineer_mem,
                    sample_id,  bamhelp.name)
                task_list.append(process)                 
        pipelineHelpers.CheckTaskStatus(
                    task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished 3')



@follows(find_roi_bam)
@files(bamgineerTasks.repair_task_list)
def repair( inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """implementing cnv module and finding reads not matching hg19 at germline SNP locations"""
    task_list = []
    log_msg = ' [re-pairing reads] ' + '[' + sample_id + '] '
    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    
    
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        sambamba = sambamba_path
        python = sys.executable
        current_path = params.GetProgramPath()
        script_path = pipelineHelpers.GetScriptPath(
                sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('med')
       
        for chr in chr_list:
        
            bamfn,sortbyname,sortbyCoord, bedfn = init_file_names('chr'+str(chr), "gain", tmpbams_path, haplotype_path)
            bamsortfn = sub('.bam$', '.sorted.bam', bamfn)
            
            if(os.path.isfile(bamsortfn)):
            
                script = open('{0}re-pair_{1}_{2}.sh'.format(script_path, 'chr'+str(chr), "gain"), 'w')
                script.write('#!/bin/bash\n\n')
                script.write('python {path}/re-pair.py {inbam} {sam} {bamba} \n'.format(inbam=bamsortfn, path=current_path, sam = samtools_path, bamba=sambamba_path ))        
                script.close()
                process = pipelineHelpers.RunTask( 
                    os.path.abspath(script.name),4, bamgineer_mem,
                    sample_id, bamhelp.name)
                task_list.append(process)
            
        pipelineHelpers.CheckTaskStatus(
                        task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished Re-pairing')


@follows(repair)
@files(bamgineerTasks.mutate_gain_task_list)
def mutate_gain(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """mutating reads according to haplotype at germline SNP locations"""
    task_list = []
    log_msg = ' [implement_cnv] ' + '[' + sample_id + '] '
    
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
        python = sys.executable
        sambamba = sambamba_path
        current_path = params.GetProgramPath()
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('med')
        for inp in inputs[0]:
            chrevent=os.path.basename(inp).strip().split("_")[0]
            chr = re.split('(\d+)',chrevent)[1]
           
            sentinel_path, results_path,haplotype_path,cancer_dir_path,tmpbams_path,finalbams_path = taskHelpers.GetProjectNamePathRunID()
            bedfn= "/".join([haplotype_path, 'gain_het_snp_chr' + chr + '.bed'])
            diffn =   "/".join([tmpbams_path,"diff.bam"])
            nonhet= "/".join([tmpbams_path, 'diff_only1_' +  os.path.basename(inp)])
            hetfn=sub('.sorted.bam$',".mutated_het.bam", inp)
            hetfnsorted = sub('.sorted.bam$',".mutated_het.sorted.bam", inp)
            mergedsortfn = sub('.sorted.bam$',".mutated_merged.sorted.bam", inp)
            mergedrenamedfn = sub('.sorted.bam$',".mutated_merged_renamed.sorted.bam", inp)
            
            script = open('{0}mutate_{1}_{2}.sh'.format(script_path, 'chr'+str(chr), "gain"), 'w')
            script.write('#!/bin/bash\n')
            script.write('#')
            script.write('#$ -cwd')
            script.write('sort -u {bf} -o {bf}\n\n'.format(bf=bedfn))
            script.write('python {path}/mutate.py {repairedbam} {bf} {happath}\n\n'.format(repairedbam=inp, bf=bedfn ,path=current_path , happath=haplotype_path))        
            script.write('{samba} sort {het} -o {hetsort}\n\n'.format(samba= sambamba_path,het=hetfn, hetsort=hetfnsorted))
            
            script.write('{bamdiff} diff --in1 {repairedbam} --in2 {hetsort} --out {dif}\n\n'.format(bamdiff= bamutil_path, repairedbam=inp, hetsort=hetfnsorted ,dif=diffn,path=current_path ))  
            script.write('{samba} merge {merged} {hetonly} {nonhetonly}\n\n'.format(samba= sambamba_path,merged=mergedsortfn,hetonly=hetfnsorted, nonhetonly= nonhet))
            script.write('rm {het} {nonhetonly}  \n\n'.format(het=hetfn,nonhetonly= nonhet))
            script.write('python {path}/rename-reads.py {inp2} {outp}\n\n'.format(inp2= mergedsortfn, outp=mergedrenamedfn, path=current_path))
            
            script.close()
            process = pipelineHelpers.RunTask( 
                os.path.abspath(script.name),4, bamgineer_mem,
                sample_id, bamhelp.name)
            task_list.append(process)
                
            pipelineHelpers.CheckTaskStatus(
                            task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished Re-pairing')
        

@follows(mutate_gain)
@files(bamgineerTasks.subsample_gain_task_list)
def subsample_gain(inputs, output_sentinel, outputs, sample_id, prev_sentinel):    
    """adjusting sample rate for Bam files"""
    task_list = []
    log_msg = ' [implement_cnv] ' + '[' + sample_id + '] '
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
        python = sys.executable
        current_path = params.GetProgramPath()
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('med')
        
        for inp in inputs[0]:
            print()
        

@follows(find_roi_bam)
@files(bamgineerTasks.mutate_gain_task_list)
def mutate_loss():
    """mutating reads according to haplotype at germline SNP locations"""
    task_list = []
    log_msg = ' [implement_cnv] ' + '[' + sample_id + '] '
    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        print()
        


def init_file_names(chr, event,tmpbams_path, haplotypedir):
    
    flist=[]
    results_path = params.GetConfigReader().get('RESULTS', 'results_path')
    split_path = "/".join([results_path, "splitbams"])
    roibam = "/".join([tmpbams_path ,chr + event +"_roi.bam"])
    sortbyname =  "/".join([split_path,  chr + '.byname.bam'])
    sortbyCoord = "/".join([split_path,  chr + '.bam'])
    hetsnp   = "/".join([haplotypedir, event+'_het_snp_' + chr + '.bed'])
    flist.extend([roibam,sortbyname,sortbyCoord,hetsnp])
    return flist

