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
import subprocess

global bases
bases = ('A','T','C','G')

log = pipelineHelpers.GetLogFile('Bamgineer')
import utils
import vcf
import gzip
import shutil
chr_list = range(1,22)
event_list=['gain','loss']

sentinel_path, results_path,haplotype_path,cancer_dir_path,tmpbams_path,finalbams_path = taskHelpers.GetProjectNamePathRunID()


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

            if (locals()[event + 'cnv']):
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

            script.write('#!/bin/bash\n')
            script.write('#\n')
            script.write('#$ -cwd \n')
            script.write('module load samtools \n')
            script.write(command +'\n')
            script.write('samtools view -bh {bam} chr{num} > '
                         '{out}\n'.format(bam=inputs[0][0],
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
        python = sys.executable
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('med')
        for outp, inp, num  in zip( outputs[0], inputs[0], chr_list):

                script = open(
                    '{0}sortbyname_chr{1}.sh'.format(script_path,
                                                         num), 'w')

                script.write('#!/bin/bash\n')
                script.write('#\n')
                script.write('#$ -cwd \n')
                script.write('module load sambamba \n')
                script.write('sambamba sort -n {inp} '
                             '-o {outp} -t 4 \n'.format(inp = inp, outp=outp))
                script.close()
                process = pipelineHelpers.RunTask(
                    os.path.abspath(script.name), 4, bamgineer_mem,
                    sample_id,  bamhelp.name)
                task_list.append(process)
        pipelineHelpers.CheckTaskStatus(
                    task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished SplitBams')


@follows(sort_by_name)
@files(bamgineerTasks.find_roi_bam_task_list)
def find_roi_bam(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """finding ROI bam for each haplotype/event/chr"""
    task_list = []
    log_msg = ' [FindRoiBam] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        python = sys.executable
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('med')
        sentinel_path, results_path,haplotype_path,cancer_dir_path,tmpbams_path,finalbams_path = taskHelpers.GetProjectNamePathRunID()
        
        for inp, op in izip(inputs[0],outputs[0]):
            opsorted=sub('.bam$',".sorted.bam", op)
            chr=os.path.basename(op).strip().split(".")[0]
            event=os.path.basename(op).strip().split(".")[1]
            exonsinroibed = "/".join([haplotype_path,   event + "_exons_in_roi_"+ str(chr) +'.bed'])

            if (os.path.isfile(exonsinroibed)):
                script = open(
                    '{0}find_roi_{1}_{2}.sh'.format(script_path,
                                                         chr, event), 'w')
                script.write('#!/bin/bash\n\n')
                script.write('#\n')
                script.write('#$ -cwd \n')
                script.write('module load bedtools \n')
                script.write('module load sambamba \n')

                script.write('sort -u {exonbed} -o {exonbed} \n'.format(exonbed=exonsinroibed))
                script.write('bedtools pairtobed -abam {inp} '
                             '-b {bf} -type either > {outp} \n'.format(inp = inp,
                                               bf=exonsinroibed, outp=op))
                script.write('sambamba sort {outp} -o '
                              '{outpsorted} \n'.format(outp=op, outpsorted= opsorted))
                script.write('rm {outp} \n'.format( outp=op))
                script.close()
                process = pipelineHelpers.RunTask(
                    os.path.abspath(script.name), 1, bamgineer_mem,
                    sample_id,  bamhelp.name)

                task_list.append(process)
        pipelineHelpers.CheckTaskStatus(
                    task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished FindROI')
    

@follows(find_roi_bam)
@files(bamgineerTasks.repair_task_list)
def repair_gain( inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """implementing cnv module and finding reads not matching hg19 at germline SNP locations"""
    task_list = []
    log_msg = ' [re-pairing reads] ' + '[' + sample_id + '] '
    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')

    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        python = sys.executable
        current_path = params.GetProgramPath()
        script_path = pipelineHelpers.GetScriptPath(
                sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('high')

        print('**** ' + str(inputs[0]))

        for inp in inputs[0]:

            chr= os.path.basename(inp).strip().split(".")[0]

            script = open('{0}re-pair_{1}_{2}.sh'.format(script_path, chr, "gain"), 'w')
            script.write('#!/bin/bash\n\n')
            script.write('module load samtools/1.2 \n')
            script.write('module load sambamba \n')
            script.write('python {path}/re-pair.py {inbam} \n'.format(inbam=inp, path=current_path ))
            script.close()
            process = pipelineHelpers.RunTask(
                os.path.abspath(script.name),4, bamgineer_mem,
                sample_id, bamhelp.name)
            task_list.append(process)

        pipelineHelpers.CheckTaskStatus(
                        task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished Re-pairing')


@follows(repair_gain)
@files(bamgineerTasks.mutate_gain_task_list)
def mutate_gain(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """mutating reads according to haplotype at germline SNP locations"""
    task_list = []
    log_msg = ' [implement_cnv] ' + '[' + sample_id + '] '

    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
        python = sys.executable

        current_path = params.GetProgramPath()
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('med')
        sentinel_path, results_path,haplotype_path,cancer_dir_path,tmpbams_path,finalbams_path = taskHelpers.GetProjectNamePathRunID()

        for inp in inputs[0]:

            chr= os.path.basename(inp).strip().split(".")[0]


            bedfn= "/".join([haplotype_path, 'gain_het_snp_' + chr + '.bed'])
            diffn =   "/".join([tmpbams_path,"diff.bam"])
            nonhet= "/".join([tmpbams_path, 'diff_only1_' +  os.path.basename(inp)])
            hetfn=sub('.gain.roi.repaired.sorted.bam$','.gain.mutated.het.bam', inp)
            hetfnsorted = sub('.gain.roi.repaired.sorted.bam$','.gain.mutated.het.sorted.bam', inp)
            mergedsortfn = sub('.gain.roi.repaired.sorted.bam$','.gain.mutated.merged.sorted.bam', inp)
            mergedrenamedfn = sub('.gain.roi.repaired.sorted.bam$','.gain.renamed.mutated.merged.sorted.bam', inp)

            script = open('{0}mutate_{1}_{2}.sh'.format(script_path, chr, "gain"), 'w')
            script.write('#!/bin/bash\n')
            script.write('#')
            script.write('#$ -cwd \n')
            script.write('module load samtools/1.2 \n')
            script.write('module load sambamba \n')
            script.write('module load bamUtil \n')

            script.write('sort -u {bf} -o {bf}\n\n'.format(bf=bedfn))
            script.write('python {path}/mutate.py {repairedbam} {bf} {happath}\n\n'.format(repairedbam=inp, bf=bedfn ,path=current_path , happath=haplotype_path))
            script.write('sambamba sort {het} -o {hetsort}\n\n'.format(het=hetfn, hetsort=hetfnsorted))
            script.write('bam diff --in1 {repairedbam} --in2 {hetsort} --out {dif}\n\n'.format(repairedbam=inp, hetsort=hetfnsorted ,dif=diffn ))
            script.write('sambamba merge {merged} {hetonly} {nonhetonly}\n\n'.format(merged=mergedsortfn,hetonly=hetfnsorted, nonhetonly= nonhet))
            script.write('rm {het} {nonhetonly}  \n\n'.format(het=hetfn,nonhetonly= nonhet))
            script.write('python {path}/rename-reads.py {inp2} {outp}\n\n'.format(inp2= mergedsortfn, outp=mergedrenamedfn, path=current_path))

            script.close()
            process = pipelineHelpers.RunTask(
                os.path.abspath(script.name),4, bamgineer_mem,
                sample_id, bamhelp.name)
            task_list.append(process)

            pipelineHelpers.CheckTaskStatus(
                            task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished Mutating')


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
            chrevent=os.path.basename(inp).strip().split("_")[0]
            chr = re.split('(\d+)',chrevent)[1]

            original_bam = sub('.repaired.mutated.merged.sorted.bam', '.sorted.bam', inp)
            sentinel_path, results_path,haplotype_path,cancer_dir_path,tmpbams_path,finalbams_path = taskHelpers.GetProjectNamePathRunID()
            GAIN_FINAL = "/".join([finalbams_path,  'CHR'+str(chr).upper() +'_GAIN.bam'])

            script = open('{0}sample_{1}_{2}.sh'.format(script_path, 'chr'+str(chr), "gain"), 'w')
            script.write('#!/bin/bash\n')
            script.write('#\n')
            script.write('#$ -cwd \n')
            script.write('module load samtools/1.2 \n')
            script.write('python {path}/subsample_gain.py {inbam} {origroi} {fg} \n'.format(path=current_path,inbam=inp, origroi=original_bam , fg=GAIN_FINAL))

            script.close()
            process = pipelineHelpers.RunTask(
                os.path.abspath(script.name),4, bamgineer_mem,
                sample_id, bamhelp.name)
            task_list.append(process)

            pipelineHelpers.CheckTaskStatus(
                            task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg+ 'Finished Sampling Gain Event')


@follows(find_roi_bam)
@files(bamgineerTasks.mutate_loss_task_list)
def mutate_loss(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """mutating reads according to haplotype at germline SNP locations"""
    task_list = []
    log_msg = ' [implement_cnv] ' + '[' + sample_id + '] '

    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
        python = sys.executable

        current_path = params.GetProgramPath()
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('med')
        sentinel_path, results_path,haplotype_path,cancer_dir_path,tmpbams_path,finalbams_path = taskHelpers.GetProjectNamePathRunID()

        for inp in inputs[0]:

            chr= os.path.basename(inp).strip().split(".")[0]

            bedfn= "/".join([haplotype_path, 'loss_het_snp_' + chr + '.bed'])

            if (os.path.isfile(inp) and os.path.isfile(bedfn)):

                diffn =   "/".join([tmpbams_path,"diff.bam"])
                nonhet= "/".join([tmpbams_path, 'diff_only1_' +  os.path.basename(inp)])
                hetfn=sub('.roi.sorted.bam$',".mutated.het.bam", inp)
                hetfnsorted = sub('.roi.sorted.bam$',".mutated.het.sorted.bam", inp)
                mergedsortfn = sub('.roi.sorted.bam$',".mutated.merged.sorted.bam", inp)

                script = open('{0}mutate_{1}_{2}.sh'.format(script_path, chr, "loss"), 'w')
                script.write('#!/bin/bash\n')
                script.write('#')
                script.write('#$ -cwd \n')
                script.write('module load samtools/1.2 \n')
                script.write('module load sambamba \n')
                script.write('module load bamUtil \n')

                script.write('sort -u {bf} -o {bf}\n\n'.format(bf=bedfn))
                script.write('python {path}/mutate.py {inp1} {bf} {happath}\n\n'.format(inp1=inp, bf=bedfn ,path=current_path , happath=haplotype_path))
                script.write('sambamba sort {het} -o {hetsort}\n\n'.format(het=hetfn, hetsort=hetfnsorted))
                script.write('bam diff --in1 {repairedbam} --in2 {hetsort} --out {dif}\n\n'.format(repairedbam=inp, hetsort=hetfnsorted ,dif=diffn ))
                script.write('sambamba merge {merged} {hetonly} {nonhetonly}\n\n'.format(merged=mergedsortfn,hetonly=hetfnsorted, nonhetonly= nonhet))
                script.write('rm {het} {nonhetonly}  \n\n'.format(het=hetfn,nonhetonly= nonhet))

                script.close()
                process = pipelineHelpers.RunTask(
                    os.path.abspath(script.name),4, bamgineer_mem,
                    sample_id, bamhelp.name)
                task_list.append(process)

                pipelineHelpers.CheckTaskStatus(
                                task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg + 'Finished Mutating')


@follows(mutate_loss)
@files(bamgineerTasks.subsample_loss_task_list)
def subsample_loss(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """adjusting sample rate for Bam files"""
    task_list = []
    log_msg = ' [subsample loss events] ' + '[' + sample_id + '] '
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):
        pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
        python = sys.executable
        current_path = params.GetProgramPath()
        script_path = pipelineHelpers.GetScriptPath(
            sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('med')

        for inp in inputs[0]:
            chrevent=os.path.basename(inp).strip().split("_")[0]
            chr = re.split('(\d+)',chrevent)[1]
            original_bam = sub('.mutated.merged.sorted.bam', '.sorted.bam', inp)
            sentinel_path, results_path,haplotype_path,cancer_dir_path,tmpbams_path,finalbams_path = taskHelpers.GetProjectNamePathRunID()
            LOSS_FINAL = "/".join([finalbams_path,  'CHR'+str(chr).upper() +'_LOSS.bam'])

            script = open('{0}sample_{1}_{2}.sh'.format(script_path, 'chr'+str(chr), "loss"), 'w')
            script.write('#!/bin/bash\n')
            script.write('#\n')
            script.write('#$ -cwd \n')
            script.write('module load samtools/1.2 \n')
            script.write('python {path}/subsample_loss.py {inbam} {fl} \n'.format(path=current_path,inbam=inp, fl=LOSS_FINAL))

            script.close()
            process = pipelineHelpers.RunTask(
                os.path.abspath(script.name),4, bamgineer_mem,
                sample_id, bamhelp.name)
            task_list.append(process)

            pipelineHelpers.CheckTaskStatus(
                            task_list, output_sentinel, log, log_msg)
    pipelineHelpers.Logging('INFO', log, log_msg+ 'Finished Sampling Loss Event')



@follows(subsample_loss)
@follows(subsample_gain)
@files(bamgineerTasks.complete_pipeline_task_list)
def complete_pipeline(inputs, output_sentinel, outputs, sample_id, prev_sentinel):
    """merge, sort, clean up """
    task_list = []
    log_msg = ' [Final merge] ' + '[' + sample_id + '] '

    pipelineHelpers.Logging('INFO', log, log_msg + 'Starting')
    if pipelineHelpers.CheckSentinel(prev_sentinel, log, log_msg):

        python = sys.executable
        current_path = params.GetProgramPath()
        script_path = pipelineHelpers.GetScriptPath(
                sample_id, bamhelp.name)
        bamgineer_mem = bamhelp.GetBamgineerMem('high')
        mergedbamname = params.GetOutputFileName()

        script = open('{0}mergesort.sh'.format(script_path), 'w')
        script.write('#!/bin/bash\n')
        script.write('#\n')
        script.write('#$ -cwd \n')
        script.write('module load sambamba \n')

        script.write('python {path}/mergesort.py '
                                     ' {mergedfinal} {finalbamdir}\n'.format(path=current_path,  mergedfinal=mergedbamname, finalbamdir=finalbams_path))

        script.close()
        process = pipelineHelpers.RunTask( os.path.abspath(script.name), 4, bamgineer_mem,
                            sample_id, bamhelp.name)
        task_list.append(process)
        pipelineHelpers.CheckTaskStatus(
                    task_list, output_sentinel, log, log_msg)


    pipelineHelpers.Logging('INFO', log, log_msg + 'COMPLETE!')
