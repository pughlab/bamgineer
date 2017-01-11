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

global bases
bases = ('A','T','C','G')
chromosome_list = ['chr1','chr2','chr3','chr4', 'chr5', 'chr6','chr7','chr8','chr9', 'chr10','chr11','chr12', 'chr13','chr14','chr15','chr16', 'chr17', 'chr18','chr19','chr20','chr21', 'chr22']
purity = [ '0.2', '0.4', '0.6','0.8','1.0']
event_list=['gain','loss']
cancer_list = ['brca','crc','gbm']

chromosome_event =   [ 'chr1_gain','chr2_gain','chr3_gain','chr4_gain', 'chr5_gain', 'chr6_gain','chr7_gain','chr8_gain','chr9_gain', 'chr10_gain','chr11_gain', 
                       'chr12_gain','chr13_gain','chr14_gain','chr15_gain','chr16_gain', 'chr17_gain', 'chr18_gain','chr19_gain','chr20_gain', 'chr21_gain', 'chr22_gain',
                       'chr1_loss','chr2_loss','chr3_loss','chr4_loss', 'chr5_loss', 'chr6_loss','chr7_loss','chr8_loss','chr9_loss', 'chr10_loss','chr11_loss',
                      'chr12_loss','chr13_loss','chr14_loss','chr15_loss','chr16_loss', 'chr17_loss', 'chr18_loss','chr19_loss','chr20_loss','chr21_loss', 'chr22_loss']  


def runCommand(cmd):
    try:
        thread = Thread(group=None, target=lambda:os.system(cmd))
        thread.run()
        # Later
        if not thread.is_alive():
            return 0
        else:
            return 1
    except OSError as e:
        
        logger.debug("Execution failed %s", e)
        return -1
        
def createDirectory(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def gzipFile(filename):
    with open(filename, 'rb') as f_in, gzip.open(filename+'.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

#filter vcf so that all SNPs are at least 20bp apart:    
def thinVCF(invcf, outvcf):
   command = " ".join(["vcftools --vcf", invcf, "--thin 20 --out", outvcf,  "--recode"])
   print("thin VCF called with command: "+command )
   runCommand(command)

def extractPairedReadfromROI(inbamfn, bedfn, outbamfn, flag = "either"):
    
    command = " ".join(["bedtools pairtobed -abam", inbamfn,"-b",bedfn, "-type", flag,">", outbamfn])
    runCommand(command)
    
def extractBAMfromROI_All(inbamfn, bedfn, outbamfn):
    command = " ".join(["samtools view -b -L", bedfn, inbamfn, ">",outbamfn])
    runCommand(command)

def dedupBam(inbamfn, outbamfn):
    command = " ".join(["/cluster/tools/software/samtools/0.1.18/samtools rmdup", inbamfn, outbamfn])
    runCommand(command)

def bamDiff(bamfn1, bamfn2,path):
    command = " ".join(["bam diff", "--in1", bamfn1, "--in2", bamfn2, "--out" ,"/".join([path,"diff.bam"])]) # ("roi.bam - vcf.bam"; seperate reads that do not overlap SNP regions from the ones that do)
    runCommand(command ) #this is in fact from BamUtil package

#phase unphased VCF into Hap1 and Hap2 phased alleles using BEAGLE
def phaseVCF(vcf, phasevcffn):
    print (" ___ phasing vcf file ___ "  )
    if(not vcf.endswith('.vcf.gz')):
        gzipFile(vcf)
        vcf = vcf+'.gz'
        
    path, vcffn = os.path.split("/".join([os.getcwd(),vcf]))
    phasevcffn = sub('.vcf.gz$', '_phased', vcffn)
    command = " ".join(["/mnt/work1/software/java/8/jdk1.8.0_45/bin/java -Xmx4g -jar /mnt/work1/users/pughlab/projects/Benchmarking/Beagle/beagle.09Nov15.d2a.jar", "gt="+vcf, "out="+ "/".join([RESULTS, phasevcffn])])
    print(command)
    runCommand (command)
    return phasevcffn

def getVCFHaplotypes(phasedvcf, hap1, hap2):
    out_hap1 = open(hap1, 'w')
    out_mat = open(hap2, 'w')
    
    
    if(phasedvcf.endswith('.vcf.gz')):
        vcfh = gzip.GzipFile(phasedvcf, 'rb')
        
                
        for line in vcfh:
            c = line.strip('\n').split("\t")
            if (len(c) == 10 ):
                if(c[9] == '0|1:1'):
                    out_hap1.write(line)
                    continue
                elif(c[9] == '1|0:1'):
                     out_mat.write(line)
                     continue
                
            elif(line.startswith('#')):
                out_hap1.write(line)
                out_mat.write(line)
                continue
        
    out_hap1.close()
    out_mat.close()
   
def getCNVPaternalMaternal(cnv, hap1, mat):
    out_hap1 = open(hap1, 'w')
    out_mat = open(mat, 'w')
     
    cnvh = open(cnv, 'r')
    
    for line in cnvh:
        c = line.strip('\n').split("\t")
        if (len(c) == 5 ):
            if(c[3] == 'hap1ernal'):
                out_hap1.write(line)
                continue
            elif(c[3] == 'maternal'):
                out_mat.write(line)
                continue
    out_hap1.close()
    out_mat.close()
 
def convertvcftobed(vcf, bed):
    
    vcfh = open(vcf, 'r')
    bedh = open(bed, 'w')
    
    for line in vcfh:
        c = line.strip('\n').split("\t")
        
        if (not line.startswith('#') and len(c) >= 5 and (len(c[3])+len(c[4]) == 2)):
           start = int(c[1]) - 1
           bedh.write(c[0]+'\t'+str(start) +'\t'+ str(c[1])+'\t' + str(c[3]) + '\t' + str(c[4]) + '\n') #chr start stop ref alt
        
    bedh.close()
    

#report exon bedfiles that are within the defined CNV region
def findExonsinCNVregion(cnvpath, exonpath, intersectfile, wa=False):
    cwd = os.path.dirname(__file__)
    cnvcompletepath = os.path.realpath(cnvpath.format(cwd))
    exoncompletepath = os.path.realpath(exonpath.format(cwd))
      
    cnvfile = pybedtools.example_bedtool(cnvcompletepath)
    exonfile = pybedtools.example_bedtool(exoncompletepath)
    
    f = open(intersectfile, 'w')
    if(wa == False):
        print >> f, exonfile.intersect(cnvfile,u=True)
    elif(wa==True):
        print >> f, exonfile.intersect(cnvfile,u=True,wa=True)
            
    f.close()
    return intersectfile

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

def MutateReads(bedfn, bamfn, outfn):
    
    logger.debug("___ mutating reads and finding reads not matching hg19 at germline SNP locations ___")
    samfile = pysam.Samfile(bamfn, "rb" )
    
    path, filename = os.path.split(bedfn)
    outbam = pysam.Samfile(outfn, 'wb', template=samfile) 
    
    refbamfn = sub('_het_alt_roi.bam$',"_REFREADS.bam", outfn)
    altbamfn = sub('_het_alt_roi.bam$',"_ALTREADS.bam", outfn)
    hapbamfn = sub('_het_alt_roi.bam$',"_HAP.bam", outfn)
       
    
    refbam = pysam.Samfile(refbamfn, 'wb', template=samfile) 
    altbam = pysam.Samfile( altbamfn, 'wb', template=samfile) 
    hapbam = pysam.Samfile(hapbamfn, 'wb', template=samfile) 
    
    sortedsamfn = sub('.bam$', '.sorted', bamfn)

    pysam.sort(bamfn, sortedsamfn)
    pysam.index(sortedsamfn+".bam")
    
    alignmentfile = pysam.AlignmentFile(sortedsamfn+".bam", "rb" )
   
    bedfile = open(bedfn, 'r')
    
    covpath = "/".join([haplotypedir, "written_coverage_het.txt"])
    covfile = open(covpath, 'w')
    
    snpratiopath = "/".join([haplotypedir, "het_snp_ratio.txt"])
    snpaltratiofile = open(snpratiopath,'w')
    
    writtenreads = []
    readscoveringsnpregion = []
    for bedline in bedfile:
       
        c = bedline.strip().split()
        if (len(c) == 6 ):
            chr2 = c[0]
            chr = c[0].strip("chr")
            start = int(c[1])
            end   = int(c[2])
            refbase = str(c[3])
            altbase = str(c[4])
            haplotype = str(c[5])
        else :
            
            continue
        
        readmappings = alignmentfile.fetch(chr2, start, end)
      
        total_num_reads_at_this_locus = 0
        num_ref_reads_at_this_locus = 0
        num_non_ref_reads_at_this_locus = 0
        num_reads_written = 0
          
        for shortread in readmappings:
            readscoveringsnpregion.append(shortread.qname)           
            if(shortread.is_paired and shortread.is_proper_pair and not shortread.is_duplicate  
               and not shortread.is_secondary and not shortread.qname in writtenreads and shortread.mapping_quality >= 30 ):
                
                try:
                    
                    index = shortread.get_reference_positions().index(start)
                    nuecleotide = shortread.seq[index]
                    mate = alignmentfile.mate(shortread)
                    
                    if (refbase in bases and nuecleotide in bases ):
                        if(haplotype == "hap1"):
                            if(nuecleotide != refbase and nuecleotide == altbase  ):
                                altbam.write(shortread)#SOROUSH
                                altbam.write(mate)
                                outbam.write(shortread)
                                outbam.write(mate)
                                hapbam.write(shortread)#SOROUSH
                                hapbam.write(mate)
                                writtenreads.append(shortread.qname)
                                num_reads_written += 2
                                continue
                                    
                            elif(nuecleotide == refbase ):
                                
                                refbam.write(shortread)#SOROUSH
                                refbam.write(mate)
                                tmpread = shortread.query_sequence
                                                                        
                                if ((abs(shortread.tlen) > 200)):
                                                             
                                    tmpread = shortread.query_sequence
                                    basetomutate = altbase
                                
                                    for i in range(0, len(tmpread) - 1):              
                                        mutated = tmpread[:index] +  basetomutate + tmpread[index + 1:]
                                    
                                    shortread.query_sequence = mutated
                                    outbam.write(shortread)
                                    outbam.write(mate)
                                    writtenreads.append(shortread.qname)
                                    num_reads_written += 2
                                    continue       
                        elif(haplotype == "hap2"  ):
                            if(nuecleotide == refbase and nuecleotide != altbase):
                                refbam.write(shortread)
                                refbam.write(mate)
                                hapbam.write(shortread)
                                hapbam.write(mate)
                                
                                outbam.write(shortread)
                                outbam.write(mate)
                                writtenreads.append(shortread.qname)
                                num_reads_written += 2
                                continue
                                    
                            elif(nuecleotide != refbase ): 
                                tmpread = shortread.query_sequence
                                altbam.write(shortread)
                                altbam.write(mate)
                                                                      
                                if ((abs(shortread.tlen) > 200)):
                                                             
                                    tmpread = shortread.query_sequence
                                    basetomutate = refbase
                                
                                    for i in range(0, len(tmpread) - 1):              
                                        mutated = tmpread[:index] +  basetomutate + tmpread[index + 1:]
                                    
                                    shortread.query_sequence = mutated
                                    outbam.write(shortread)
                                    outbam.write(mate)
                                    writtenreads.append(shortread.qname)
                                    num_reads_written += 2
                                    continue       
                        
                except (KeyError,ValueError) as e :
                    pass  
        
        if(float(total_num_reads_at_this_locus) > 0):    
            ratio2 =  float(num_reads_written) / float(total_num_reads_at_this_locus)            
            covfile.write(str(num_reads_written) + '\t' +str(total_num_reads_at_this_locus) + '\t' + 'ratio: '+ str(ratio2) + '\n')
    
    for read in alignmentfile:
        readsaroundsnp = 0 
        if (read.is_proper_pair and read.is_paired and readsaroundsnp < num_non_ref_reads_at_this_locus and
            not read.is_secondary and not read.qname in readscoveringsnpregion and not read.qname in writtenreads):
            
            outbam.write(read)
            readsaroundsnp +=1
            writtenreads.append(read.qname)
            try:
                outbam.write(alignmentfile.mate(read))
                writtenreads.append(read.qname)
                                                                             
            except (KeyError,ValueError) as e :
                print("no mate found for "+read.qname)
                pass  
                               
    outbam.close()
    refbam.close()
    altbam.close()
    hapbam.close()
    covfile.close()
    snpaltratiofile.close()
    os.remove(bamfn)



def call_subprocess(cmd):
    try:
        proc = subprocess.Popen(cmd, shell=False)
        out, err = proc.communicate()
    except:    
        logger.exception("Exception in call_subprocess " , sys.exc_info()[0])
        return
    
 
def runBamgineer( chromosome_event ):
    
    chr,event = chromosome_event .split("_")
    if (args.phase):
         
        hetroibam = "/".join([splittmpbams ,chr + event +"_het_roi.bam"])
        hetaltroibam = "/".join([splittmpbams, chr + event + "_het_alt_roi.bam"])
        nonhetroibam = "/".join([splittmpbams, chr + event +  "_non_het_roi.bam"])
        overlapsbam  = "/".join([splittmpbams, chr + event + "_overlaps.bam"])
        
        sortbyname =  "/".join([splitbams,  chr + '.byname.bam'])
        sortbyCoord = "/".join([splitbams,  chr + '.bam'])
        
        hetaltroibamsorted = sub('.bam$','.sorted',hetaltroibam)
      
        #For deletion
        refbam = "/".join([splittmpbams, chr + event +"_REFREADS.bam"])
        altbam = "/".join([splittmpbams, chr + event +"_ALTREADS.bam"])
        refbamsorted = sub('.bam$','.sorted', refbam)
        altbamsorted = sub('.bam$','.sorted', altbam)
        
        hapbam = "/".join([splittmpbams, chr + event +"_HAP.bam"])
        hapbamsorted = sub('.bam$','.sorted', hapbam)
       
        
        HET_ALT =  "/".join([splittmpbams, chr + event +'_HET_ALT.bam'])
        nonhetroibamsorted = sub('.bam$','.sorted',nonhetroibam)
        NONHET = "/".join([splittmpbams,  chr + event +'_NONHET.bam'])
       
        comlog =  "/".join([haplotypedir, event + chr + '.LOG'])
        
        
        try:
            if not terminating.is_set():
             
                if (os.path.isfile("/".join([haplotypedir, event+'_het_' + chr + '.bed']))):
                    
                    extractPairedReadfromROI(sortbyname,"/".join([haplotypedir, event+'_het_' + chr + '.bed']), hetroibam)
                      
                    if (os.path.isfile("/".join([haplotypedir, event + '_het_snp_' + chr +  '.bed']))):
                       
                        if(os.path.isfile(hetroibam)):
                            
                            
                            MutateReads("/".join([haplotypedir,  event + '_het_snp_' + chr +  '.bed']), hetroibam, hetaltroibam)
                            pysam.sort(hetaltroibam, hetaltroibamsorted )
                            os.remove(hetaltroibam)
                            HET_ALT =  hetaltroibamsorted+'.bam'
                            
                        else:
                            logger.debug('hetroibam not existing ')
                                    
                              
            if(os.path.isfile("/".join([haplotypedir,event + '_non_het_'+ chr + '.bed']))):
                
                extractPairedReadfromROI(sortbyname, "/".join([haplotypedir,event + '_non_het_'+ chr + '.bed']), nonhetroibam)
                if (os.path.isfile("/".join([haplotypedir,  event+'_het_' + chr + '.bed'])) and os.path.isfile(nonhetroibam)):
                   pysam.sort(nonhetroibam,nonhetroibamsorted)
                   os.remove(nonhetroibam)
                   bamDiff(nonhetroibamsorted+'.bam' , hetaltroibamsorted+'.bam' ,splittmpbams)
                
                   if (os.path.isfile("/".join([splittmpbams,  'diff_only1_'+ chr+ event + '_non_het_roi.sorted.bam']))):
                        os.rename("/".join([splittmpbams, 'diff_only1_'+ chr+ event + '_non_het_roi.sorted.bam']), NONHET )
                        logger.debug('diff_only1_'+ chr+ event + '_non_het_roi.sorted.bam' + '  **** exists and NONHET is '+ NONHET )
                else:
                  
                    NONHET = "/".join([splittmpbams, chr+ event + '_non_het_roi.sorted.bam'])
                    logger.debug('diff_only1_'+ chr+ event + '_non_het_roi.sorted.bam' + '  #### does not exists and NONHET is '+ NONHET )
                   
            
            if (event == 'gain'):
                logger.debug('in gain event ')
                gain_HET_ALT_RE_PAIR =  "/".join([splittmpbams,  str(chr) + 'gain_HET_ALT_RE_PAIR.bam'])
                gain_HET_ALT_RE_PAIR_SAMPLED =  "/".join([splittmpbams, chr + event +'_HET_ALT_RE_PAIR_SAMPLED.bam'])
                gain_HET_ALT_RE_PAIR_SAMPLED_RENAMED =  "/".join([splittmpbams ,  str(chr) +'gain_HET_ALT_RE_PAIR_SAMPLED_RENAMED.bam'])
                
                
                gain_NONHET_RE_PAIR = "/".join([splittmpbams, chr + event +'_NONHET_REPAIR.bam'])
                gain_NONHET_RE_PAIR_SAMPLED = "/".join([splittmpbams, chr + event +'_NONHET_SAMPLED.bam'])
                gain_NONHET_RE_PAIR_SAMPLED_RENAMED = "/".join([splittmpbams, str(chr) + 'gain_NONHET_RE_PAIR_SAMPLED_RENAMED.bam'])
                
                gain_NONHET_FINAL = "/".join([finalbams,  str(chr).upper() +'_GAIN_NH'])
                gain_HET_FINAL = "/".join([finalbams,  str(chr).upper() +'_GAIN_H'])
                
                if(os.path.isfile(HET_ALT)):
                    
                    samapleratehet = splitAndRePair(HET_ALT,  gain_HET_ALT_RE_PAIR, chr , "het")
                    if(samapleratehet < 0.9):
                        subsample(gain_HET_ALT_RE_PAIR,gain_HET_ALT_RE_PAIR_SAMPLED, str(samapleratehet*1.1)) # we need to keep a bit more (by 15-20%)# this factor may change 
                    else:
                        gain_HET_ALT_RE_PAIR_SAMPLED = gain_HET_ALT_RE_PAIR
                        
                    renamereads(gain_HET_ALT_RE_PAIR_SAMPLED, gain_HET_ALT_RE_PAIR_SAMPLED_RENAMED )
                    pysam.sort(gain_HET_ALT_RE_PAIR_SAMPLED_RENAMED,  gain_HET_FINAL)
                    logger.debug('adjusted sampling rate for HET'+chr +' = '+ str(samapleratehet*1.1))
                       
                if(os.path.isfile(NONHET)):
                    samapleratenonhet = splitAndRePair(NONHET, gain_NONHET_RE_PAIR, chr, "nonhet")
                    if(samapleratenonhet < 1.0):
                        subsample(gain_NONHET_RE_PAIR, gain_NONHET_RE_PAIR_SAMPLED,  str(samapleratenonhet))
                    else:
                        gain_NONHET_RE_PAIR_SAMPLED = gain_NONHET_RE_PAIR
                    renamereads(gain_NONHET_RE_PAIR_SAMPLED, gain_NONHET_RE_PAIR_SAMPLED_RENAMED)
                    pysam.sort(gain_NONHET_RE_PAIR_SAMPLED_RENAMED, gain_NONHET_FINAL)
                    
                
            elif (event == 'loss'):
                inbam_deletion = "/".join([finalbams , str(chr).upper() + 'LOSS_FINAL.bam'])
                if(not os.path.isfile(NONHET) and not os.path.isfile(HET_ALT)): 
                    os.symlink(sortbyCoord, inbam_deletion)
                    logger.debug('No loss event defined for chromosome: ' + chr +' creating symlink for ')
                
                if(os.path.isfile(NONHET) and os.path.isfile(HET_ALT)):
                    loss_NONHET_SAMPLED =   "/".join([splittmpbams,  chr + event +'_NONHET_SAMPLED.bam'])
                    loss_HET_SAMPLED =   "/".join([splittmpbams,  chr + event +'_HET_SAMPLED.bam'])
                    
                    
                    logger.debug("In LOSS module NONHET is " + NONHET )
                    subsample(NONHET ,loss_NONHET_SAMPLED , str(0.5))
                    inbam_minus_loss_NONHET_ALT = "/".join([splittmpbams, chr + event+ '_MINUS_NONHET_ALT.bam'])
                    
                    bamDiff(sortbyCoord , loss_NONHET_SAMPLED, splittmpbams)
                    os.rename("/".join([splittmpbams,  'diff_only1_' + chr + '.bam']), inbam_minus_loss_NONHET_ALT)
                    logger.debug(" loss module renaming:   "+ "/".join([splittmpbams, 'diff_only1_' + chr + '.bam'])+
                                 '   '+inbam_minus_loss_NONHET_ALT)
                      
                    pysam.sort(hapbam, hapbamsorted)
                    bamDiff(inbam_minus_loss_NONHET_ALT, hapbamsorted+'.bam', splittmpbams)
                    os.rename("/".join([splittmpbams,'diff_only1_'+ chr + event+ '_MINUS_NONHET_ALT.bam']),inbam_deletion)
                       
        except (KeyboardInterrupt):
            logger.error('Exception Crtl+C pressed in the child process  in Bamgineer for chr ' + chr)
            terminating.set()
            return
        except Exception as e:   
            logger.exception("Exception in Runbamgineer %s" ,e )
            terminating.set()
            return
        return        


def splitAndRePair(inbamfn, outbamfn, chr,param=None):
    print(" calling splitAndRePair non_het version" )
    splitfn1 = '/'.join([splittmpbams,chr+'sp1_nh.bam'])
    splitfn2 = '/'.join([splittmpbams,chr+'sp2_nh.bam'])
    
    splitPairs(inbamfn, splitfn1, splitfn2 )
    inbam = pysam.Samfile(inbamfn, 'rb')
    splt1sortedfn = sub('.bam$', '.sorted', splitfn1)
    splt2sortedfn = sub('.bam$', '.sorted', splitfn2)
    
    pysam.sort(splitfn1 , splt1sortedfn )
    pysam.sort(splitfn2 , splt2sortedfn)
    
    splt1 = pysam.Samfile(splt1sortedfn + ".bam", 'rb') 
    splt2 = pysam.Samfile(splt2sortedfn + ".bam", 'rb')
    spltcount = pysam.Samfile(splt1sortedfn + ".bam", 'rb')
    outbam = pysam.Samfile(outbamfn, 'wb', template=inbam)  

    itr1 = splt1.fetch(until_eof=True)
    itr2 = splt2.fetch(until_eof=True)
    
    num_reads_to_write = 0
    for row in spltcount:
        num_reads_to_write+=1
    
    writtencount = 0
    samplerate = 0
    start = True
    for read1, read2 in  izip(itr1, itr2):
              
        try:
           
            if(read2.qname != read1.qname and start and read2.pos < read1.pos ):
                read2 = itr2.next()
                start = False
                continue
            read1next = itr1.next()
            read2next = itr2.next()
        
        except StopIteration:
            break
        
        if(read2.rnext == read2.tid and read1.rnext == read1.tid and read1.qname != read2.qname and read2.tid == read1.tid and
           read2next.rnext == read2next.tid and read1next.rnext == read1next.tid and read1next.qname != read2next.qname and read2next.tid == read1next.tid):
            
            delta =   abs(read2.pos - read1.pos) + 1
            deltanext =   abs(read2next.pos - read1next.pos) + 1
            
            if(abs(read1.pos - read2next.pos) < 12100 and read1.mapping_quality >= 20 and read2next.mapping_quality >= 20):
               
                tlen1 = -20000
                if(read1.tlen > 0 and read2next.tlen < 0):
                    tlen1 = abs(read2next.pos - read1.pos) + abs(read2next.qlen)
                elif(read1.tlen < 0 and read2next.tlen > 0):
                    tlen1 = -abs(read2next.pos - read1.pos) - abs(read2next.qlen)
                
                if(tlen1 != -20000):
                    read1.tlen = tlen1
                    read2next.tlen = -tlen1
                  
                    read1.pnext = read2next.pos
                    read2next.pnext = read1.pos
                    read2next.qname = read1.qname
                    outbam.write(read1)
                    outbam.write(read2next)
                    writtencount = writtencount + 1
        
            if(abs(read1next.pos - read2.pos) < 12100 and read2.mapping_quality >= 20 and read1next.mapping_quality >= 20):
               
                tlen1 = -20000
                if(read1next.tlen > 0 and read2.tlen < 0):
                    tlen1 = abs(read2.pos - read1next.pos) + abs(read2.qlen)
                elif(read1next.tlen < 0 and read2.tlen > 0):
                    tlen1 = -abs(read2.pos - read1next.pos) -abs(read2.qlen)
               
                if(tlen1 != -20000):    
                    read1next.tlen = tlen1
                    read2.tlen = -tlen1
                    
                    read2.pnext = read1next.pos
                    read1next.pnext = read2.pos
                    read2.qname = read1next.qname
                    
                    outbam.write(read1next)
                    outbam.write(read2)
                    writtencount = writtencount + 1
   
    inbam.close()         
    splt1.close()
    splt2.close()
    outbam.close()
    
    if(num_reads_to_write > 0):
        percentkept = float(writtencount)/float(num_reads_to_write)
        print(' % of kept reads for '+param+ ' : '+str(percentkept))
        return(0.5/percentkept)
 
def removeReadsOverlappingHetRegion(inbamfn, bedfn,outbamfn,path):
    print "___ removing reads overlapping heterozygous region ___"
    inbamsorted =  sub('.bam$','.sorted',inbamfn)
    pysam.sort(inbamfn, inbamsorted)
    pysam.index(inbamsorted+'.bam')
    
    alignmentfile = pysam.AlignmentFile(inbamsorted+'.bam', "rb" )
    outbam = pysam.Samfile(outbamfn, 'wb', template=alignmentfile )
    
    bedfile = open(bedfn, 'r')
    
    for bedline in bedfile:
        c = bedline.strip().split()
        
        if (len(c) == 3 ):
            chr2 = c[0]
            chr = c[0].strip("chr")
            start = int(c[1])
            end   = int(c[2])
        else :
            continue
        
        try:
            readmappings = alignmentfile.fetch(chr2, start, end)
        except  ValueError as e:
            print("problem fetching the read ")
        
        
        for shortread in readmappings:
            try:
                outbam.write(shortread)
            except ValueError as e:
                print ("problem removing read :" + shortread.qname)
    outbamsorted =  sub('.bam$','.sorted',outbamfn)            
    pysam.sort(outbamfn, outbamsorted)
    print(' bamdiff ' + str(inbamsorted+'.bam') + ' ' + str(outbamsorted +'.bam'))
    bamDiff(inbamsorted+'.bam', outbamsorted +'.bam', path )
    outbam.close()           

def getMeanSTD(inbam):
    """ awk '{ if ($9 > 0) { N+=1; S+=$9; S2+=$9*$9 }} END { M=S/N; print "n="N", mean="M", stdev="sqrt ((S2-M*M*N)/(N-1))}' """
    command = " ".join([awk, inbam])
    runCommand(command )


def splitPairs(inbamfn,pair1fn, pair2fn):
    command1 = " ".join(["samtools view -u -h -f 0x0043", inbamfn, ">", pair1fn])
    command2 = " ".join(["samtools view -u -h -f 0x0083", inbamfn, ">", pair2fn])
    runCommand(command1)
    runCommand(command2)

def renamereads(inbamfn, outbamfn):
    
    inbam = pysam.Samfile(inbamfn, 'rb')
    outbam = pysam.Samfile(outbamfn, 'wb', template=inbam)

    paired = {}

    n = 0
    p = 0
    u = 0
    w = 0
    m = 0
    
    for read in inbam.fetch(until_eof=True):
        n += 1
        if read.is_paired:
            p += 1
            if read.qname in paired:
                uuid = paired[read.qname]
                del paired[read.qname]
                read.qname = uuid
                outbam.write(read)
                w += 1
                m += 1
            else:
                newname = str(uuid4())
                paired[read.qname] = newname
                read.qname = newname
                outbam.write(read)
                w += 1
        else:
            u += 1
            read.qname = str(uuid4())
            outbam.write(read)
            w += 1

    outbam.close()
    inbam.close()

def mergeBamFiles(mergedbamfn, inbamfn1, inbamfn2, inbamfn3):
    command = " ".join(["samtools merge", mergedbamfn, inbamfn1, inbamfn2, inbamfn3])
    runCommand(command)
    sortedmerged = sub('.bam$', '.sorted', mergedbamfn)
    pysam.sort(mergedbamfn, sortedmerged)
    pysam.index(sortedmerged+".bam")
    
    #remove the unsorted bam
    try:
        os.remove(mergedbamfn)
    except OSError:
        pass


def subsample(bamfn1, bamfn2, samplingrate = 0.5):
    command = " ".join(["samtools view -s", samplingrate ,"-b", bamfn1, ">", bamfn2])
    runCommand(command)
  
def sortByName(inbamfn, outbamfn):
    command = " ".join(["sambamba sort -n", inbamfn, "-o", outbamfn])
    print(command)
    runCommand(command)

def splitBed(bedfn, event):
    path, filename = os.path.split(bedfn)
    command=  "".join(["""awk '($1 ~ "chr"){print $0 >> """ ,'"{}"'.format(event), """$1".bed"}' """, bedfn])
    os.chdir(path)
    runCommand(command)
 
def generatePurities(purity):
    
    try:
         if not terminating.is_set():
            purityDir = createDirectory("/".join([finalbams_hap, str(purity)]))
            os.chdir(purityDir)
            subsample
            
    
    except (KeyboardInterrupt):
    
        logger.error('Exception Crtl+C pressed in the child process  in generation purities' )
        terminating.set()
        return
    
    except:    
    
        logger.exception("message")
        terminating.set()
        return
    return          

            
class MultiProcessingLogHandler(logging.Handler):
    def __init__(self, handler, queue, child=False):
        logging.Handler.__init__(self)

        self._handler = handler
        self.queue = queue

        # we only want one of the loggers to be pulling from the queue.
        # If there is a way to do this without needing to be passed this
        # information, that would be great!
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
        #print "receive on"
        while (self.shutdown == False) or (self.queue.empty() == False):
            # so we block for a short period of time so that we can
            # check for the shutdown cases.
            try:
                record = self.queue.get(True, self.polltime)
                self._handler.emit(record)
            except Queue.Empty, e:
                pass

    def send(self, s):
        # send just puts it in the queue for the server to retrieve
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


def initPool(queue, level, terminating_):
    """
    This causes the logging module to be initialized with the necessary info
    in pool threads to work correctly.
    """
    logging.getLogger('').setLevel(level)
    global terminating
    terminating = terminating_

def sortByNamemp(chr):

    try:
        if not terminating.is_set():
            logger.debug('sorting bam by name for '+str(chr))
            splitfile =  "/".join([splitbams, str(chr) + '.bam'])
            path, filename = os.path.split(splitfile)
        
            sortbyname = "/".join([splitbams, sub('.bam$','.byname',filename)])
            command = " ".join(["samtools view -b" ,inbam, chr ,">" , splitfile])
            runCommand(command)
            pysam.sort("-n", splitfile, sortbyname )
            
    except (KeyboardInterrupt):
        logger.error('Exception Crtl+C pressed in the child process  in sortByName for chr ' + chr)
        terminating.set()
        return
    except:    
        logger.exception("message")
    return 


def createHaplotypes(hetsnp_orig_bed, hetsnp_hap1_bed ):
    try:
        inbedh = open(hetsnp_orig_bed, 'r')
        inbedh2 = open(hetsnp_hap1_bed, 'r')  
        for line in inbedh:
            c = line.strip('\n').split("\t")
            c2 = ""
            
            while (c2 != c):
                 c2 = inbedh2.readline.strip('\n').split("\t")
                 outbedh.write(c2 +'\t'+'hap1' +'\n')
            
            outbedh.write(c2 +'\t'+'hap2' +'\n')
            
        outbedh.close()
    except:
        print('exception')
   

 
def initialize(tmpdir):
    try:
        
        if(args.phase):
            ############# THE SAME FOR ALL CANCERS #############
            logger.debug(' --- Initialization called --- ')
            vpath, vcf = os.path.split(vcfpath)
            phasedvcf = "/".join([RESULTS, sub('.vcf$', '_phased.vcf.gz', vcf)])
            vcftobed =  "/".join([RESULTS, sub('.vcf$', '.bed', vcf)])
            hap1vcf = "/".join([RESULTS,"hap1_het.vcf"])
            hap2vcf = "/".join([RESULTS, "hap2_het.vcf"])
            hap1vcffiltered = "/".join([RESULTS, "hap1_het_filtered"])
            hap2vcffiltered = "/".join([RESULTS, "hap2_het_filtered"])
            hap1vcffilteredtobed = "/".join([RESULTS, "hap1_het_filtered.bed"])
            hap2vcffilteredtobed = "/".join([RESULTS, "hap2_het_filtered.bed"])
            phased_bed =  "/".join([RESULTS, "PHASED.BED"])
              
            phaseVCF(vcfpath, phasedvcf)
            getVCFHaplotypes(phasedvcf, hap1vcf, hap2vcf)
            thinVCF(hap1vcf, hap1vcffiltered)
            thinVCF(hap2vcf, hap2vcffiltered)
            convertvcftobed(hap1vcffiltered+".recode.vcf", hap1vcffilteredtobed)
            convertvcftobed(hap2vcffiltered+".recode.vcf", hap2vcffilteredtobed)
           
            cmd1 = """sed -i 's/$/\thap1/' """+ hap1vcffilteredtobed
            cmd2 = """sed -i 's/$/\thap2/' """+ hap2vcffilteredtobed
            cmd3 = "cat " + hap1vcffilteredtobed + " " + hap2vcffilteredtobed + " > " + 'tmp.bed'
            cmd4 = "sort -V -k1,1 -k2,2 tmp.bed > " + phased_bed  
                
            runCommand(cmd1)
            runCommand(cmd2)
            runCommand(cmd3)
            runCommand(cmd4)
            os.remove('tmp.bed')  
            
            for  event in event_list: 
                roibed = "/".join([tmpdir,  event + "_roi.bed"])
                exonsinroibed = "/".join([tmpdir,   event + "_exons_in_roi.bed"])
                exonsinhetbed = "/".join([tmpdir,  event + "_exons_in_het.bed"])
                nonhetbed = "/".join([tmpdir, event + "_non_het.bed"])
                hetbed = "/".join([tmpdir, event + "_het.bed"])
                hetsnpbed = "/".join([tmpdir,  event + "_het_snp.bed"])
                
                logger.debug("Calling initialization for " + str(event))
                
                intersectBed( exons, globals()[event + 'cnv'], exonsinroibed, wa=True)
                intersectBed( exonsinroibed, phased_bed, exonsinhetbed, wa=True)  
                intersectBed(phased_bed, exonsinroibed, hetsnpbed, wa=True)
                 
                subtractBeds(exonsinroibed, phased_bed , nonhetbed) # so that 
                intersectBed( phased_bed, exonsinroibed, hetbed, wa=True)
                
                splitBed(hetsnpbed, event+'_het_snp_')
                splitBed(hetbed, event+'_het_')
                splitBed(nonhetbed, event+'_non_het_')             
            
    except:
        logger.error('Initialization error ', sys.exc_info()[0])
        raise
    return

        

def mergeSortBamFiles(mergedBamfn, finalbamdir):
    command = ""
    os.chdir(finalbamdir)
    matches = []
    
    for root,dirnames, filenames in os.walk(finalbamdir):
        for filename in fnmatch.filter(filenames, '*.bam'):
            
            path = os.path.join(root, filename)
            if os.path.islink(path):
                path = os.path.realpath(path)
                
            if (not matches.__contains__(path)):
                matches.append(path)
                command = " ".join([path, command])
            
    command2 = " ".join(["sambamba merge", mergedBamfn, command ])
    print(command2)
    runCommand(command2)
    
def main(args):
    global inbam, tmpdir, tmpbams, splitbams, finalbams, vcffiltered, vcffilteredtobed, vcffilteredtobed_padded, exons, gaincnv, losscnv, ref, logger,splittmpbams ,RESULTS,cancerDir, vcfpath, cancerType, haplotypedir
    
    vcfpath = args.vcfFile
    inbam = args.inbamFile
    exons = args.exonsFile
    ref = args.refFastaFile
    outbam = args.outBamFile
    cancerType = args.cancerType    
    gaincnv = args.cnvAmpFile
    losscnv = args.cnvDelFile
        
    path, filename = os.path.split(inbam)
    RESULTS = "/".join([os.path.abspath(os.path.dirname(__file__)), 'RESULTS'])
    createDirectory(RESULTS)
    
    #for c in cancerType:    
    cancerDir = "/".join([RESULTS, cancerType.upper()])
    createDirectory(cancerDir)
    
    path, filename = os.path.split(inbam)
    haplotypedir = "/".join([cancerDir, "haplotypedir"])
    logdir = "/".join([cancerDir, "logs"])
    splitbams = "/".join([os.path.abspath(os.path.dirname(__file__)), "splitbams"])
    tmpbams = "/".join([cancerDir, "tmpbams"])
    splittmpbams = "/".join([tmpbams,"splittmpbams"])
    finalbams = "/".join([cancerDir, "finalbams"])
    
    createDirectory(haplotypedir)
    createDirectory(logdir)
    createDirectory(tmpbams)
    createDirectory(splitbams)
    createDirectory(splittmpbams)
    createDirectory(finalbams)
        
    if( args.phase):    
        terminating = multiprocessing.Event()
        result = []
        logfile = "/".join([logdir,  "DEBUG_LOG.log"])
       
        stream = StringIO.StringIO()
        logger = logging.getLogger('')
        logger.setLevel(logging.DEBUG)
        logQueue = multiprocessing.Queue(16)
        
        filehandler = MultiProcessingLogHandler(logging.FileHandler(logfile), logQueue)
        logger.addHandler(filehandler)
        filehandler.setLevel(logging.DEBUG)
        
        initialize(haplotypedir)
        pool1 = multiprocessing.Pool(processes=16, initializer=initPool, initargs=[logQueue, logger.getEffectiveLevel(), terminating] ) 
        
        try:
            result2 = pool1.map_async(runBamgineer,  chromosome_event ).get(9999999) 
            pool1.close()
           
        except KeyboardInterrupt:  
            logger.debug('You cancelled the program!')
            pool1.terminate()
            
        except Exception as e: 
            
            logger.exception("Exception in main %s" , e)
            pool1.terminate()
            
        
        finally:
            pool1.join()
            
        time.sleep(.1)
        
        mergeSortBamFiles(outbam , finalbams )
        logging.shutdown()
        shutil.rmtree(tmpbams)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='adds CN spikes to reads, outputs modified reads as .bam along with mates')
    
    parser.add_argument('-outbam', dest='outBamFile', required=True,
                        help='.bam file name for output')  
    parser.add_argument('-ref', dest='refFastaFile', required=True,
                       help='hg19 file ')
    parser.add_argument('-vcf', dest='vcfFile', required=True,
                        help='input vcf file ')
    #without phasing
    parser.add_argument('-cnv_amp', dest='cnvAmpFile', required=False,
                        help='CNV amplification .bed file name')
    parser.add_argument('-cnv_del', dest='cnvDelFile', required=False,
                        help='CNV deletion .bed file name')
    
    parser.add_argument('-inbam', dest='inbamFile', required=False,
                        help='sam/bam file from which to obtain reads')
    parser.add_argument('-exons', dest='exonsFile', required=True,
                        help='Exon .bed file name')
    parser.add_argument('-cancertype', dest='cancerType', required=False,
                        help='acronyms for cancer type')
    
    parser.add_argument('-phase',dest= 'phase', action="store_true")
    
    
    args = parser.parse_args()
    t0 = time.time()
    main(args)
    t1 = time.time()
    print '***** Multi-processing phase took  %f ' %(t1 - t0) +' seconds to finish ******'
    logger.debug(' ***** Multi-processing phase took %f ' %(t1 - t0) +' seconds to finish ***** ')
  
