#!/usr/bin/env python
import sys
import pybedtools
import pysam
import os
import subprocess
from uuid import uuid4
from re import sub
from itertools import izip

import gzip
import shutil
import traceback
import multiprocessing
from multiprocessing import Pool
from contextlib import closing
from pathos.multiprocessing import ProcessingPool 
import signal
import itertools
import ntpath
import fnmatch

from threading import Thread
from helpers import parameters as params


def phaseVCF(vcfpath, phasevcfpath):
    print (" ___ phasing vcf file ___ "  )    
    if(not vcfpath.endswith('.vcf.gz')):
        gzipFile(vcfpath)
        vcfpath = vcfpath+'.gz'
    
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    path, vcffn = os.path.split(vcfpath)
    path2, vcffn2 = os.path.split(phasevcfpath)
    phasevcffn = sub('.vcf.gz$', '_phased', vcffn)
    command = " ".join([java_path,"-Xmx4g -jar", beagle_path, "gt="+vcfpath, "out="+"/".join([path2, phasevcffn]), "2> beagle.log"])
    
    runCommand(command)
    return phasevcffn

def create_chr_bam_list():
    chrom_event= []
    for c in range(1,23):
        for e in ['nonhet','mutated']:
            chev = "_".join(['chr'+str(c), e])
            chrom_event.append(chev)
    return chrom_event

def gzipFile(filename):
    with open(filename, 'rb') as f_in, gzip.open(filename+'.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def thinVCF(invcf, outvcf):
   java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
   command = " ".join([vcftools_path,"--vcf", invcf, "--thin 50 --out", outvcf,  "--recode"])
   runCommand(command)

def extractPairedReadfromROI(inbamfn, bedfn, outbamfn, flag = "either"):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    command = " ".join([bedtools_path,"pairtobed -abam", inbamfn,"-b",bedfn, "-type", flag,">", outbamfn, "2> bedtool.log"])
    runCommand(command)
    
def extractBAMfromROI_All(inbamfn, bedfn, outbamfn):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    command = " ".join([samtool_path, "view -b -L", bedfn, inbamfn, ">",outbamfn])
    runCommand(command)

def dedupBam(inbamfn, outbamfn):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    command = " ".join([samtool_path, "rmdup", inbamfn, outbamfn])
    runCommand(command)

def runCommand(cmd):
    try:
        thread = Thread(group=None, target=lambda:os.system(cmd))
        thread.run()
        if not thread.is_alive():
            return 0
        else:
            return 1
    except OSError as e:
        sys.exit(1)


def getVCFHaplotypes(phasedvcf, hap1, hap2):
    out_hap1 = open(hap1, 'w')
    out_hap2 = open(hap2, 'w')
    
    if(phasedvcf.endswith('.vcf.gz')):
        vcfh = gzip.GzipFile(phasedvcf, 'rb')
          
        for line in vcfh:
            c = line.strip('\n').split("\t")
            if (len(c) == 10 ):
                if(c[9] == '0|1:1'):
                    out_hap1.write(line)
                    continue
                elif(c[9] == '1|0:1'):
                     out_hap2.write(line)
                     continue
                
            elif(line.startswith('#')):
                out_hap1.write(line)
                out_hap2.write(line)
                continue
        
    out_hap1.close()
    out_hap2.close()
   

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

def bamDiff(bamfn1, bamfn2, path):
    bamutil_path = params.GetConfigReader().get('SOFTWARE', 'bamutil_path')
    command = " ".join([bamutil_path,"diff", "--in1", bamfn1, "--in2", bamfn2, "--out" ,"/".join([path,"diff.bam"])]) 
    runCommand(command ) 

def intersectBed(bed1fn, bed2fn, intersectfile, wa=False):
    cwd = os.path.dirname(__file__)
    bed1fncompletepath = os.path.realpath(bed1fn.format(cwd))
    bed2fncompletepath = os.path.realpath(bed2fn.format(cwd))
      
    bed1 = pybedtools.example_bedtool(bed1fncompletepath)
    bed2 = pybedtools.example_bedtool(bed2fncompletepath)
      
    f = open(intersectfile, 'w')
    if(wa == False):
        
        print >> f, bed1.intersect(bed2,u=True)
    elif(wa==True):
        print >> f, bed1.intersect(bed2,u=True,wa=True)
            
    f.close()
    return intersectfile

def call_subprocess(cmd):
    try:
        proc = subprocess.Popen(cmd, shell=False)
        out, err = proc.communicate()
    except:    
        logger.exception("Exception in call_subprocess " , sys.exc_info()[0])
        return
    
 

#readname function from: https://github.com/adamewing/bamsurgeon    
def renamereads(inbamfn, outbamfn):
    
    inbam = pysam.Samfile(inbamfn, 'rb')
    outbam = pysam.Samfile(outbamfn, 'wb', template=inbam)

    paired = {}
    n = 0;p = 0;u = 0;w = 0;m = 0
    
    for read in inbam.fetch(until_eof=True):
        n += 1
        if read.is_paired:
            p += 1
            if read.qname in paired:
                uuid = paired[read.qname]
                del paired[read.qname]
                read.qname = uuid
                outbam.write(read)
                w += 1;m += 1
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


def subsample(bamfn1, bamfn2, samplingrate = 0.5):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    command = " ".join([samtools_path,"view -s", samplingrate ,"-b", bamfn1, ">", bamfn2])
    runCommand(command)

def splitBamByChr(inbamfn, path,chr):
    java_path, beagle_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([samtools_path, "view -bh", inbamfn, chr, ">",  "/".join[path,chr+".bam"]])
    #sortByName("/".join[path,chr+".bam"], "/".join[path,chr+".byname.bam"])
    print(command)
    runCommand(command)

def sortByName(inbamfn, outbamfn):
    java_path, beagle_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([samtools_path, "sort -n", inbamfn, "-o", outbamfn])
    print(command)
    runCommand(command)
    
def sortIndexBam(inbamfn, outbamfn):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    command = " ".join([sambamba_path, "sort", inbamfn, "-o", outbamfn])
    command2 = " ".join([sambamba_path,"index", outbamfn])
    
    runCommand(command)
    runCommand(command2)

def sortBam(inbamfn, outbamfn):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    command = " ".join([sambamba_path, "sort", inbamfn, "-o", outbamfn])
    runCommand(command)
  
def getProperPairs(inbamfn, outbamfn):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    command = " ".join([samtools_path,"view -u -h -f 0x0003", inbamfn ,">", outbamfn])
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

def merge_bams(bamfn1, bamfn2, mergefn):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    command = " ".join([sambamba_path, "merge", mergefn, bamfn1, bamfn2 ])
    runCommand(command)

def mergeSortBamFiles(mergedBamfn, finalbamdir):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
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
            
    command2 = " ".join([sambamba_path, "merge", mergedBamfn, command ])
    runCommand(command2)
    
def getMeanSTD(inbam):
    """ awk '{ if ($9 > 0) { N+=1; S+=$9; S2+=$9*$9 }} END { M=S/N; print "n="N", mean="M", stdev="sqrt ((S2-M*M*N)/(N-1))}' """
    command = " ".join([awk, inbam])
    runCommand(command )


def splitPairs(inbamfn):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    pair1fn =  sub('.bam$', '_read1.bam', inbamfn)
    pair2fn =  sub('.bam$', '_read2.bam', inbamfn)
    command1 = " ".join([samtools_path, "view -u -h -f 0x0043", inbamfn, ">", pair1fn])
    command2 = " ".join([samtools_path ,"view -u -h -f 0x0083", inbamfn, ">", pair2fn])
    runCommand(command1)
    runCommand(command2)

def getStrands(inbamfn):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    outbamfn_forward =  sub('.bam$', '_forward.bam', inbamfn)
    outbamfn_reverse =  sub('.bam$', '_reverse.bam', inbamfn)
    command1 = " ".join([samtools_path,"view -F 0x10", inbamfn, ">",outbamfn_forward])
    command2 = " ".join([samtools_path,"view -f 0x10", inbamfn, ">",outbamfn_reverse])    
    runCommand(command1)
    runCommand(command2)
    
    
def countReads(inbamfn):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    cmd = " ".join([samtools_path, "view", inbamfn, "|wc -l"])
    out,err= subprocess.Popen(cmd, 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE, 
                            stdin=subprocess.PIPE, shell=True).communicate()
    return "".join(out.split())

def find_unpaired_reads(inbamfn):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath()
    unpairedfn =  sub('.bam$', '.unpairedfn.bam', inbamfn)
    command1 = " ".join([samtools_path ,"view -u -h -f  0x0004", inbamfn, ">", unpairedfn])
    runCommand(command1)

def splitPairAndStrands(inbamfn):
    
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    read1_strand1sortfn =  sub('.bam$', '.read1_pos.bam', inbamfn)
    read1_strand2sortfn =  sub('.bam$', '.read1_neg.bam', inbamfn)
    read2_strand1sortfn =  sub('.bam$', '.read2_pos.bam', inbamfn)
    read2_strand2sortfn =  sub('.bam$', '.read2_neg.bam', inbamfn)
    
    read_comp =  sub('.bam$', '.complementary.bam', inbamfn)
    mapped_all  =  sub('sorted.bam$', 'mapped_all.bam', inbamfn)
    
    command1 = " ".join([samtools_path ,"view -u -h -f 0x0061", inbamfn, ">", read1_strand1sortfn])
    command2 = " ".join([samtools_path ,"view -u -h -f 0x0051", inbamfn, ">", read1_strand2sortfn])
    command3 = " ".join([samtools_path ,"view -u -h -f 0x0091", inbamfn, ">", read2_strand1sortfn])
    command4 = " ".join([samtools_path ,"view -u -h -f 0x00A1", inbamfn, ">", read2_strand2sortfn])
    

    command5 = " ".join([sambamba_path ,"""view -F "unmapped or mate_is_unmapped or secondary_alignment or not (paired) or duplicate" """, "-f bam" ,inbamfn ,">", read_comp])
   
    runCommand(command1)
    runCommand(command2)
    runCommand(command3)
    runCommand(command4)
    
    runCommand(command5)


def extract_proper_paired_reads(inbamfn, properfn):
    #properfn = sub('.bam$', '_proper.bam', inbamfn)
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath() 
    command = " ".join([samtools_path,"view -f 0x03 -bSq 30", inbamfn, ">", properfn])
    runCommand(command)
    os.remove(inbamfn)
    
def removeIfEmpty(bamdir,file):
    java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path  = params.GetSoftwarePath()  
    if file.endswith(".bam"):
       command = " ".join([samtools_path, "view", "/".join([bamdir, file]), "| less | head -1 | wc -l" ])
       nline= subprocess.check_output(command, shell = True)  
       if (os.path.isfile( "/".join([bamdir, file])) and (int(nline) == 0)):
               os.remove("/".join([bamdir, file]))
                      
     
    
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
   
#chr 21 and 22 for test, change it to 1
def create_chr_event_list():
    chrom_event= []
    for c in range(1,23):
        for e in ['gain','loss']:
            chev = "_".join(['chr'+str(c), e])
            chrom_event.append(chev)
    return chrom_event








