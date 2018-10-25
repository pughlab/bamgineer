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
#from pathos.multiprocessing import ProcessingPool
import signal
import itertools
import ntpath
import fnmatch
from helpers import handlers as handle

from threading import Thread
from helpers import parameters as params
import pandas as pd
from collections import defaultdict



def phaseVCF(vcfpath, phasevcfpath):
    print (" ___ phasing vcf file ___ ")
    if not vcfpath.endswith('.vcf.gz'):
        gzipFile(vcfpath)
        vcfpath = vcfpath + '.gz'

    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()

    path, vcffn = os.path.split(vcfpath)
    path2, vcffn2 = os.path.split(phasevcfpath)
    phasevcffn = sub('.vcf.gz$', '_phased', vcffn)
    command = " ".join([java_path, "-Xmx4g -jar", beagle_path, "gt=" + vcfpath, "out=" + "/".join([path2, phasevcffn]),
                        "2> beagle.log"])

    runCommand(command)
    return phasevcffn


def gzipFile(filename):
    with open(filename, 'rb') as f_in, gzip.open(filename + '.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def thinVCF(invcf, outvcf):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([vcftools_path, "--vcf", invcf, "--thin 50 --out", outvcf, "--recode"])
    runCommand(command)


def extractPairedReadfromROI(inbamfn, bedfn, outbamfn, flag="either"):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join(
        [bedtools_path, "pairtobed -abam", inbamfn, "-b", bedfn, "-type", flag, ">", outbamfn, "2> bedtool.log"])
    runCommand(command)

def extractAllReadsfromROI(inbamfn, bedfn, outbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join(
        [bedtools_path, "intersect -abam", inbamfn, "-b", bedfn, ">", outbamfn, "2> bedtool.log"])
    runCommand(command)

def extractPairedBAMfromROI(inbamfn, bedfn, outbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([samtools_path, "view -b -f 0x0001 -L", bedfn, inbamfn, ">", outbamfn])
    runCommand(command)


def dedupBam(inbamfn, outbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([samtool_path, "rmdup", inbamfn, outbamfn])
    runCommand(command)

#def insertSizeMetrics(inbamfn, metricsfn, histfn):
#    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
#    command = " ".join([java_path, "-Xmx8g -jar", picard_path, "CollectInsertSizeMetrics", "I=" + inbamfn, "O=" + metricsfn, "H=" + histfn,"M=0.5"])
#    print("*****Collecting Insert Size Metrics****")
#    runCommand(command)

   # return metricsfn

def removeDupSambamba(bamrepairedfinalsortfn, tmpbams_path=''):
    bamrepairedfinalmarkedfn = sub('.sorted.bam$', ".marked.bam", bamrepairedfinalsortfn)
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([sambamba_path, "markdup","--remove-duplicates", "--nthreads", str(4), bamrepairedfinalsortfn, bamrepairedfinalmarkedfn])
    runCommand(command)
    return bamrepairedfinalmarkedfn

def removeDupPicard(bamrepairedfinalsortfn, tmpbams_path=''):
    print (" ___ removing repaired duplicates ___ ")

    bamrepairedfinalmarkedfn = sub('.re_paired_final.sorted.bam$', ".re_paired_final.marked.bam", bamrepairedfinalsortfn)
    markedmetricsfn = sub('.re_paired_final.sorted.bam$', ".marked_metrics.txt", bamrepairedfinalsortfn) 
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([java_path, "-Xmx8g -jar", picard_path, "MarkDuplicates", "I=" + bamrepairedfinalsortfn, "O=" + bamrepairedfinalmarkedfn, "M=" + markedmetricsfn, "REMOVE_DUPLICATES=true", "ASSUME_SORTED=true"])
    runCommand(command)
    return bamrepairedfinalmarkedfn

def runCommand(cmd):
    try:
        thread = Thread(group=None, target=lambda: os.system(cmd))
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

    if phasedvcf.endswith('.vcf.gz'):
        vcfh = gzip.GzipFile(phasedvcf, 'rb')

        for line in vcfh:
            c = line.strip('\n').split("\t")
            if (len(c) == 10):
                if (c[9] == '0|1:1'):
                    out_hap1.write(line)
                    continue
                elif (c[9] == '1|0:1'):
                    out_hap2.write(line)
                    continue

            elif (line.startswith('#')):
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

        if not line.startswith('#') and len(c) >= 5 and (len(c[3]) + len(c[4]) == 2):
            start = int(c[1]) - 1
            bedh.write(c[0] + '\t' + str(start) + '\t' + str(c[1]) + '\t' + str(c[3]) + '\t' + str(
                c[4]) + '\n')  # chr start stop ref alt

    bedh.close()


# report exon bedfiles that are within the defined CNV region
def findExonsinCNVregion(cnvpath, exonpath, intersectfile, wa=False):
    cwd = os.path.dirname(__file__)
    cnvcompletepath = os.path.realpath(cnvpath.format(cwd))
    exoncompletepath = os.path.realpath(exonpath.format(cwd))

    cnvfile = pybedtools.example_bedtool(cnvcompletepath)
    exonfile = pybedtools.example_bedtool(exoncompletepath)

    f = open(intersectfile, 'w')
    if wa == False:
        print >> f, exonfile.intersect(cnvfile, u=True)
    elif wa == True:
        print >> f, exonfile.intersect(cnvfile, u=True, wa=True)

    f.close()
    return intersectfile


def subtractBeds(bedfn1, bedfn2, diffn):
    cwd = os.path.dirname(__file__)
    bed1fullpath = os.path.realpath(bedfn1.format(cwd))
    bed2fullpath = os.path.realpath(bedfn2.format(cwd))

    bed1 = pybedtools.example_bedtool(bed1fullpath)
    bed2 = pybedtools.example_bedtool(bed2fullpath)

    print (bed2fullpath + "\n" + bed2fullpath + "\n" + diffn)
    f = open(diffn, 'w')
    print >> f, bed1.subtract(bed2, A=True)
    f.close()


def bamDiff(bamfn1, bamfn2, path):
    bamutil_path = params.GetConfigReader().get('SOFTWARE', 'bamutil_path')
    command = " ".join([bamutil_path, "diff", "--in1", bamfn1, "--in2", bamfn2, "--out", "/".join([path, "diff.bam"])])
    runCommand(command)


def intersectBed(bed1fn, bed2fn, intersectfile, wa=False, wb=False):
    cwd = os.path.dirname(__file__)
    bed1fncompletepath = os.path.realpath(bed1fn.format(cwd))
    bed2fncompletepath = os.path.realpath(bed2fn.format(cwd))

    bed1 = pybedtools.example_bedtool(bed1fncompletepath)
    bed2 = pybedtools.example_bedtool(bed2fncompletepath)

    f = open(intersectfile, 'w')
    if wa == False:

        print >> f, bed1.intersect(bed2, u=True)
    elif wa == True and wb == False:
        print >> f, bed1.intersect(bed2, u=True, wa=True)
    elif wa == True and wb == True:
        print >> f, bed1.intersect(bed2, u=True, wa=True, wb=True)

    f.close()
    return intersectfile


def call_subprocess(cmd):
    try:
        proc = subprocess.Popen(cmd, shell=False)
        out, err = proc.communicate()
    except:
        logger.exception("Exception in call_subprocess ", sys.exc_info()[0])
        return


# readname function from: https://github.com/adamewing/bamsurgeon
def renamereads(inbamfn, outbamfn):
    inbam = pysam.Samfile(inbamfn, 'rb')
    outbam = pysam.Samfile(outbamfn, 'wb', template=inbam)

    paired = {}
    n = 0;
    p = 0;
    u = 0;
    w = 0;
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
                w += 1;
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


def subsample(bamfn1, bamfn2, samplingrate=0.5):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([samtools_path, "view -s", samplingrate, "-b", bamfn1, ">", bamfn2])
    runCommand(command)


def splitBamByChr(inbamfn, path, chr):
    if (chr is not None):
        java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
        command = " ".join([samtools_path, "view -bh", inbamfn, str(chr), ">", "/".join([path, str(chr) + ".bam"])])
        runCommand(command)


def sortByName(inbamfn, outbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()

    if (inbamfn is not None):
        command = " ".join([sambamba_path, "sort -n", inbamfn, "-o", outbamfn])
        print(command)
        runCommand(command)


def sortIndexBam(inbamfn, outbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([sambamba_path, "sort", inbamfn, "-o", outbamfn])
    command2 = " ".join([sambamba_path, "index", outbamfn])

    runCommand(command)
    runCommand(command2)


def sortBam(inbamfn, outbamfn, tmpbams_path=''):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([sambamba_path, "sort", inbamfn, "-o", outbamfn, '--tmpdir=', tmpbams_path])
    runCommand(command)


def getProperPairs(inbamfn, outbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([samtools_path, "view -u -h -f 0x0003", inbamfn, ">", outbamfn])
    runCommand(command)

def splitBed0(bedfn, postfix):
    path, filename = os.path.split(bedfn)
    command = "".join(["""awk '($1 ~ "chr"){print $0 >> $1 """, '"{}"'.format(postfix), """".bed"}' """, bedfn])
    os.chdir(path)
    runCommand(command)

def merge_bams(bamfn1, bamfn2, mergefn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([sambamba_path, "merge", mergefn, bamfn1, bamfn2, "--nthreads", str(4)])
    runCommand(command)

def merge_final(mergefn, finalbamdir):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    mergemarkedfn = sub('.bam$', ".marked.bam", mergefn)
    os.chdir(finalbamdir)
    command1 = " ".join([sambamba_path, "merge", mergefn, "*.bam", "--nthreads", str(4)])
    command2 = " ".join([sambamba_path, "markdup","--remove-duplicates", "--nthreads", str(4), mergefn, mergemarkedfn])
    runCommand(command1)
    print (" ___ removing merged duplicates near breakpoints ___ ")
    runCommand(command2)
    os.remove(mergefn)
    os.remove(mergefn + '.bai')
    os.rename(mergemarkedfn, mergefn)
    os.rename(mergemarkedfn + '.bai', mergefn + '.bai')

def mergeSortBamFiles(mergedBamfn, finalbamdir):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = ""
    os.chdir(finalbamdir)
    matches = []
    num_files = 0

    for root, dirnames, filenames in os.walk(finalbamdir):
        for filename in fnmatch.filter(filenames, '*.bam'):

            path = os.path.join(root, filename)
            if os.path.islink(path):
                path = os.path.realpath(path)

            if not matches.__contains__(path):
                matches.append(path)
                command = " ".join([path, command])
                num_files = num_files + 1

    if num_files > 1:
        command2 = " ".join([sambamba_path, "merge", mergedBamfn, command, "--nthreads", str(4)])
        runCommand(command2)
    elif num_files == 1:

        if str(command.strip()).endswith("GAIN.bam"):
            path, fname = os.path.split(str(command.strip()))
            inbam_original = '/'.join([params.GetSplitBamsPath(), sub('_gain', '', fname.lower())])

            command2 = " ".join([sambamba_path, "merge", mergedBamfn, command, inbam_original, "--nthreads", str(4)])
            runCommand(command2)

        elif str(command.strip()).endswith("LOSS.bam"):

            outbam = sub('.bam$', '.sort.bam', str(command.strip()))
            sortBam(command, outbam, finalbamdir)
            os.remove(str(command.strip()))


def splitPairs(inbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    pair1fn = sub('.bam$', '_read1.bam', inbamfn)
    pair2fn = sub('.bam$', '_read2.bam', inbamfn)
    command1 = " ".join([samtools_path, "view -u -h -f 0x0043", inbamfn, ">", pair1fn])
    command2 = " ".join([samtools_path, "view -u -h -f 0x0083", inbamfn, ">", pair2fn])
    runCommand(command1)
    runCommand(command2)


def getStrands(inbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    outbamfn_forward = sub('.bam$', '_forward.bam', inbamfn)
    outbamfn_reverse = sub('.bam$', '_reverse.bam', inbamfn)
    command1 = " ".join([samtools_path, "view -F 0x10", inbamfn, ">", outbamfn_forward])
    command2 = " ".join([samtools_path, "view -f 0x10", inbamfn, ">", outbamfn_reverse])
    runCommand(command1)
    runCommand(command2)


def countReads(inbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    cmd = " ".join([samtools_path, "view", inbamfn, "|wc -l"])
    out, err = subprocess.Popen(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                stdin=subprocess.PIPE, shell=True).communicate()
    return "".join(out.split())


def find_unpaired_reads(inbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    unpairedfn = sub('.bam$', '.unpairedfn.bam', inbamfn)
    command1 = " ".join([samtools_path, "view -u -h -f  0x0004", inbamfn, ">", unpairedfn])
    runCommand(command1)


def splitPairAndStrands(inbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    read1_strand1sortfn = sub('.bam$', '.read1_pos.bam', inbamfn)
    read1_strand2sortfn = sub('.bam$', '.read1_neg.bam', inbamfn)
    read2_strand1sortfn = sub('.bam$', '.read2_pos.bam', inbamfn)
    read2_strand2sortfn = sub('.bam$', '.read2_neg.bam', inbamfn)

    mapped_all = sub('sorted.bam$', 'mapped_all.bam', inbamfn)

    command1 = " ".join([samtools_path, "view -u -h -f 0x0061", inbamfn, ">", read1_strand1sortfn])
    command2 = " ".join([samtools_path, "view -u -h -f 0x0051", inbamfn, ">", read1_strand2sortfn])
    command3 = " ".join([samtools_path, "view -u -h -f 0x0091", inbamfn, ">", read2_strand1sortfn])
    command4 = " ".join([samtools_path, "view -u -h -f 0x00A1", inbamfn, ">", read2_strand2sortfn])

    runCommand(command1)
    runCommand(command2)
    runCommand(command3)
    runCommand(command4)


def splitStrands(inbamfn):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    read_strand1sortfn = sub('.bam$', '.read_pos.bam', inbamfn)
    read_strand2sortfn = sub('.bam$', '.read_neg.bam', inbamfn)

    mapped_all = sub('sorted.bam$', 'mapped_all.bam', inbamfn)

    command1 = " ".join([samtools_path, "view -u -h -f 33", inbamfn, ">", read_strand1sortfn])
    command2 = " ".join([samtools_path, "view -u -h -f 17", inbamfn, ">", read_strand2sortfn])

    runCommand(command1)
    runCommand(command2)

def extract_proper_paired_reads(inbamfn, properfn):
    # properfn = sub('.bam$', '_proper.bam', inbamfn)
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    command = " ".join([samtools_path, "view -f 0x03 -bSq 30", inbamfn, ">", properfn])
    runCommand(command)
    os.remove(inbamfn)


def removeIfEmpty(bamdir, file):
    java_path, beagle_path, picard_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    if file.endswith(".bam"):
        command = " ".join([samtools_path, "view", "/".join([bamdir, file]), "| less | head -1 | wc -l"])
        nline = subprocess.check_output(command, shell=True)
        if os.path.isfile("/".join([bamdir, file])) and (int(nline) == 0):
            os.remove("/".join([bamdir, file]))


def createHaplotypes(hetsnp_orig_bed, hetsnp_hap1_bed):
    try:
        inbedh = open(hetsnp_orig_bed, 'r')
        inbedh2 = open(hetsnp_hap1_bed, 'r')
        for line in inbedh:
            c = line.strip('\n').split("\t")
            c2 = ""

            while c2 != c:
                c2 = inbedh2.readline.strip('\n').split("\t")
                outbedh.write(c2 + '\t' + 'hap1' + '\n')

            outbedh.write(c2 + '\t' + 'hap2' + '\n')

        outbedh.close()
    except:
        print('exception')


def create_chr_event_list(cnv_list, chr_list):
    chrom_event = []
    for c in chr_list:
        for cnv_path in cnv_list:
            e = os.path.splitext(ntpath.basename(cnv_path))[0]
            chev = "_".join([str(c), str(e)])
            chrom_event.append(chev)
    return chrom_event


def createEventBedFiles(cnv_dir, bedfn):
    cnv_number_list = []
    event_list = []
    df = pd.read_csv(bedfn, header=None, sep='\t')
    #cnv_number_list = list(set(df[df.columns[-1]].tolist()))

    df.columns = ['chr', 'start_pos', 'end_pos', 'hap_type', 'abs_cn']
    cnv_number_list = list(zip(df.abs_cn, df.hap_type, df.start_pos))
    #print cnv_number_list
    for num in cnv_number_list:
        cn = int(num[0])
	hap = str(num[1])
	start = int(num[2])
        dfi = df.loc[(df['abs_cn'] == cn) & (df['hap_type'] == hap) & (df['start_pos'] == start)]

        if cn == 0:
            fn = 'deepdel.bed'
        elif cn == 1 and len(hap) == cn:
            fn = 'loss'+hap+str(start)+'.bed'
        elif cn == 2 and len(hap) == cn:
            fn = 'loh'+hap+str(start)+'.bed'
        elif cn == 3 and len(hap) == cn:
            fn = 'gain'+hap+str(start)+'.bed'
        elif cn == 4 and len(hap) == cn:
            fn = 'amp4'+hap+str(start)+'.bed'
        elif cn == 5 and len(hap) == cn:
            fn = 'amp5'+hap+str(start)+'.bed'
        elif cn == 6 and len(hap) == cn:
            fn = 'amp6'+hap+str(start)+'.bed'
        elif cn == 7 and len(hap) == cn:
            fn = 'amp7'+hap+str(start)+'.bed'
        elif cn == 8 and len(hap) == cn:
            fn = 'amp8'+hap+str(start)+'.bed'
        elif cn == 9 and len(hap) == cn:
            fn = 'amp9'+hap+str(start)+'.bed'
        elif cn == 10 and len(hap) == cn:
            fn = 'amp10'+hap+str(start)+'.bed'
        elif cn == 11 and len(hap) == cn:
            fn = 'amp11'+hap+str(start)+'.bed'
        elif cn == 12 and len(hap) == cn:
            fn = 'amp12'+hap+str(start)+'.bed'
        elif cn == 13 and len(hap) == cn:
            fn = 'amp13'+hap+str(start)+'.bed'
        elif cn == 14 and len(hap) == cn:
            fn = 'amp14'+hap+str(start)+'.bed'
        elif cn == 15 and len(hap) == cn:
            fn = 'amp15'+hap+str(start)+'.bed'
        elif cn > 15:
            print('CNV number must be smaller than 15')
	elif len != cn:
	    print('Allelic ratio must match CNV# (i.e. AAB = CNV of 3)')

        dfi.to_csv("/".join([cnv_dir, fn]), sep='\t', header=None, encoding='utf-8', index=False)

    return

def splitBedByChr(cnvbedfn, hap_dir):
    chr_bedlist = []
    #event_list = []
    df = pd.read_csv(cnvbedfn, header=None, sep='\t')
    chr_bedlist = list(set(df[df.columns[0]].tolist()))
    df.columns = ['chr', 'start', 'end'] #, 'hap_type', 'abs_cn']
    for chromosome in chr_bedlist:
        dfi = df.loc[(df['chr'] == chromosome)]
	fn = chromosome+'_non_roi.bed'
	#else statement for len != cn
        dfi.to_csv("/".join([hap_dir, fn]), sep='\t', header=None, encoding='utf-8', index=False)

    return

def removeIfEmptyBed(cnvbedfn):
    if os.path.getsize(cnvbedfn) == 0:
    	os.remove(cnvbedfn)	

def generatePhasedBed(hap1vcffilteredtobed, hap2vcffilteredtobed, phased_bed):
    print (" ___ generating phased bed ___ ")
    df1 = pd.read_csv(hap1vcffilteredtobed, header=None, sep='\t')
    df2 = pd.read_csv(hap2vcffilteredtobed, header=None, sep='\t')

    df1[len(df1.columns)] = "hap1"
    df2[len(df2.columns)] = "hap2"

    df = df1.append(df2)

    df = df.sort_values([0, 1], ascending=True)
    df.to_csv(phased_bed, sep='\t', header=None, encoding='utf-8', index=False)


def filterColumns(inp, outp, cols):

    print (" ___ filtering bed file columns for " + os.path.basename(inp)+" ___")
    df1 = pd.read_csv(inp, header=None, sep='\t')
    last_col = len(df1.columns) - 1
    cols.append(last_col)
    df2 = df1.filter(cols)
    df2.to_csv(outp, sep='\t', header=None, encoding='utf-8', index=False)

def splitBed(bedfn, postfix):

    path, filename = os.path.split(bedfn)
    df = pd.read_csv(bedfn, header=None, sep='\t')

    chroms = df[0].unique()
    for chr in chroms:
        outfilename = str(chr) + postfix + '.bed'
        outp = '/'.join([path, outfilename])
        df[df[0] == chr].to_csv(outp, sep='\t', header=None, encoding='utf-8', index=False)

def read_pair_generator(bam, chromosome, start_position, end_position):
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(chromosome, start_position, end_position):
        if not read.is_paired:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]
