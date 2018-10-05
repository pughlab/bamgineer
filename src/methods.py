import csv
import glob
import logging
import time
import random
from uuid import uuid4
from helpers import handlers as handle
from helpers import parameters as params
from helpers import bamgineerHelpers as bamhelp
from utils import *



global bases
bases = ('A', 'T', 'C', 'G')


def initPool(queue, level, terminating_):
    # This causes the logging module to be initialized with the necessary info in pool threads
    logging.getLogger('').setLevel(level)
    global terminating
    terminating = terminating_


def initialize0(results_path, cancer_dir_path):
    """
    Initialize paths for files provided by user.
    Phasing is also performed if requested by user.
    """
    try:
        vcf_path = bamhelp.GetVCF()
        exons_path = bamhelp.GetExons()
        reference_path = bamhelp.GetRef()
        bedtools_path = bamhelp.GetBedtoolsPath()
        vpath, vcf = os.path.split(vcf_path)

        if params.GetPhase():
            phasedvcf = "/".join([results_path, sub('.vcf$', '_phased.vcf.gz', vcf)])
            vcftobed = "/".join([results_path, sub('.vcf$', '.bed', vcf)])
            hap1vcf = "/".join([results_path, "hap1_het.vcf"])
            hap2vcf = "/".join([results_path, "hap2_het.vcf"])
            hap1vcffiltered = "/".join([results_path, "hap1_het_filtered"])
            hap2vcffiltered = "/".join([results_path, "hap2_het_filtered"])
            hap1vcffilteredtobed = "/".join([results_path, "hap1_het_filtered.bed"])
            hap2vcffilteredtobed = "/".join([results_path, "hap2_het_filtered.bed"])
            phased_bed = "/".join([results_path, "PHASED.BED"])

            phaseVCF(vcf_path, phasedvcf)
            getVCFHaplotypes(phasedvcf, hap1vcf, hap2vcf)
            convertvcftobed(hap1vcf, hap1vcffilteredtobed)
            convertvcftobed(hap2vcf, hap2vcffilteredtobed)
            generatePhasedBed(hap1vcffilteredtobed, hap2vcffilteredtobed, phased_bed)

    except:

        logger.exception("Initialization error !")
        raise

    return


def initialize_pipeline(phase_path, haplotype_path, cnv_path):
    """
    Intersect CNV region with exons and then intersect with phased VCF bed.
    """
    exons_path = bamhelp.GetExons()
    cnv_bed = params.GetCNV()

    event, extension = os.path.splitext(os.path.basename(cnv_path))

    phased_bed = "/".join([phase_path, "PHASED.BED"])
    nonroibedfn = "/".join([haplotype_path, "non_roi.bed"])
    bedtools_path = bamhelp.GetBedtoolsPath()

    try:
        logger.debug(' --- Initializing input files  --- ')
        exonsinroibed = "/".join([haplotype_path, "exons_in_roi" + str(event) + ".bed"])
        nonhetbed = "/".join([haplotype_path, "non_het" + str(event) + ".bed"])
        hetbed = "/".join([haplotype_path, "het" + str(event) + ".bed"])
        hetsnpbed = "/".join([haplotype_path, "het_snp" + str(event) + ".bed"])
        tmp1 = "/".join([haplotype_path, str(event) + "_tmp1.bed"])
        tmp2 = "/".join([haplotype_path, str(event) + "_tmp2.bed"])
    
        command = " ".join([bedtools_path, "intersect -a", cnv_path, "-b", exons_path, " > ", exonsinroibed])
        runCommand(command)
        splitBed(exonsinroibed, '_exons_in_roi' + str(event))
        command = " ".join([bedtools_path, "intersect -a", phased_bed, "-b", exonsinroibed, "-wa -wb >", tmp2])
        runCommand(command)
        removeIfEmptyBed(tmp2)
        if os.path.isfile(tmp2):
            filterColumns(tmp2, hetsnpbed, [0,1,2,3,4,5,9])#[i for i in range(0, 8)])
            splitBed(hetsnpbed, '_het_snp' + str(event))
    
        # non-roi bed:
        if not os.path.isfile(nonroibedfn):
            command = " ".join([bedtools_path, "subtract -a", exons_path, "-b", cnv_bed, ">", nonroibedfn])
            runCommand(command)
            removeIfEmptyBed(nonroibedfn)
            if os.path.isfile(nonroibedfn):
                splitBedByChr(nonroibedfn, haplotype_path) 
         
    except:
        logger.exception("Initialization error !")
        raise
    logger.debug("--- initialization complete ---")
    return

def init_file_names(chr, tmpbams_path, haplotypedir, event):
    """
    Initialize file names for: 
    ROI bam, chromosome sorted by name, chromosome sorted by coordinate, and heterozygous SNPs bed.
    """
    flist = []

    roibam = "/".join([tmpbams_path, chr + "_roi" + event + ".bam"])
    splitbams = params.GetSplitBamsPath()
    hetsnp = "/".join([haplotypedir, chr + '_het_snp' + event + '.bed'])

    if not splitbams:
        splitbams = "/".join([res_path, 'splitbams'])

    sortbyname = "/".join([splitbams, chr + '.byname.bam'])
    sortbyCoord = "/".join([splitbams, chr + '.bam'])

    flist.extend([roibam, sortbyname, sortbyCoord, hetsnp])
    return flist


def find_roi_bam(chromosome_event):
    """
    Extract paired reads from original bam using the user specificed ROI (cnv.bed). 
    """
    chr, event = chromosome_event.split("_")
    roi, sortbyname, sortbyCoord, hetsnp = init_file_names(chr, tmpbams_path, haplotype_path, event)
    exonsinroibed = "/".join([haplotype_path, chr + "_exons_in_roi" + event + ".bed"])
 
    success = False
    try:
        if not terminating.is_set():
            roisort = sub('.bam$', '.sorted', roi)
            if os.path.isfile(exonsinroibed):
                cmd = " ".join(["sort -u", exonsinroibed, "-o", exonsinroibed]);
                runCommand(cmd)
                print(" ___ extracting roi bams  ___")
                extractPairedReadfromROI(sortbyname, exonsinroibed, roi)
                removeIfEmpty(tmpbams_path, ntpath.basename(roi))
                pysam.sort(roi, roisort)
                pysam.index(roisort + '.bam')
                os.remove(roi)
                success = True

            else:
                logger.debug(exonsinroibed + ' does not exist!')
                return

    except (KeyboardInterrupt):
        logger.error('Exception Crtl+C pressed in the child process  in find_roi_bam for chr ' + chr + event)
        terminating.set()
        success = False
        return

    except Exception as e:
        logger.exception("Exception in find_roi_bam %s", e)
        terminating.set()
        success = False
        return

    if (success):
        logger.debug("find_roi_bam complete successfully for " + chr + event)

    return

def find_non_roi_bam(chr_list):
    """
    Extract paired reads from original bam using generated non-ROI bed. 
    """
    chr = chr_list
    splitbams = params.GetSplitBamsPath()
    
    sortbyname = "/".join([splitbams, chr + '.byname.bam'])
    sortbyCoord = "/".join([splitbams, chr + '.bam'])
    nonroi = "/".join([finalbams_path, chr + "_non_roi.bam"])
    exonsnonroibed = "/".join([haplotype_path, chr + "_non_roi.bed"])

    success = False
    try:
        if not terminating.is_set():
            nonroisort = sub('.bam$', '.sorted', nonroi)
            if os.path.isfile(nonroisort):
                success = True

            else:
                if os.path.isfile(exonsnonroibed):

                    cmd = " ".join(["sort -u", exonsnonroibed, "-o", exonsnonroibed]);
                    runCommand(cmd)
                    print(" ___ extracting non-roi bams  ___")
                    extractAllReadsfromROI(sortbyCoord, exonsnonroibed, nonroi)
                    removeIfEmpty(finalbams_path, ntpath.basename(nonroi))
                    pysam.sort(nonroi, nonroisort)
                    pysam.index(nonroisort + '.bam')
                    os.remove(nonroi)
                    success = True

                else:
                    logger.debug(exonsnonroibed + ' does not exist!')
                    return

    except (KeyboardInterrupt):
        logger.error('Exception Crtl+C pressed in the child process  in find_roi_bam for chr ' + chr)
        terminating.set()
        success = False
        return
    
    except Exception as e:
        logger.exception("Exception in find_non_roi_bam %s", e)
        terminating.set()
        success = False
        return
    
    if (success):
        logger.debug("find_non_roi_bam complete successfully for " + chr)
    
    return


def modify_hap(bamsortfn, hap1_bamsortfn, hap2_bamsortfn):
    """
    Subfunction for split_hap function.
    Merges  hap1 SNP reads and hap2 SNP reads.
    Finds other SNPs (homozygous) not found in hap1 and hap2. 
    Splits the homozygous SNP reads 50/50.
    """
    hap12_bamfn = sub('.sorted.bam$', ".hap12.bam", bamsortfn)
    hap12_bamsortfn = sub('.sorted.bam$', ".hap12.sorted", bamsortfn)


    intersnpbamfn = sub('.sorted.bam$', ".intersnp.bam", bamsortfn)
    intersnpbamsortfn = sub('.sorted.bam$', ".intersnp.sorted", bamsortfn)


    hap1_intersnpbamfn = sub('.sorted.bam$', ".hap1_intersnp.bam", bamsortfn)
    hap2_intersnpbamfn = sub('.sorted.bam$', ".hap2_intersnp.bam", bamsortfn) 
    hap1_intersnpbamsortfn = sub('.sorted.bam$', ".hap1_intersnp.sorted", bamsortfn)
    hap2_intersnpbamsortfn = sub('.sorted.bam$', ".hap2_intersnp.sorted", bamsortfn)        

    
    # merge hap1 and hap2 into hap12 and sort 
    merge_bams(hap1_bamsortfn + '.bam', hap2_bamsortfn + '.bam', hap12_bamfn)
    sortBam(hap12_bamfn, hap12_bamsortfn + '.bam', tmpbams_path)
    
    # find difference between normal bam and hap12 (intersnps) and sort
    bamDiff(bamsortfn, hap12_bamsortfn + '.bam', tmpbams_path)
    sortBam("/".join([tmpbams_path, 'diff_only1_' + os.path.basename(bamsortfn)]), intersnpbamsortfn + '.bam', tmpbams_path)

    # subsample 50% of reads from inter_snps and assign to hap1 and sort
    subsample(intersnpbamsortfn + '.bam', hap1_intersnpbamfn, str(0.5))
    sortBam(hap1_intersnpbamfn, hap1_intersnpbamsortfn + '.bam', tmpbams_path)
 
    # find difference between hap1 inter_snp bam and total inter_snp bam and sort 
    bamDiff(intersnpbamsortfn + '.bam', hap1_intersnpbamsortfn + '.bam', tmpbams_path)
    sortBam("/".join([tmpbams_path, 'diff_only1_' + os.path.basename(intersnpbamsortfn + '.bam')]), hap2_intersnpbamsortfn + '.bam', tmpbams_path)
    
    # remove intermediate files
    os.remove(hap12_bamfn)
    os.remove(hap12_bamfn + '.bai')
    os.remove(hap12_bamsortfn + '.bam')
    os.remove(hap12_bamsortfn + '.bam.bai')
    os.remove(intersnpbamsortfn + '.bam')
    os.remove(intersnpbamsortfn + '.bam.bai')
    os.remove(hap1_intersnpbamfn)
    os.remove("/".join([tmpbams_path, 'diff_only1_' +  os.path.basename(intersnpbamsortfn + '.bam')]))
    os.remove("/".join([tmpbams_path, 'diff_only1_' +  os.path.basename(bamsortfn)]))
    #os.remove(intersnpbamsortfn + '.bam')
    #os.remove('diff_only1_' + intersnpbamsortfn + '.bam.bai')

    return hap1_intersnpbamsortfn, hap2_intersnpbamsortfn

def merge_haps(bamsortfn, hap1_bamsortfn, hap2_bamsortfn, hap1_intersnpbamsortfn, hap2_intersnpbamsortfn):
    """
    Subfunction for split_hap function.
    Merges hap1 and hap2 to each 50% split of homozygous SNP reads.
    """

    hap1_finalbamfn = sub('.sorted.bam$', ".hap1_final.bam", bamsortfn)
    hap2_finalbamfn = sub('.sorted.bam$', ".hap2_final.bam", bamsortfn)
    hap1_finalbamsortfn = sub('.sorted.bam$', ".hap1_final.sorted", bamsortfn)
    hap2_finalbamsortfn = sub('.sorted.bam$', ".hap2_final.sorted", bamsortfn)
    
    # merge hap1 with hap1 intersnps and hap2 with hap2 intersnps
    merge_bams(hap1_bamsortfn + '.bam', hap1_intersnpbamsortfn + '.bam', hap1_finalbamfn)
    merge_bams(hap2_bamsortfn + '.bam', hap2_intersnpbamsortfn + '.bam', hap2_finalbamfn)

    # sort final bams
    sortBam(hap1_finalbamfn, hap1_finalbamsortfn + '.bam', tmpbams_path)
    sortBam(hap2_finalbamfn, hap2_finalbamsortfn + '.bam', tmpbams_path)
    
    # remove intermediate files
    #os.remove(hap1_intersnpbamsortfn + '.bam')
    #os.remove(hap2_intersnpbamsortfn + '.bam')
    os.remove(hap1_finalbamfn)
    os.remove(hap2_finalbamfn)
    os.remove(hap1_finalbamfn + '.bai')
    os.remove(hap2_finalbamfn + '.bai')
 
    return hap1_finalbamsortfn, hap2_finalbamsortfn

def split_hap(bamsortfn, chr, event):
    """
    Phased VCF bed file read in specifying haplotype and ref/var bases at SNP position.
    Fetch reads overlapping every het. SNP position found in phased bed.
    Write reads to hap1 and hap2 bams based on base at SNP position on read: 
    Variant base at SNP position sent to specified haplotype and read with reference base at SNP position sent to other haplotype.
    Only write reads if read id does not exist in read 1 set or read 2 set.
    """
    print(" ___ splitting original bam into hap1 and hap2 ___")

    fn, sortbyname, sortbyCoord, bedfn = init_file_names(chr, tmpbams_path, haplotype_path, event)

    hap1_bamfn = sub('.sorted.bam$', ".hap1.bam", bamsortfn)
    hap2_bamfn = sub('.sorted.bam$', ".hap2.bam", bamsortfn)
    hap1_bamsortfn = sub('.sorted.bam$', ".hap1.sorted", bamsortfn)
    hap2_bamsortfn = sub('.sorted.bam$', ".hap2.sorted", bamsortfn)
    hap1_finalbamsortfn = sub('.sorted.bam$', ".hap1_final.sorted", bamsortfn)
    hap2_finalbamsortfn = sub('.sorted.bam$', ".hap2_final.sorted", bamsortfn)


    if not terminating.is_set():
        if not os.path.isfile(bamsortfn):
            raise ValueError('Could not find file bamsortfn')

        if (os.path.isfile(bedfn)):
            cmd = " ".join(["sort -u", bedfn, "-o", bedfn]);
            runCommand(cmd)

            samfile = pysam.Samfile(bamsortfn, "rb")
            alignmentfile = pysam.AlignmentFile(bamsortfn, "rb")
            outbam1 = pysam.Samfile(hap1_bamfn, 'wb', template=samfile)
            outbam2 = pysam.Samfile(hap2_bamfn, 'wb', template=samfile)

            bedfile = open(bedfn, 'r')
            readids1 = set()
            readids2 = set()

            for bedline in bedfile:
                c = bedline.strip().split()
    
                if len(c) == 8:
                    chr2 = c[0]
                    chr = c[0].strip("chr")
                    start = int(c[1])
                    end = int(c[2])
                    refbase = str(c[3])
                    altbase = str(c[4])
                    haplotype = str(c[5])
                    copy_number = int(c[7])
                else:
                    continue
        
                readmappings = alignmentfile.fetch(chr2, start, end)

                for shortread in readmappings:
                    shortreadid = shortread.qname 
                    if shortread.is_read1 and shortreadid not in readids1 or shortread.is_read2 and shortreadid not in readids2:  

                        problem_with_read = False

                        try:
                            index = shortread.get_reference_positions(full_length=True).index(start)
                            tmpread = shortread.query_sequence
                            qual = shortread.query_qualities
                            tmpread_index = tmpread[index]

                            if shortread.is_read1 and shortreadid not in readids1:
                                readids1.add(shortreadid)
                            elif shortread.is_read2 and shortreadid not in readids2:
                                readids2.add(shortreadid)
             
                            if tmpread_index == altbase and haplotype == "hap1":
                                outbam1.write(shortread)
                            elif tmpread_index == refbase and haplotype == "hap1":
                                outbam2.write(shortread)
                            elif tmpread_index == altbase and haplotype == "hap2":
                                outbam2.write(shortread)
                            elif tmpread_index == refbase and haplotype == "hap2":
                                outbam1.write(shortread)
                            else: 
                                problem_with_read = True
                 
                            shortread.query_qualities = qual

                        except Exception as e:
                            problem_with_read = True
                            pass
                         
            outbam1.close()
            outbam2.close() 

            # sort hap1 and hap2
            sortBam(hap1_bamfn, hap1_bamsortfn + '.bam', tmpbams_path)
            sortBam(hap2_bamfn, hap2_bamsortfn + '.bam', tmpbams_path)
            
            # see modify_hap and merge_haps functions for details       
            hap1_intersnpbamsortfn, hap2_intersnpbamsortfn = modify_hap(bamsortfn, hap1_bamsortfn, hap2_bamsortfn)
            hap1_finalbamsortfn, hap2_finalbamsortfn = merge_haps(bamsortfn, hap1_bamsortfn, hap2_bamsortfn, hap1_intersnpbamsortfn, hap2_intersnpbamsortfn)

            
            # remove intermediate files
            os.remove(hap1_intersnpbamsortfn + '.bam')
            os.remove(hap2_intersnpbamsortfn + '.bam')
            os.remove(hap1_intersnpbamsortfn + '.bam.bai')
            os.remove(hap2_intersnpbamsortfn + '.bam.bai')
            #os.remove(hap1_bamfn)
            #os.remove(hap2_bamfn)
            os.remove(hap1_bamsortfn + '.bam')
            os.remove(hap2_bamsortfn + '.bam')
            os.remove(hap1_bamsortfn + '.bam.bai')
            os.remove(hap2_bamsortfn + '.bam.bai')
     
        else:    
            subsample(bamsortfn, hap1_bamfn, str(0.5))
            sortBam(hap1_bamfn, hap1_finalbamsortfn + '.bam', tmpbams_path)
         
            bamDiff(bamsortfn, hap1_finalbamsortfn + '.bam', tmpbams_path)
            sortBam("/".join([tmpbams_path, 'diff_only1_' + os.path.basename(bamsortfn)]), hap2_finalbamsortfn + '.bam', tmpbams_path)
            
            # remove intermediate files
            os.remove(hap1_bamfn)
            os.remove("/".join([tmpbams_path, 'diff_only1_' +  os.path.basename(bamsortfn)]))

     
    return hap1_finalbamsortfn, hap2_finalbamsortfn

def readBamStrand(bamsortfn, strand):
    """
    Subfunction for rePair1/rePair2 functions.
    Splits ROI bam into read1 and read2 as well as forward(pos) and reverse(neg) strands.
    Fetches all reads in each file (read1_pos, read1_neg, read2_pos, read2_neg).
    """
    read1fn = sub('.bam$', '.read1_' + strand + '.bam', bamsortfn)
    read2fn = sub('.bam$', '.read2_' + strand + '.bam', bamsortfn)

    if not os.path.isfile(read1fn) or not os.path.isfile(read2fn):
        splitPairAndStrands(bamsortfn)
    
    pysam.index(read1fn)
    pysam.index(read2fn)
    
    splt1 = pysam.Samfile(read1fn, 'rb')
    splt2 = pysam.Samfile(read2fn, 'rb')
    
    itrA = splt1.fetch(until_eof=True)
    itrB = splt2.fetch(until_eof=True)
    
    return read1fn, read2fn, itrA, itrB, splt1, splt2

def defineSearchSpace(readX, strand, direction):
    """
    Subfunction for rePair1/rePair2 functions.
    Defines search space to look for potential read pairs. 
    Dependent on strand specificity as well as direction of repairing.
    """
    if (strand == 'neg' and direction == 'back') or (strand == 'pos' and direction == 'forw'):
        insert_size = readX.tlen - readX.qlen
        minpos = readX.pos + 75 + insert_size
        maxpos = readX.pos + 150 + insert_size
    
    elif (strand == 'pos' and direction == 'back') or (strand == 'neg' and direction == 'forw'):
        insert_size = abs(readX.tlen) - readX.qlen
        maxpos = readX.pos - 75 - insert_size
        minpos = readX.pos - 150 - insert_size
    return insert_size, minpos, maxpos

def generateReadPairs(tmpA, tmpB, strand, direction):
    """
    Subfunction for rePair1/rePair2 functions.
    Generates new pairs using existing reads.
    New read ID is a universally unique identifier (uuid). 
    Dependent on strand specificity as well as direction of repairing.
    """
    tlenFR = tmpB.pos - tmpA.pos + tmpB.qlen
    tlenRF = tmpA.pos - tmpB.pos + tmpA.qlen
    
    tmpqname = str(uuid4())
    
    tmpA.pnext = tmpB.pos
    tmpA.qname = tmpqname
    tmpB.pnext = tmpA.pos
    tmpB.qname = tmpqname
    
    if (strand == 'neg' and direction == 'back') or (strand == 'pos' and direction == 'forw'):
        tmpA.tlen = tlenFR
        tmpB.tlen = -tlenFR
    elif (strand == 'pos' and direction == 'back') or (strand == 'neg' and direction == 'forw'):
        tmpA.tlen = -tlenRF
        tmpB.tlen = tlenRF
    
    return tmpA, tmpB

def rePair1(bamsortfn):
    """
    Subfunction for repairReads function.
    Create new pairs from read1 -> read2.
    """
    # Throws an error if bamsortfn is not found
    if not os.path.isfile(bamsortfn):
        raise ValueError('Could not find file bamsortfn')
    bamrepairedfn = sub('.bam$', ".re_paired.bam", bamsortfn)
    errorbamrepairedfn = sub('.bam$', ".ERROR.re_paired.bam", bamsortfn)
    bamrepairedsortfn = sub('.bam$', ".re_paired.sorted.bam", bamsortfn)
 
    inbam = pysam.Samfile(bamsortfn, 'rb')
    outbam = pysam.Samfile(bamrepairedfn, 'wb', template=inbam)
    errorbam = pysam.Samfile(errorbamrepairedfn, 'wb', template=inbam)

    writtencount = 0
    strands = ['pos', 'neg']

    for strand in strands:
        # Takes bamsortfn and splits it based on Read 1/2 and Strand
        read1fn, read2fn, itrA, itrB, splt1, splt2 = readBamStrand(bamsortfn, strand)
        counter = 0 
        
        while (True):
            try:
                counter += 1
                direction='forw'
                readRef = itrA.next() 
                chromosome = readRef.reference_name

                    # Defines the search space for the other read in the opposite splt
                insert_size, minpos, maxpos = defineSearchSpace(readRef, strand, direction)
                if (insert_size > 0 and insert_size < 1000):
                    itrTarget = splt2.fetch(chromosome, minpos, maxpos)
                        
                    listTarget = []
                    itrs_list = list(itrTarget)
                        
                        
                    if len(itrs_list) <= 5: # Takes all target reads
                        listTarget = itrs_list
                    
                    elif len(itrs_list) > 5: # Takes a random sample of target reads
                        listTarget = [i for i in random.sample(itrs_list, 5)]
                
                # Loops through all target reads
                    for i in range(len(listTarget)):
                        readTarget = listTarget[i]
                    # If the read IDs dont match, create a new read-pair by altering the description of the read and output
                        if readRef.qname != readTarget.qname:
                            tmpA, tmpB = generateReadPairs(readRef, readTarget, strand, direction)
                            if counter % 2 != 0:
                                outbam.write(tmpA)
                                outbam.write(tmpB)
                elif (insert_size < 0 or insert_size > 1000): 
                    errorbam.write(readRef)

            except StopIteration:
                break
                    
        splt1.close()
        splt2.close()
        os.remove(read1fn)
        os.remove(read2fn)
        os.remove(read1fn + '.bai')
        os.remove(read2fn + '.bai')

    inbam.close()
    outbam.close()
    errorbam.close()

    bamrepairedsortfn = sub('sorted.re_paired', 're_paired', bamrepairedsortfn)
    sortBam(bamrepairedfn, bamrepairedsortfn, tmpbams_path)
    os.remove(bamrepairedfn)

    return bamrepairedsortfn


def rePair2(bamsortfn):
    """
    Subfunction for repairReads function.
    Create new pairs from read2 -> read1.
    """
    # Throws an error if bamsortfn is not found
    if not os.path.isfile(bamsortfn):
        raise ValueError('Could not find file bamsortfn')

    bamrepaired2fn = sub('.bam$', ".re_paired2.bam", bamsortfn)
    errorbamrepaired2fn = sub('.bam$', ".ERROR.re_paired2.bam", bamsortfn)
    bamrepaired2sortfn = sub('.bam$', ".re_paired2.sorted.bam", bamsortfn)
 
    inbam = pysam.Samfile(bamsortfn, 'rb')
    outbam2 = pysam.Samfile(bamrepaired2fn, 'wb', template=inbam)
    errorbam2 = pysam.Samfile(errorbamrepaired2fn, 'wb', template=inbam)

    writtencount = 0
    strands = ['pos', 'neg']

    for strand in strands:
        # Takes bamsortfn and splits it based on Read 1/2 and Strand
        read1fn, read2fn, itrA, itrB, splt1, splt2 = readBamStrand(bamsortfn, strand)
        counter = 0 
        
        while (True):
            try:
                counter += 1
                direction='back'
                readRef = itrB.next()
                chromosome = readRef.reference_name

                    # Defines the search space for the other read in the opposite splt
                insert_size, minpos, maxpos = defineSearchSpace(readRef, strand, direction)
                if (insert_size > 0 and insert_size < 1000):
                    itrTarget = splt1.fetch(chromosome, minpos, maxpos)
                    listTarget = []
                    itrs_list = list(itrTarget)
    
                    if len(itrs_list) <= 5: # Takes all target reads
                        listTarget = itrs_list
                    elif len(itrs_list) > 5: # Takes a random sample of target reads
                        listTarget = [i for i in random.sample(itrs_list, 5)]
                
                # Loops through all target reads
                    for i in range(len(listTarget)):
                        readTarget = listTarget[i]

                    # If the read IDs dont match, create a new read-pair by altering the description of the read and output  
                        if readRef.qname != readTarget.qname:
                            tmpA, tmpB = generateReadPairs(readRef, readTarget, strand, direction)
                            
                            if counter % 2 == 0:
                                outbam2.write(tmpA)
                                outbam2.write(tmpB)

                elif (insert_size < 0 or insert_size > 1000): 
                    errorbam2.write(readRef)
            except StopIteration:
                break
            
        splt1.close()
        splt2.close()
        os.remove(read1fn)
        os.remove(read2fn)
        os.remove(read1fn + '.bai')
        os.remove(read2fn + '.bai')

    inbam.close()
    outbam2.close()
    errorbam2.close()

    bamrepaired2sortfn = sub('sorted.re_paired', 're_paired', bamrepaired2sortfn)
    sortBam(bamrepaired2fn, bamrepaired2sortfn, tmpbams_path)
    os.remove(bamrepaired2fn)

    return bamrepaired2sortfn

def merge_pairs(bamsortfn):
    """
    Subfunction for repairReads function.
    Perform previous repairing subfunctions with hap1 or hap2 bam from split_hap function.
    Merge repaired reads from rePair1 and rePair2.
    Mark and remove duplicates of the same pairs generated by rePair1 and rePair2.    
    """
    
    bamrepairedfinalfn = sub('.sorted.bam$', ".re_paired_final.bam", bamsortfn)
    bamrepairedfinalsortfn = sub('.sorted.bam$', ".re_paired_final.sorted.bam", bamsortfn)
    bamrepairedfinalmarkedfn = sub('.sorted.bam$', ".re_paired_final.marked.bam", bamsortfn)
    bamrepairedfinalsortmarkedfn = sub('.sorted.bam$', ".re_paired_final.marked.sorted.bam", bamsortfn)

    bamrepairedsortfn = rePair1(bamsortfn)
    bamrepaired2sortfn = rePair2(bamsortfn)

    merge_bams(bamrepairedsortfn, bamrepaired2sortfn, bamrepairedfinalfn)
    sortBam(bamrepairedfinalfn, bamrepairedfinalsortfn, tmpbams_path)
    
    print(" ___ removing repaired duplicates ___")
    removeDupSambamba(bamrepairedfinalsortfn, tmpbams_path)
    sortBam(bamrepairedfinalmarkedfn, bamrepairedfinalsortmarkedfn, tmpbams_path)

    # remove intermediate files
    os.remove(bamrepairedsortfn)
    os.remove(bamrepaired2sortfn)
    os.remove(bamrepairedfinalfn)
    os.remove(bamrepairedfinalsortfn)
    os.remove(bamrepairedfinalmarkedfn)
    os.remove(bamrepairedsortfn + '.bai')
    os.remove(bamrepaired2sortfn + '.bai')
    os.remove(bamrepairedfinalfn + '.bai')
    os.remove(bamrepairedfinalsortfn + '.bai')
    os.remove(bamrepairedfinalmarkedfn + '.bai')
    
    
    return bamrepairedfinalsortmarkedfn

def repairReads(bamsortfn):
    """
    Uses split_hap output to generate hap1 and hap2 repaired bams.
    """

    hap1_finalbamsortfn = sub('.sorted.bam$', ".hap1_final.sorted.bam", bamsortfn)
    hap2_finalbamsortfn = sub('.sorted.bam$', ".hap2_final.sorted.bam", bamsortfn)
    
    print(" ___ re-pairing hap1 bam reads ___") 
    hap1_bamrepairedfinalsortmarkedfn = merge_pairs(hap1_finalbamsortfn)
    
    print(" ___ re-pairing hap2 bam reads ___") 
    hap2_bamrepairedfinalsortmarkedfn = merge_pairs(hap2_finalbamsortfn)
    
    #os.remove(hap1_finalbamsortfn)
    #os.remove(hap2_finalbamsortfn) 
    #os.remove(hap1_finalbamsortfn + '.bai')
    #os.remove(hap2_finalbamsortfn + '.bai')

    return hap1_bamrepairedfinalsortmarkedfn, hap2_bamrepairedfinalsortmarkedfn


def split_bam_by_chr(chr):
    inbam = params.GetInputBam()
    spltbams_path = "/".join([res_path, 'splitbams'])

    try:
        if not terminating.is_set():
            logger.debug("___ spliting bam by chromosome ___")
            splitBamByChr(inbam, spltbams_path, str(chr))
            sortByName("/".join([spltbams_path, str(chr) + ".bam"]),
                       "/".join([spltbams_path, str(chr) + ".byname.bam"]))

    except KeyboardInterrupt:
        logger.error('Exception Crtl+C pressed in the child process  in split_bam_by_chr')
        terminating.set()
        return False

    except Exception as e:
        logger.exception("Exception in split_bam_by_chr %s", e)
        terminating.set()
        return False

    return


def implement_cnv(chromosome_event):
    """
    Performs all previous subfunctions.
    Calculates coverage ratio based on repaired read count compared to ROI read count. 
    Downsamples to the requested allelic ratio specificed by user in cnv.bed.
    """
    chr, event = chromosome_event.split("_")
    logger.debug("___ Bamgineer main engine started ___")
    success = True

    try:
        if not terminating.is_set():
            bamfn, sortbyname, sortbyCoord, bedfn = init_file_names(chr, tmpbams_path, haplotype_path, event)
            exonsinroibed = "/".join([haplotype_path, chr + "_exons_in_roi" + event + ".bed"])
            bamsortfn = sub('.bam$', '.sorted.bam', bamfn)
            hap1_finalbamsortfn = sub('.sorted.bam$', ".hap1_final.sorted.bam", bamsortfn)
            hap2_finalbamsortfn = sub('.sorted.bam$', ".hap2_final.sorted.bam", bamsortfn)

            if ((not os.path.isfile(bedfn)) and (os.path.isfile(exonsinroibed))):
                bedfn = exonsinroibed

            if os.path.isfile(bedfn):
                if os.path.samefile(bedfn, exonsinroibed):
                    fn = list(csv.reader(open(bedfn, 'rb'), delimiter='\t'))
                    copy_number = int(fn[0][4])
                    hap_type = str(fn[0][3])
        
                else:
                    fn = list(csv.reader(open(bedfn, 'rb'), delimiter='\t'))
                    copy_number = int(fn[0][7])
                    hap_type = str(fn[0][6])

                if event.startswith('amp') or event.startswith('gain') or event.startswith('loh') or event.startswith('loss'):
                    
                    bamrepairedsortfn = sub('.sorted.bam$', ".re_paired.sorted.bam", bamsortfn)
                    hap1_bamrepairedfinalsortmarkedfn = sub('.sorted.bam$', ".hap1_final.re_paired_final.marked.sorted.bam", bamsortfn)
                    hap2_bamrepairedfinalsortmarkedfn = sub('.sorted.bam$', ".hap2_final.re_paired_final.marked.sorted.bam", bamsortfn)
                    hap1_final = sub('.sorted.bam$', ".hap1_FINAL.sorted.bam", bamsortfn)
                    hap2_final = sub('.sorted.bam$', ".hap2_FINAL.sorted.bam", bamsortfn)
                    hap1_finalmarked = sub('.sorted.bam$', ".hap1.marked.bam", bamsortfn)
                    hap2_finalmarked = sub('.sorted.bam$', ".hap2.marked.bam", bamsortfn)
                    hap1_finalbamsortfn = sub('.sorted.bam$', ".hap1_final.sorted.bam", bamsortfn)
                    hap2_finalbamsortfn = sub('.sorted.bam$', ".hap2_final.sorted.bam", bamsortfn)
                    HAP1_FINAL = "/".join([tmpbams_path, str(chr).upper()+ str(event).upper()+ '_HAP1.bam'])
                    HAP2_FINAL = "/".join([tmpbams_path, str(chr).upper()+ str(event).upper()+ '_HAP2.bam'])
                    HAP12_FINAL = "/".join([finalbams_path, str(chr).upper()+ str(event).upper()+ '_HAP12.bam'])

                    if os.path.isfile(bamsortfn):

                        split_hap(bamsortfn, chr, event)
                        repairReads(bamsortfn)
                        merge_bams(hap1_finalbamsortfn, hap1_bamrepairedfinalsortmarkedfn, hap1_final)
                        merge_bams(hap2_finalbamsortfn, hap2_bamrepairedfinalsortmarkedfn, hap2_final)
                        
                        print("___ removing hap1 merged normal duplicates ___")
                        removeDupSambamba(hap1_final, tmpbams_path)
                        print("___ removing hap2 merged normal duplicates ___")
                        removeDupSambamba(hap2_final, tmpbams_path)
                        
                        #coverageratio = (float(countReads(hap1_bamrepairedfinalsortmarkedfn))+ float(countReads(hap2_bamrepairedfinalsortmarkedfn))) / float(countReads(bamsortfn))
                        coverageratio = (float(countReads(hap1_finalmarked))+ float(countReads(hap2_finalmarked))) / float(countReads(bamsortfn))
                        logger.debug("+++ coverage ratio for: " + ntpath.basename(bamsortfn) + ": " + str(coverageratio))

                        alleleA = hap_type.count('A')
                        alleleB = hap_type.count('B')
                        samplerate1 = float((alleleA/coverageratio)+1) # 1 is random seed
                        samplerate2 = float((alleleB/coverageratio)+1) # 1 is random seed

                        if coverageratio < copy_number/2:
                            logger.error('not enough reads or repairing search space is too small for ' + ntpath.basename(bamsortfn))
                            success = False
                            return

                        elif alleleA + alleleB != copy_number:
                            logger.error('allelic ratio adds up incorrectly (correct: AAB is CN=3)')
                            success = False
                            return
			
			elif alleleA > coverageratio or alleleB > coverageratio:
                            logger.error('requested individual allelic ratio is greater than available repaired reads')
			    success = False
                            return

                        elif alleleB == 0:
			    if alleleA > coverageratio:
                                logger.error('requested individual allelic ratio is greater than available repaired reads')
			        success = False
                                return
                            coverageratio = float(countReads(hap1_finalmarked)) / float(countReads(hap1_finalbamsortfn))
                            samplerate1 = float((alleleA/coverageratio)+1) # 1 is random seed
                            subsample(hap1_finalmarked, HAP1_FINAL, str(samplerate1))
                            sortBam(HAP1_FINAL, HAP12_FINAL, tmpbams_path)
                            success = True
                        
                        elif alleleA == 0:
			    if alleleB > coverageratio:
                                logger.error('requested individual allelic ratio is greater than available repaired reads')
			        success = False
                                return
                            coverageratio = float(countReads(hap2_finalmarked)) / float(countReads(hap2_finalbamsortfn))
                            samplerate2 = float((alleleB/coverageratio)+1) # 1 is random seed
                            subsample(hap2_finalmarked, HAP2_FINAL, str(samplerate2))
                            sortBam(HAP2_FINAL, HAP12_FINAL, tmpbams_path)
                            success = True

                        else:
                            subsample(hap1_finalmarked, HAP1_FINAL, str(samplerate1))
                            subsample(hap2_finalmarked, HAP2_FINAL, str(samplerate2))
                            merge_bams(HAP1_FINAL, HAP2_FINAL, HAP12_FINAL)
                            success = True

            
            else:
                logger.debug(bedfn + ' does not exist!')
                success = False

    except KeyboardInterrupt:
        logger.error('Exception Crtl+C pressed in the child process  in find_roi_bam for chr ' + chr + event)
        terminating.set()
        success = False
        return

    except Exception as e:
        logger.exception("Exception in find_roi_bam %s", e)
        terminating.set()
        success = False
        return

    if success:
        logger.debug("implement_cnv complete successfully for " + chr + event)

    return


def removeReadsOverlappingHetRegion(inbamfn, bedfn, outbamfn, path):
    print " ___ removing reads overlapping heterozygous region ___"
    inbamsorted = sub('.bam$', '.sorted', inbamfn)
    pysam.sort(inbamfn, inbamsorted)
    pysam.index(inbamsorted + '.bam')

    alignmentfile = pysam.AlignmentFile(inbamsorted + '.bam', "rb")
    outbam = pysam.Samfile(outbamfn, 'wb', template=alignmentfile)

    bedfile = open(bedfn, 'r')

    for bedline in bedfile:
        c = bedline.strip().split()

        if len(c) == 3:
            chr2 = c[0]
            chr = c[0].strip("chr")
            start = int(c[1])
            end = int(c[2])
        else:
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
    
    outbamsorted = sub('.bam$', '.sorted', outbamfn)
    pysam.sort(outbamfn, outbamsorted)
    bamDiff(inbamsorted + '.bam', outbamsorted + '.bam', path)
    outbam.close()


def run_pipeline(results_path):
    """
    Pools parallel processes using multiprocessing.
    Final merge of non-roi bams and bamgineered bams.
    """
    print(results_path)
    global haplotype_path, cancer_dir_path, tmpbams_path, finalbams_path, log_path, logfile, terminating, logger, logQueue, res_path
    res_path = results_path
    haplotype_path, cancer_dir_path, tmpbams_path, finalbams_path, log_path, logfile = handle.GetProjectPaths(results_path)
    terminating, logger, logQueue = handle.GetLoggings(logfile)

    chr_list = ['chr' + str(x) for x in range(1, 23)]
    chr_list.extend(['chrX', 'chrY'])

    t0 = time.time()
    outbamfn = params.GetOutputFileName()

    cnv_list = glob.glob("/".join([params.GetCNVDir(), '*.*']))
    chromosome_event = create_chr_event_list(cnv_list, chr_list)

    logger.debug('pipeline started!')

    phase_path = '/'.join([results_path, 'phasedvcfdir'])
    
    if not os.path.exists('/'.join([results_path, 'phasedvcfdir'])):
        os.makedirs(phase_path)

    initialize0(phase_path, cancer_dir_path)

    for cnv_path in cnv_list:
        initialize_pipeline(phase_path, haplotype_path, cnv_path)

    pool1 = multiprocessing.Pool(processes=12, initializer=initPool, initargs=[logQueue, logger.getEffectiveLevel(), terminating])
    try:

        if not params.GetSplitBamsPath():

            if not os.path.exists("/".join([res_path, 'splitbams'])):
                os.makedirs("/".join([res_path, 'splitbams']))
                params.SetSplitBamsPath("/".join([res_path, 'splitbams']))

            result0 = pool1.map_async(split_bam_by_chr, chromosome_event).get(9999999)

        result1 = pool1.map_async(find_roi_bam, chromosome_event).get(9999999)
        result2 = pool1.map_async(implement_cnv, chromosome_event).get(9999999)
        result3 = pool1.map_async(find_non_roi_bam, chr_list).get(9999999)
        pool1.close()

    except KeyboardInterrupt:
        logger.debug('You cancelled the program!')
        pool1.terminate()

    except Exception as e:
        logger.exception("Exception in main %s", e)
        pool1.terminate()

    finally:
        pool1.join()

    time.sleep(.1)
    merge_final(outbamfn, finalbams_path)
    #mergeSortBamFiles(outbamfn, finalbams_path)
    t1 = time.time()
    #shutil.rmtree(tmpbams_path)
    logger.debug(' ***** pipeline finished in ' + str(round((t1 - t0) / 60.0, 1)) + ' minutes ***** ')
    logging.shutdown()
