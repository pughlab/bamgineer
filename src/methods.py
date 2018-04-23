import csv
import glob
import logging
import time
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
            thinVCF(hap1vcf, hap1vcffiltered)
            thinVCF(hap2vcf, hap2vcffiltered)
            convertvcftobed(hap1vcffiltered + ".recode.vcf", hap1vcffilteredtobed)
            convertvcftobed(hap2vcffiltered + ".recode.vcf", hap2vcffilteredtobed)

            generatePhasedBed(hap1vcffilteredtobed, hap2vcffilteredtobed, phased_bed)

    except:

        raise

    return


def initialize_pipeline(phase_path, haplotype_path, cnv_path):
    exons_path = bamhelp.GetExons()

    event, extension = os.path.splitext(os.path.basename(cnv_path))

    phased_bed = "/".join([phase_path, "PHASED.BED"])
    bedtools_path = bamhelp.GetBedtoolsPath()

    try:
        logger.debug(' --- Initializing input files  --- ')
        exonsinroibed = "/".join([haplotype_path, "exons_in_roi" + str(event) + ".bed"])

        nonhetbed = "/".join([haplotype_path, "non_het" + str(event) + ".bed"])
        hetbed = "/".join([haplotype_path, "het" + str(event) + ".bed"])
        hetsnpbed = "/".join([haplotype_path, "het_snp" + str(event) + ".bed"])

        tmp = "/".join([haplotype_path, str(event) + "_tmp.bed"])
        command = " ".join([bedtools_path, "intersect -a", exons_path, "-b", cnv_path, "-wa -wb > ", tmp])
        runCommand(command)

        filterColumns(tmp, exonsinroibed, [0, 1, 2])

        splitBed(exonsinroibed, '_exons_in_roi' + str(event))
        command = " ".join([bedtools_path, "intersect -a", phased_bed, "-b", exonsinroibed, "-wa -wb >", tmp])
        runCommand(command)

        filterColumns(tmp, hetsnpbed, [i for i in range(1, 6)])

        splitBed(hetsnpbed, '_het_snp' + str(event))
        os.remove(tmp)
    except:
        logger.exception("Initialization error !")
        raise
    logger.debug("--- initialization complete ---")
    return


def init_file_names(chr, tmpbams_path, haplotypedir, event=''):
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


def re_pair_reads(bamsortfn, copy_number):
    if os.path.isfile(bamsortfn):

        bamrepairedfn = sub('.bam$', ".re_paired.bam", bamsortfn)
        bamrepairedsortfn = sub('.bam$', ".re_paired.sorted.bam", bamsortfn)

        inbam = pysam.Samfile(bamsortfn, 'rb')
        outbam = pysam.Samfile(bamrepairedfn, 'wb', template=inbam)

        writtencount = 0
        strands = ['pos', 'neg']

        for strand in strands:
            read1fn = sub('.bam$', '.read1_' + strand + '.bam', bamsortfn)
            read2fn = sub('.bam$', '.read2_' + strand + '.bam', bamsortfn)

            if not os.path.isfile(read1fn) or not os.path.isfile(read2fn):
                splitPairAndStrands(bamsortfn)

            splt1 = pysam.Samfile(read1fn, 'rb')
            splt2 = pysam.Samfile(read2fn, 'rb')

            itrA = splt1.fetch(until_eof=True)
            itrB = splt2.fetch(until_eof=True)

            if (params.GetctDNA()):
                sigma = 40
                coff = 2
                block_size = int(copy_number)
            else:
                sigma = 85
                coff = 5
                block_size = int(copy_number) * 4

            writtenreads = []

            while (True):

                try:
                    lista = []
                    listb = []
                    readsinblock = []
                    poslist = []

                    while (True):

                        readA = itrA.next()
                        lista.append(readA)
                        readsinblock.append(readA.qname)
                        poslist.append(readA.pos)

                        if (len(lista) == block_size):
                            break

                    minpos = min(poslist) - 1000
                    maxpos = max(poslist) + 1000

                    while (True):

                        readB = itrB.next()

                        if readB.qname in readsinblock or (readB.pos < maxpos or readB.pos > minpos):
                            listb.append(readB)

                        if len(listb) == block_size:
                            break

                    for i in range(0, block_size):
                        readA = lista[i]
                        tmpA = readA
                        for j in range(0, block_size):

                            readB = listb[j]
                            tmpB = readB

                            tlenFR = tmpB.pos - tmpA.pos + tmpB.qlen
                            tlenRF = tmpA.pos - tmpB.pos + tmpA.qlen

                            if readA.qname != readB.qname:

                                tmpqname = str(uuid4())

                                if strand == 'pos':
                                    if tlenFR >= abs(tmpB.tlen) - coff * sigma and tlenFR < abs(tmpB.tlen) + coff * sigma:
                                        tmpA.tlen = tlenFR
                                        tmpB.tlen = -tlenFR
                                        tmpA.pnext = tmpB.pos
                                        tmpB.pnext = tmpA.pos
                                        tmpA.qname = tmpqname
                                        tmpB.qname = tmpqname
                                        outbam.write(tmpA)
                                        outbam.write(tmpB)
                                        writtenreads.append(tmpB.qname)

                                elif strand == 'neg':
                                    if tlenRF >= abs(tmpB.tlen) - coff * sigma and tlenRF < abs(tmpB.tlen) + coff * sigma:
                                        tmpA.tlen = -tlenRF
                                        tmpB.tlen = tlenRF
                                        tmpA.pnext = tmpB.pos
                                        tmpB.pnext = tmpA.pos
                                        tmpA.qname = tmpqname
                                        tmpB.qname = tmpqname
                                        outbam.write(tmpA)
                                        outbam.write(tmpB)
                                        writtenreads.append(tmpB.qname)

                except StopIteration:
                    break

            os.remove(read1fn)
            os.remove(read2fn)

        splt1.close()
        splt2.close()

        inbam.close()
        outbam.close()

        bamrepairedsortfn = sub('sorted.re_paired', 're_paired', bamrepairedsortfn)
        sortBam(bamrepairedfn, bamrepairedsortfn, tmpbams_path)
        os.remove(bamrepairedfn)

    return


def mutate_reads(bamsortfn, chr, event=''):
    fn, sortbyname, sortbyCoord, bedfn = init_file_names(chr, tmpbams_path, haplotype_path, event)
    cmd = " ".join(["sort -u", bedfn, "-o", bedfn]);
    runCommand(cmd)
    hetbamfn = sub('.sorted.bam$', ".mutated_het.bam", bamsortfn)
    hetbamfnsorted = sub('.sorted.bam$', ".mutated_het.sorted", bamsortfn)
    allreadsfn = sub('.sorted.bam$', ".all.reads.bam", bamsortfn)
    allreadssortfn = sub('.sorted.bam$', ".all.reads.sorted", bamsortfn)
    mergedsortfn = sub('.sorted.bam$', ".mutated_merged.sorted.bam", bamsortfn)
    try:
        if not terminating.is_set():

            if (os.path.isfile(bamsortfn) and os.path.isfile(bedfn)):

                samfile = pysam.Samfile(bamsortfn, "rb")
                alignmentfile = pysam.AlignmentFile(bamsortfn, "rb")
                outbam = pysam.Samfile(hetbamfn, 'wb', template=samfile)
                allreads = pysam.Samfile(allreadsfn, 'wb', template=samfile)

                bedfile = open(bedfn, 'r')
                covpath = "/".join([haplotype_path, "written_coverage_het.txt"])
                snpratiopath = "/".join([haplotype_path, "het_snp_ratio.txt"])

                num_reads_written = 0
                num_total_reads = 0

                for bedline in bedfile:
                    c = bedline.strip().split()

                    if len(c) == 7:
                        chr2 = c[0]
                        chr = c[0].strip("chr")
                        start = int(c[1])
                        end = int(c[2])
                        refbase = str(c[3])
                        altbase = str(c[4])
                        haplotype = str(c[5])
                        copy_number = int(c[6])
                    else:
                        continue

                    readmappings = alignmentfile.fetch(chr2, start, end)

                    # sex chromosome
                    if params.GetXY() and (chr == 'chrX' or chr == 'chrY'):
                        haplotype = 'hap1'
                        print('sex chromosome ' + str(chr))

                    for shortread in readmappings:

                        allreads.write(shortread)
                        num_total_reads += 1
                        problem_with_read = False

                        try:
                            index = shortread.get_reference_positions(full_length=True).index(start)
                            tmpread = shortread.query_sequence
                            qual = shortread.query_qualities
                            mutated_hap1 = tmpread[:index] + altbase + tmpread[index + 1:]
                            mutated_hap2 = tmpread[:index] + refbase + tmpread[index + 1:]
                            if haplotype == "hap1":
                                shortread.query_sequence = mutated_hap1

                            elif haplotype == "hap2":
                                shortread.query_sequence = mutated_hap2

                            shortread.query_qualities = qual

                        except Exception as e:
                            problem_with_read = True
                            pass

                        if not problem_with_read:
                            outbam.write(shortread)
                            num_reads_written += 1

                outbam.close()
                allreads.close()

                sortBam(hetbamfn, hetbamfnsorted + '.bam', tmpbams_path)
                sortBam(allreadsfn, allreadssortfn + '.bam', tmpbams_path)

                os.remove(hetbamfn)
                os.remove(allreadsfn)

                # ratio of het reads to nonhet reads, we need to adjust the coverage
                ratio = float(num_reads_written) / float(num_total_reads)
                bamsortfnsampled = sub('.sorted.bam$', ".sampled.nh.bam", bamsortfn)
                subsample(bamsortfn, bamsortfnsampled, str(ratio))
                bamDiff(bamsortfnsampled, allreadssortfn + '.bam', tmpbams_path)

                if "/".join([tmpbams_path, 'diff_only1_' + os.path.basename(bamsortfnsampled)]):
                    merge_bams("/".join([tmpbams_path, 'diff_only1_' + os.path.basename(bamsortfnsampled)]),
                               hetbamfnsorted + '.bam', mergedsortfn)
                    os.remove("/".join([tmpbams_path, 'diff_only1_' + os.path.basename(bamsortfnsampled)]))

                os.remove(bamsortfnsampled)
                os.remove(allreadssortfn + '.bam')
                os.remove(allreadssortfn + '.bam.bai')

                os.remove(hetbamfnsorted + '.bam')
                os.remove(hetbamfnsorted + '.bam.bai')

    except KeyboardInterrupt:
        logger.error('Exception Crtl+C pressed in the child process  in mutaute_reads')
        terminating.set()
        return
    except Exception as e:
        logger.exception("Exception in mutate_reads %s", e)
        terminating.set()
        return
    return


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


# cn change is 1 for CN=1,2,...,8
def calculate_sample_rate(inbam, outbam, cnchange, purity):
    logger.debug("___ adjusting sample rate ___")


def implement_cnv(chromosome_event):
    chr, event = chromosome_event.split("_")

    logger.debug("___ Bamgineer main engine started ___")
    success = True
    try:
        if not terminating.is_set():
            bamfn, sortbyname, sortbyCoord, bedfn = init_file_names(chr, tmpbams_path, haplotype_path, event)
            bamsortfn = sub('.bam$', '.sorted.bam', bamfn)

            if os.path.isfile(bedfn):
                fn = list(csv.reader(open(bedfn, 'rb'), delimiter='\t'))
                copy_number = int(fn[0][6])

                if not params.GetXY() or (chr != 'chrX' and chr != 'chrY'):

                    if copy_number == 2:
                        event = 'loh'
                    elif copy_number == 3:
                        event = 'gain'
                    elif copy_number > 3:
                        event = 'amp'

                else:

                    logger.debug("*** handling single sex chromosome for: " + ntpath.basename(bamsortfn))
                    if copy_number == 1:
                        event = 'loh'
                    elif copy_number == 2:
                        event = 'gain'
                    elif copy_number > 2:
                        event = 'amp'

                if event.startswith('amp') or event.startswith('gain'):

                    bamrepairedsortfn = sub('.sorted.bam$', ".re_paired.sorted.bam", bamsortfn)
                    mergedsortfn = sub('.sorted.bam$', ".mutated_merged.sorted.bam", bamrepairedsortfn)
                    GAIN_FINAL = "/".join([finalbams_path, str(chr).upper() + '_GAIN.bam'])

                    if os.path.isfile(bamsortfn):

                        re_pair_reads(bamsortfn, copy_number)
                        mutate_reads(bamrepairedsortfn, chr, event)
                        coverageratio = float(countReads(mergedsortfn)) / float(countReads(bamsortfn))
                        logger.debug(
                            "+++ coverage ratio for: " + ntpath.basename(bamsortfn) + ": " + str(coverageratio))

                        if coverageratio < copy_number - 2:
                            logger.error('not enough reads for ' + ntpath.basename(bamsortfn))
                            return
                        else:
                            samplerate = float(copy_number - 2) / coverageratio
                            subsample(mergedsortfn, GAIN_FINAL, str(samplerate))

                elif event == 'loss':

                    inbam_deletion = "/".join([finalbams_path, str(chr).upper() + '_LOSS.bam'])

                    if os.path.isfile(bamsortfn):

                        mutate_reads(bamsortfn, chr, 'loss')
                        mergedsortfn = sub('.sorted.bam$', ".mutated_merged.sorted.bam", bamsortfn)
                        mergedsortsampledfn = sub('.sorted.bam$', ".mutated_merged.sampled.sorted.bam", bamsortfn)

                        ratio_kept = float(countReads(mergedsortfn)) / float(countReads(bamsortfn))
                        samplerate = round(0.5 / ratio_kept, 2)
                        LOSS_FINAL = "/".join([finalbams_path, str(chr).upper() + '_LOSS.bam'])
                        logger.debug("ratios kept for:" + ntpath.basename(bamsortfn) + ": " + str(ratio_kept))
                        subsample(mergedsortfn, mergedsortsampledfn, str(samplerate))
                        bamDiff(sortbyCoord, mergedsortsampledfn, tmpbams_path)
                        os.rename("/".join([tmpbams_path, 'diff_only1_' + chr + '.bam']), LOSS_FINAL)

                    elif (not os.path.isfile(inbam_deletion) and os.path.isfile(
                            sortbyCoord)):  # if it exists from previous runs

                        os.symlink(sortbyCoord, inbam_deletion)

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
    print "___ removing reads overlapping heterozygous region ___"
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
    print(results_path)
    global haplotype_path, cancer_dir_path, tmpbams_path, finalbams_path, log_path, logfile, terminating, logger, logQueue, res_path
    res_path = results_path
    haplotype_path, cancer_dir_path, tmpbams_path, finalbams_path, log_path, logfile = handle.GetProjectPaths(
        results_path)
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

    sampool1 = multiprocessing.Pool(processes=12, initializer=initPool,
                                 initargs=[logQueue, logger.getEffectiveLevel(), terminating])
    try:

        if not params.GetSplitBamsPath():

            if not os.path.exists("/".join([res_path, 'splitbams'])):
                os.makedirs("/".join([res_path, 'splitbams']))
                params.SetSplitBamsPath("/".join([res_path, 'splitbams']))

            result0 = pool1.map_async(split_bam_by_chr, chromosome_event).get(9999999)

        result1 = pool1.map_async(find_roi_bam, chromosome_event).get(9999999)
        result2 = pool1.map_async(implement_cnv, chromosome_event).get(9999999)
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
    mergeSortBamFiles(outbamfn, finalbams_path)
    t1 = time.time()
    shutil.rmtree(tmpbams_path)
    logger.debug(' ***** pipeline finished in ' + str(round((t1 - t0) / 60.0, 1)) + ' minutes ***** ')
    logging.shutdown()