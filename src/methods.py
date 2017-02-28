import pysam
import ntpath
from helpers import parameters as params
from helpers import handlers as handle
from helpers import bamgineerHelpers as bamhelp
from utils import *
import logging, sys
import random

global bases
bases = ('A','T','C','G')


def initPool(queue, level, terminating_):
    """
    This causes the logging module to be initialized with the necessary info
    in pool threads to work correctly.
    """
    logging.getLogger('').setLevel(level)
    global terminating
    terminating = terminating_

def initialize(results_path,haplotype_path,cancer_dir_path):
    
    try:
        event_list=['gain','loss']
        gaincnv = params.GetGainCNV()
        losscnv = params.GetLossCNV()
        
        logger.debug(' --- Initializing input files  --- ')
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
          
        phaseVCF(vcf_path, phasedvcf)
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
            roibed = "/".join([haplotype_path,  event + "_roi.bed"])
            exonsinroibed = "/".join([haplotype_path,   event + "_exons_in_roi.bed"])
            nonhetbed = "/".join([haplotype_path, event + "_non_het.bed"])
            hetbed = "/".join([haplotype_path, event + "_het.bed"])
            hetsnpbed = "/".join([haplotype_path,  event + "_het_snp.bed"])
            
            intersectBed( exons_path, locals()[event + 'cnv'], exonsinroibed, wa=True)
            intersectBed(phased_bed, exonsinroibed, hetsnpbed, wa=True)
            splitBed(exonsinroibed, event+'_exons_in_roi_')
            splitBed(hetsnpbed, event+'_het_snp_')
            
        
    except:  
        logger.exception("Initialization error !")
        raise
    
    logger.debug("--- initialization complete ---")    
    return 

def init_file_names(chr, event,tmpbams_path, haplotypedir):
    
    flist=[]
    splitbams = params.GetSplitBamsPath()
    roibam = "/".join([tmpbams_path ,chr + event +"_roi.bam"])
    sortbyname =  "/".join([splitbams,  chr + '.byname.bam'])
    sortbyCoord = "/".join([splitbams,  chr + '.bam'])
    hetsnp   = "/".join([haplotypedir, event+'_het_snp_' + chr + '.bed'])
    flist.extend([roibam,sortbyname,sortbyCoord,hetsnp])
    return flist

def find_roi_bam(chromosome_event):
    chr,event = chromosome_event .split("_")
    roi,sortbyname,sortbyCoord, hetsnp = init_file_names(chr, event, tmpbams_path, haplotype_path)
    exonsinroibed = "/".join([haplotype_path,   event + "_exons_in_roi_"+ chr +'.bed'])
    success = True
    try:
        if not terminating.is_set():
            roisort = sub('.bam$', '.sorted', roi)
            if(os.path.isfile(exonsinroibed)):
                 
                 cmd=" ".join(["sort -u", exonsinroibed, "-o", exonsinroibed]); runCommand(cmd)
                 extractPairedReadfromROI(sortbyname, exonsinroibed, roi)
                 removeIfEmpty(tmpbams_path,ntpath.basename(roi))
               
                 getProperPairs(roi, roi+'.tmp.bam')
                 pysam.sort(roi+'.tmp.bam',roisort )
                 pysam.index(roisort+'.bam')
                 os.remove(roi+'.tmp.bam')
              
    except (KeyboardInterrupt):
        logger.error('Exception Crtl+C pressed in the child process  in find_roi_bam for chr ' +chr + event)
        terminating.set()
        success=False
        return
    except Exception as e:   
        logger.exception("Exception in find_roi_bam %s" ,e )
        terminating.set()
        success=False
        return
    if(success):
        logger.debug("find_roi_bam complete successfully for "+chr + event) 
    return           
    
def mutate_reads(bamsortfn,chr, event):
    
    fn,sortbyname,sortbyCoord, bedfn = init_file_names(chr, event, tmpbams_path, haplotype_path)
    cmd=" ".join(["sort -u", bedfn, "-o", bedfn]); runCommand(cmd)
    outbamfn = sub('.sorted.bam$',".mutated_het.bam", bamsortfn)
    outbamsortfn = sub('.sorted.bam$',".mutated_het.sorted", bamsortfn)
    
    mergedsortfn = sub('.sorted.bam$',".mutated_merged.sorted.bam", bamsortfn)
    try:
        if not terminating.is_set():
            
            if(os.path.isfile(bamsortfn) and os.path.isfile(bedfn) ):
                samfile = pysam.Samfile(bamsortfn, "rb" )
                alignmentfile = pysam.AlignmentFile(bamsortfn, "rb" )
                outbam = pysam.Samfile(outbamfn, 'wb', template=samfile) 
                bedfile = open(bedfn, 'r')
                covpath = "/".join([haplotype_path, "written_coverage_het.txt"])
                covfile = open(covpath, 'w')
                snpratiopath = "/".join([haplotype_path, "het_snp_ratio.txt"])
                snpaltratiofile = open(snpratiopath,'w')
                writtenreads = []
                
                for bedline in bedfile:
                    c = bedline.strip().split()
                    if (len(c) == 6 ):
                        chr2 = c[0]; chr = c[0].strip("chr"); start = int(c[1]);end = int(c[2])
                        refbase = str(c[3]); altbase = str(c[4]); haplotype = str(c[5])
                    else:
                        continue
                    
                    readmappings = alignmentfile.fetch(chr2, start, end)          
                    num_reads_written = 0
                    for shortread in readmappings:
                        try:
                            mate = alignmentfile.mate(shortread)
                        except:
                            continue
                       
                        if(shortread.is_paired and shortread.is_proper_pair and not shortread.is_duplicate  
                           and not shortread.is_secondary and not shortread.qname in writtenreads and shortread.mapping_quality >= 30
                           and mate.mapping_quality >= 30 and not mate.is_duplicate and mate.is_proper_pair and not mate.is_secondary):

                            try:
                                index = shortread.get_reference_positions().index(start)
                                tmpread = shortread.query_sequence
                                mutated_hap1 = tmpread[:index] +  altbase + tmpread[index + 1:]
                                mutated_hap2 = tmpread[:index] +  refbase + tmpread[index + 1:]
                                if(haplotype == "hap1"):
                                    shortread.query_sequence = mutated_hap1 
                                elif(haplotype == "hap2"):
                                    shortread.query_sequence = mutated_hap2
                            except:
                                continue
                             
                            try:
                                index_mate = mate.get_reference_positions().index(start)
                                nuecleotide_mate = mate.seq[index_mate]
                                tmpread_mate= mate.query_sequence
                                mutated_mate_hap1 = tmpread_mate[:index_mate] +  altbase + tmpread_mate[index_mate + 1:]
                                mutated_mate_hap2 = tmpread_mate[:index_mate] +  refbase + tmpread_mate[index_mate + 1:]
                                if(haplotype == "hap1"):
                                     mate.query_sequence = mutated_mate_hap1
                                elif(haplotype == "hap2"):
                                     mate.query_sequence = mutated_mate_hap2
                            except (KeyError,ValueError) as e :
                                pass

                            outbam.write(shortread)
                            outbam.write(mate)
                            writtenreads.append(shortread.qname)
                            num_reads_written += 2
                            continue
                outbam.close()
                covfile.close()
                snpaltratiofile.close()         
                sortBam(outbamfn,outbamsortfn+'.bam')
                bamDiff(bamsortfn,outbamsortfn+'.bam', tmpbams_path )
                
                merge_bams("/".join([tmpbams_path, 'diff_only1_' +  os.path.basename(bamsortfn)]), outbamsortfn+'.bam', mergedsortfn)
                os.remove("/".join([tmpbams_path,  'diff_only1_' +  os.path.basename(bamsortfn)]))
                os.remove(outbamfn)
                os.remove(outbamsortfn+'.bam')
    
    except (KeyboardInterrupt):
        logger.error('Exception Crtl+C pressed in the child process  in mutaute_reads')
        terminating.set()
        return
    except Exception as e:   
        logger.exception("Exception in mutate_reads %s" ,e )
        terminating.set()
        return
    return        

#cn change is 1 for CN=1,3 and 2 for CN=0,4
def calculate_sample_rate(inbam, outbam, cnchange, purity):
    logger.debug("___ adjusting sample rate ___")

def implement_cnv(chromosome_event):
    chr,event = chromosome_event .split("_")
    
    logger.debug("___ Bamgineer main engine started ___")
    success = True
    try:
        if not terminating.is_set():
            bamfn,sortbyname,sortbyCoord, bedfn = init_file_names(chr, event, tmpbams_path, haplotype_path)
            bamsortfn = sub('.bam$', '.sorted.bam', bamfn)
          
            if(event== 'gain'):
                    bamrepairedsortfn = sub('.sorted.bam$', ".re_paired.sorted.bam", bamsortfn)
                    mergedsortfn = sub('.sorted.bam$',".mutated_merged.sorted.bam", bamrepairedsortfn)
                    mergedrenamedfn = sub('.sorted.bam$',".mutated_merged_renamed.sorted.bam", bamrepairedsortfn)
            
                    GAIN_FINAL = "/".join([finalbams_path,  str(chr).upper() +'_GAIN.bam'])
                    if(os.path.isfile(bamsortfn)):
                        re_pair_reads(bamsortfn)
                        mutate_reads(bamrepairedsortfn, chr, 'gain')
                        renamereads(mergedsortfn, mergedrenamedfn)
                        
                        ratio_kept1 = float(countReads(bamsortfn))/float(countReads(bamfn))
                        ratio_kept2 = float(countReads(bamrepairedsortfn))/float(countReads(bamsortfn))
                        
                        samplerate= round(0.5/(ratio_kept1*ratio_kept2*0.98),2)
                        #logger.debug("ratios kept for:"+ ntpath.basename(bamsortfn)+ ": "+ str(ratio_kept1) + "  "+ str(ratio_kept2))
                        os.remove(bamfn)
                        if(samplerate < 1.0):
                            subsample(mergedrenamedfn, GAIN_FINAL,str(samplerate)) #calculate it later
                            logger.debug("___ sampling rate for " + ntpath.basename(bamsortfn)  +" : "+ str(samplerate))
                        elif(samplerate > 1.0 and samplerate< 1.1):
                            os.rename(mergedrenamedfn, GAIN_FINAL)
                        else:
                            logger.error('not enough reads for '+ntpath.basename(bamsortfn)+ 'rate: '+str(samplerate) )
                            success = False
                            return
            elif(event== 'loss'):
               
                inbam_deletion = "/".join([finalbams_path , str(chr).upper() + '_LOSS.bam'])
                if(os.path.isfile(bamsortfn)):
                    
                    mutate_reads(bamsortfn, chr, 'loss')
                    mergedsortfn = sub('.sorted.bam$',".mutated_merged.sorted.bam", bamsortfn)
                    mergedsortsampledfn = sub('.sorted.bam$',".mutated_merged.sampled.sorted.bam", bamsortfn)
                    
                    ratio_kept = float(countReads(bamsortfn))/float(countReads(bamfn))
                    samplerate= round(0.5/(ratio_kept*0.98),2)
                    LOSS_FINAL = "/".join([finalbams_path,  str(chr).upper() +'_LOSS.bam'])
                    logger.debug("ratios kept for:"+ ntpath.basename(bamsortfn)+ ": "+ str(ratio_kept))
                    subsample(mergedsortfn, mergedsortsampledfn,str(samplerate)) 
                    bamDiff(sortbyCoord, mergedsortsampledfn, tmpbams_path)
                    os.rename("/".join([tmpbams_path,  'diff_only1_' + chr + '.bam']), LOSS_FINAL)

                elif(not os.path.isfile(inbam_deletion) and os.path.isfile(sortbyCoord) ):# if it exists from previous runs 
                    
                    os.symlink(sortbyCoord, inbam_deletion)
                    
    except (KeyboardInterrupt):
        logger.error('Exception Crtl+C pressed in the child process  in find_roi_bam for chr ' +chr + event)
        terminating.set()
        success=False
        return
    except Exception as e:   
        logger.exception("Exception in find_roi_bam %s" ,e )
        terminating.set()
        success=False
        return
    if(success):
        logger.debug("implement_cnv complete successfully for "+chr + event) 
    return           

  
def re_pair_reads(bamsortfn):    
    try:
        if not terminating.is_set():
            logger.debug(" calling  re-pair-reads version" )
            bamrepairedfn = sub('.sorted.bam$',  ".re_paired.bam", bamsortfn)
            bamrepairedsortfn = sub('.sorted.bam$', ".re_paired.sorted.bam", bamsortfn)
            
            if(os.path.isfile(bamsortfn)):
                  
                inbam = pysam.Samfile(bamsortfn, 'rb')
                outbam = pysam.Samfile(bamrepairedfn, 'wb', template=inbam)  
    
                writtencount = 0
             
                #positive & negative strands
                strands=['pos','neg']
                
                for strand in strands :
                    read1fn= sub('.bam$', '.read1_'+strand+'.bam', bamsortfn)
                    read2fn= sub('.bam$', '.read2_'+strand+'.bam', bamsortfn)
                    
                    if(not os.path.isfile(read1fn) or not os.path.isfile(read2fn)):
                        splitPairAndStrands(bamsortfn)
                    
                    splt1 = pysam.Samfile(read1fn , 'rb')
                    splt2 = pysam.Samfile(read2fn , 'rb')
                    itr1 =   splt1.fetch(until_eof=True)
                    itr2 =   splt2.fetch(until_eof=True)
                    start = True
                    
                    
                    for read1, read2 in  izip(itr1, itr2):                   
                       
                        try:
                            if(read2.qname != read1.qname and start):
                                read2 = itr2.next()
                                start = False
                                continue
                            
                            read1next=itr1.next()
                            read2next=itr2.next()
                
                            if(strand == 'pos'):
                                
                                tlenabs1 = read2next.pos - read1.pos + abs(read2next.qlen)
                                tlenabs2 =  read2.pos - read1next.pos  + abs(read2.qlen)  
                                tlenmean = (abs(read1.tlen) + abs(read1next.tlen))/2
                                
                                if(tlenabs1 > 0.2*tlenmean and tlenabs1 < 5*tlenmean and read2next.qname != read1.qname and tlenabs1 > 0 and
                                   not read1.is_duplicate and not read1.is_secondary and not read2next.is_duplicate and not read2next.is_secondary):
                                    
                                    read1.tlen = tlenabs1
                                    read2next.tlen = -tlenabs1
                                    read1.pnext = read2next.pos
                                    read2next.pnext = read1.pos
                                    read2next.qname = read1.qname 
                                    outbam.write(read1)
                                    outbam.write(read2next)
                                    writtencount = writtencount + 1
                                    
                                if(tlenabs2 > 0.2*tlenmean and tlenabs2 < 5*tlenmean and read1next.qname != read2.qname and tlenabs2 > 0 and
                                   not read2.is_duplicate and not read2.is_secondary and not read1next.is_duplicate and not read1next.is_secondary ):
                              
                                    read1next.tlen = tlenabs2
                                    read2.tlen = -tlenabs2 
                                    read2.pnext = read1next.pos
                                    read1next.pnext = read2.pos
                                    read2.qname = read1next.qname
                                    outbam.write(read1next)
                                    outbam.write(read2)
                                    writtencount = writtencount + 1  
                            
                            elif(strand== 'neg'):
                                
                                tlenabs1 = read1.pos - read2next.pos + abs(read1.qlen)
                                tlenabs2 = read1next.pos -read2.pos + abs(read1next.qlen)
                                tlenmean = (abs(read1.tlen) + abs(read1next.tlen))/2
                                
                                if(tlenabs1 > 0.2*tlenmean and tlenabs1 < 5*tlenmean and read2next.qname != read1.qname and tlenabs1 > 0 and
                                   not read1.is_duplicate and not read1.is_secondary and not read2next.is_duplicate and not read2next.is_secondary):
                                    
                                    read1.tlen = -tlenabs1
                                    read2next.tlen = tlenabs1
                                    read1.pnext = read2next.pos
                                    read2next.pnext = read1.pos
                                    read2next.qname = read1.qname
                                    outbam.write(read1)
                                    outbam.write(read2next)
                                    writtencount = writtencount + 1
                                if(tlenabs2 > 0.2*tlenmean and tlenabs2 < 5*tlenmean and read1next.qname != read2.qname and tlenabs2 > 0 and
                                    not read2.is_duplicate and not read2.is_secondary and not read1next.is_duplicate and not read1next.is_secondary):
                            
                                    read1next.tlen = -tlenabs2
                                    read2.tlen = tlenabs2
                                    read2.pnext = read1next.pos
                                    read1next.pnext = read2.pos
                                    read2.qname = read1next.qname
                                    outbam.write(read1next)
                                    outbam.write(read2)
                                    writtencount = writtencount + 1

                        except StopIteration:
                            break        
               
                    splt1.close();splt2.close()
                    os.remove(read1fn)
                    os.remove(read2fn)
                    
                inbam.close()
                outbam.close() 
                
                sortBam(bamrepairedfn, bamrepairedsortfn)
                os.remove(bamrepairedfn)   
                
            
    except (KeyboardInterrupt):
        logger.error('Exception Crtl+C pressed in the child process  in re_pair_reads')
        terminating.set()
        return False
    except Exception as e:   
        logger.exception("Exception in re_pair_reads %s" ,e )
        terminating.set()
        return False
    return             
  
      
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
    bamDiff(inbamsorted+'.bam', outbamsorted +'.bam', path )
    outbam.close()           

def removeIfEmpty(bamdir,file):
    try:
        if not terminating.is_set():   
            if file.endswith(".bam"):
               command = " ".join(["samtools view", "/".join([bamdir, file]), "| less | head -1 | wc -l" ])
               nline= subprocess.check_output(command, shell = True)  
               if (os.path.isfile( "/".join([bamdir, file])) and (int(nline) == 0)):
                       os.remove("/".join([bamdir, file]))
                       logger.debug(' removing ' + "/".join([bamdir, file]))
    except (KeyboardInterrupt):
        logger.error('Exception Crtl+C pressed in the child process  in removeIfEmpty ')
        terminating.set()
        return
    except Exception as e:   
        logger.exception("Exception in removeIfEmpty %s" ,e )
        terminating.set()
        return
    return                       

def run_pipeline(results_path):
    global haplotype_path,cancer_dir_path,tmpbams_path, finalbams_path,log_path, logfile ,terminating,logger,logQueue
    haplotype_path,cancer_dir_path,tmpbams_path, finalbams_path,log_path, logfile = handle.GetProjectPaths(results_path)
    terminating,logger,logQueue = handle.GetLoggings(logfile)
    
    t0 = time.time()
    outbamfn=params.GetOutputFileName() 
    chromosome_event=create_chr_event_list()
    chromosomes_bamfiles = create_chr_bam_list()
    logger.debug('pipeline started!')
    
    initialize(results_path,haplotype_path,cancer_dir_path)
    pool1 = multiprocessing.Pool(processes=4, initializer=initPool, initargs=[logQueue, logger.getEffectiveLevel(), terminating] ) 
    try:
        result1 = pool1.map_async(find_roi_bam, chromosome_event ).get(9999999)
        result2 = pool1.map_async(implement_cnv, chromosome_event ).get(9999999)
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
    mergeSortBamFiles(outbamfn, finalbams_path )
    t1 = time.time()
    shutil.rmtree(tmpbams_path)
    logger.debug(' ***** pipeline finished in ' + str(round((t1 - t0)/60.0, 1)) +' minutes ***** ')
    logging.shutdown()
