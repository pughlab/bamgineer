import pysam
import ntpath
from helpers import parameters as params
from helpers import handlers as handle
from helpers import bamgineerHelpers as bamhelp
from utils import *

import multiprocessing, threading, logging, sys, traceback,  StringIO, Queue

global bases
bases = ('A','T','C','G')

##temporary
purity = [ '0.2', '0.4', '0.6','0.8','1.0']
cancer_list = ['brca','crc','gbm']

def initPool(queue, level, terminating_):
    """
    This causes the logging module to be initialized with the necessary info
    in pool threads to work correctly.
    """
    logging.getLogger('').setLevel(level)
    terminating = terminating_

#chr 21 and 22 for test, change it to 1
def create_chr_event_list():
    chrom_event= []
    for c in range(21,23):
        for e in ['gain','loss']:
            chev = "_".join(['chr'+str(c), e])
            chrom_event.append(chev)
    return chrom_event        

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
            exonsinhetbed = "/".join([haplotype_path,  event + "_exons_in_het.bed"])
            nonhetbed = "/".join([haplotype_path, event + "_non_het.bed"])
            hetbed = "/".join([haplotype_path, event + "_het.bed"])
            hetsnpbed = "/".join([haplotype_path,  event + "_het_snp.bed"])
            intersectBed( exons_path, locals()[event + 'cnv'], exonsinroibed, wa=True)
            intersectBed( exonsinroibed, phased_bed, exonsinhetbed, wa=True)  
            intersectBed(phased_bed, exonsinroibed, hetsnpbed, wa=True)
             
            subtractBeds(exonsinroibed, phased_bed , nonhetbed) # so that 
            intersectBed( phased_bed, exonsinroibed, hetbed, wa=True)
            
            splitBed(hetsnpbed, event+'_het_snp_')
            splitBed(hetbed, event+'_het_')
            splitBed(nonhetbed, event+'_non_het_')             
        
    except:
        logger.error('Initialization error !')
        raise
    
    logger.debug("--- initialization complete ---")    
    return 

def mutate_reads(bedfn, bamfn, outfn, haplotypedir):

    logger.debug("___ mutating reads and finding reads not matching hg19 at germline SNP locations ___")
    
    try:
        if not terminating.is_set():
            samfile = pysam.Samfile(bamfn, "rb" )
            
            path, filename = os.path.split(bedfn)
            outbam = pysam.Samfile(outfn, 'wb', template=samfile) 
            
            refbamfn = sub('_het_alt_roi.bam$',"_ref_reads.bam", outfn)
            altbamfn = sub('_het_alt_roi.bam$',"_alt_reads.bam", outfn)
            hapbamfn = sub('_het_alt_roi.bam$',"_hap.bam", outfn)
               
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
                                        
                                        refbam.write(shortread)
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
                                        
                                        elif(shortread.is_read1 and shortread.pos + index < mate.pos):
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
                                        
                                        elif(shortread.is_read2 and shortread.pos + index > mate.pos+ 101):
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
                                        elif(shortread.is_read1 and shortread.pos + index < mate.pos):
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
                                        elif(shortread.is_read2 and shortread.pos + index > mate.pos+ 101):    
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
                        pass  
                                       
            outbam.close()
            refbam.close()
            altbam.close()
            hapbam.close()
            covfile.close()
            snpaltratiofile.close()
    except (KeyboardInterrupt):
        logger.error('Exception Crtl+C pressed in the child process  in mutaute_reads')
        terminating.set()
        return
    except Exception as e:   
        logger.exception("Exception in mutate_reads %s" ,e )
        terminating.set()
        return
    return        

    
def init_file_names(chr, event,tmpbams_path, haplotypedir):
    
    flist=[]
    splitbams = params.GetSplitBamsPath()
    hetroibam = "/".join([tmpbams_path ,chr + event +"_het_roi.bam"])
    hetaltroibam = "/".join([tmpbams_path, chr + event + "_het_alt_roi.bam"])
    nonhetroibam = "/".join([tmpbams_path, chr + event +  "_non_het_roi.bam"])
    overlapsbam  = "/".join([tmpbams_path, chr + event + "_overlaps.bam"])
    
    sortbyname =  "/".join([splitbams,  chr + '.byname.bam'])
    sortbyCoord = "/".join([splitbams,  chr + '.bam'])
    hetaltroibamsorted = sub('.bam$','.sorted',hetaltroibam)
    refbam = "/".join([tmpbams_path, chr + event +"_ref_reads.bam"])
    altbam = "/".join([tmpbams_path, chr + event +"_alt_reads.bam"])
    refbamsorted = sub('.bam$','.sorted', refbam)
    altbamsorted = sub('.bam$','.sorted', altbam)
    
    hapbam = "/".join([tmpbams_path, chr + event +"_hap.bam"])
    hapbamsorted = sub('.bam$','.sorted', hapbam)
    hetroisortbam=sub('.bam$','.sorted',hetroibam)
    nonhetroibamsorted = sub('.bam$','.sorted',nonhetroibam)
    hetbed   = "/".join([haplotypedir, event+'_het_' + chr + '.bed'])
    nonhetbed   = "/".join([haplotypedir, event+'_non_het_' + chr + '.bed'])
    hetsnp   = "/".join([haplotypedir, event+'_het_snp_' + chr + '.bed'])
    
    flist.extend([hetroibam,hetroisortbam, hetaltroibam,nonhetroibam, overlapsbam,sortbyname,sortbyCoord,hetaltroibamsorted,refbam,altbam,refbamsorted, altbamsorted,hapbam,hapbamsorted,nonhetroibamsorted,hetbed,nonhetbed,hetsnp])
    return flist

def find_roi_bam(chromosome_event):
    chr,event = chromosome_event .split("_")
    
    hetroi,hetroisort,hetaltroi,nonhetroi, overlaps,sortbyname,sortbyCoord,hetaltsort,ref,alt,refsort, altsort,hap,hapsort,nonhetsort,hetbed,nonhetbed,hetsnp = init_file_names(chr, event, tmpbams_path, haplotype_path)
    HET,NHET = handle.GetHetBamPaths(tmpbams_path,chr, event)
    success = True
    try:
        if not terminating.is_set():
            
            extractPairedReadfromROI(sortbyname,hetbed, hetroi)
            extractPairedReadfromROI(sortbyname, nonhetbed, nonhetroi)
            removeIfEmpty(tmpbams_path,ntpath.basename(hetroi))
            if (os.path.isfile(hetroi)):
                mutate_reads(hetsnp, hetroi, hetaltroi,haplotype_path)
                pysam.sort(hetaltroi,hetaltsort)
                os.remove(hetaltroi)
            
            removeIfEmpty(tmpbams_path,ntpath.basename(nonhetroi))
            if (os.path.isfile(nonhetroi) and os.path.isfile(hetbed)):    
                pysam.sort(nonhetroi,nonhetsort)
                pysam.index(nonhetsort+'.bam')
                os.remove(nonhetroi)
            
                bamDiff(nonhetsort+'.bam' , hetaltsort+'.bam' , tmpbams_path )
                os.rename(hetaltsort+'.bam', HET)
                
                if (os.path.isfile("/".join([tmpbams_path,  'diff_only1_'+ chr+ event + '_non_het_roi.sorted.bam']))):
                    os.rename("/".join([tmpbams_path, 'diff_only1_'+ chr+ event + '_non_het_roi.sorted.bam']), NHET )
                else:
                    os.rename(nonhetsort+'.bam',NHET)
        
                   
            
    except (KeyboardInterrupt):
        logger.error('Exception Crtl+C pressed in the child process  in find_roi_bam for chr ' + chr)
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


def implement_gain_loss( chromosome_event ):
    chr,event = chromosome_event .split("_")

    hetroi,hetroisort,hetaltroi,nonhetroi, overlaps,sortbyname,sortbyCoord,hetaltsort,ref,alt,refsort, altsort,hap,hapsort,nonhetsort,hetbed,nonhetbed,hetsnp = init_file_names(chr, event, tmpbams_path, haplotype_path)
    HET,NHET = handle.GetHetBamPaths(tmpbams_path,chr, event)
    logger.debug(str(HET))
   
    try:
        if not terminating.is_set():
            if (event == 'gain'):
                logger.debug('in gain event ')
                gain_HET_ALT_RE_PAIR =  "/".join([tmpbams_path,  str(chr) + '_G_H_R.bam'])
                gain_HET_ALT_RE_PAIR_SAMPLED =  "/".join([tmpbams_path, str(chr) + '_G_H_R_S.bam'])
                gain_HET_ALT_RE_PAIR_SAMPLED_RENAMED =  "/".join([tmpbams_path ,  str(chr) +'_G_H_R_S_Re.bam'])
                
                gain_NONHET_RE_PAIR = "/".join([tmpbams_path, str(chr) +'_G_NH_R.bam'])
                gain_NONHET_RE_PAIR_SAMPLED = "/".join([tmpbams_path, str(chr) + '_G_NH_R_S.bam'])
                gain_NONHET_RE_PAIR_SAMPLED_RENAMED = "/".join([tmpbams_path, str(chr) + '_G_NH_R_S_Re.bam'])
                
                gain_NONHET_FINAL = "/".join([finalbams_path,  str(chr).upper() +'_GAIN_NH'])
                gain_HET_FINAL = "/".join([finalbams_path,  str(chr).upper() +'_GAIN_H'])
                
                if(os.path.isfile(HET)):
                    
                    samapleratehet = splitAndRePair(HET,  gain_HET_ALT_RE_PAIR, chr , "het")
                    if(samapleratehet < 0.9):
                        subsample(gain_HET_ALT_RE_PAIR,gain_HET_ALT_RE_PAIR_SAMPLED, str(samapleratehet*1.1))  
                    else:
                        gain_HET_ALT_RE_PAIR_SAMPLED = gain_HET_ALT_RE_PAIR
                        
                    renamereads(gain_HET_ALT_RE_PAIR_SAMPLED, gain_HET_ALT_RE_PAIR_SAMPLED_RENAMED )
                    pysam.sort(gain_HET_ALT_RE_PAIR_SAMPLED_RENAMED,  gain_HET_FINAL)
                    logger.debug('adjusted sampling rate for HET'+chr +' = '+ str(samapleratehet*1.1))
                       
                if(os.path.isfile(NHET)):
                    samapleratenonhet = splitAndRePair(NHET, gain_NONHET_RE_PAIR, chr, "nonhet")
                    if(samapleratenonhet < 1.0):
                        subsample(gain_NONHET_RE_PAIR, gain_NONHET_RE_PAIR_SAMPLED,  str(samapleratenonhet))
                    else:
                        gain_NONHET_RE_PAIR_SAMPLED = gain_NONHET_RE_PAIR
                    renamereads(gain_NONHET_RE_PAIR_SAMPLED, gain_NONHET_RE_PAIR_SAMPLED_RENAMED)
                    pysam.sort(gain_NONHET_RE_PAIR_SAMPLED_RENAMED, gain_NONHET_FINAL)
                    
            elif (event == 'loss'):
                inbam_deletion = "/".join([finalbams_path , str(chr).upper() + 'LOSS_FINAL.bam'])
                if(not os.path.isfile(NHET) and not os.path.isfile(HET)): 
                    os.symlink(sortbyCoord, inbam_deletion)
                    #logger.debug('No loss event defined for chromosome: ' + chr +' creating symlink for ')
                
                if(os.path.isfile(NHET) and os.path.isfile(HET)):
                    loss_NONHET_SAMPLED =   "/".join([tmpbams_path,  chr + event +'_NONHET_SAMPLED.bam'])
                    loss_HET_SAMPLED =   "/".join([tmpbams_path,  chr + event +'_HET_SAMPLED.bam'])
                    
                    #logger.debug("In LOSS module NONHET is " + NHET )
                    subsample(NHET ,loss_NONHET_SAMPLED , str(0.5))
                    inbam_minus_loss_NONHET_ALT = "/".join([tmpbams_path, chr + event+ '_MINUS_NONHET_ALT.bam'])
                    
                    bamDiff(sortbyCoord , loss_NONHET_SAMPLED, tmpbams_path)
                    os.rename("/".join([tmpbams_path,  'diff_only1_' + chr + '.bam']), inbam_minus_loss_NONHET_ALT)
                    #logger.debug(" loss module renaming:   "+ "/".join([tmpbams_path, 'diff_only1_' + chr + '.bam'])+
                    #             '   '+inbam_minus_loss_NONHET_ALT)
                    pysam.sort(hap, hapsort)
                    bamDiff(inbam_minus_loss_NONHET_ALT, hapsort+'.bam', tmpbams_path)
                    os.rename("/".join([tmpbams_path,'diff_only1_'+ chr + event+ '_MINUS_NONHET_ALT.bam']),inbam_deletion)
                       
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
    splitfn1 = '/'.join([tmpbams_path,chr+'sp1_nh.bam'])
    splitfn2 = '/'.join([tmpbams_path,chr+'sp2_nh.bam'])
    
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
    print('running pipeline')
    global haplotype_path,cancer_dir_path,tmpbams_path, finalbams_path,log_path, logfile
    global terminating,logger,logQueue
    
    haplotype_path,cancer_dir_path,tmpbams_path, finalbams_path,log_path, logfile = handle.GetProjectPaths(results_path)
    terminating,logger,logQueue = handle.GetLoggings(logfile)
    t0 = time.time()
    outbamfn=params.GetOutputFileName() 
    chromosome_event=create_chr_event_list()
    initialize(results_path,haplotype_path,cancer_dir_path)
    
    pool1 = multiprocessing.Pool(processes=16, initializer=initPool, initargs=[logQueue, logger.getEffectiveLevel(), terminating] ) 
    try:
        result1 = pool1.map_async(find_roi_bam, chromosome_event ).get(9999999)
        result2 = pool1.map_async(implement_gain_loss, chromosome_event ).get(9999999)
      
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
    logging.shutdown()
    t1 = time.time()
    #shutil.rmtree(.tmpbams)
    logger.debug(' ***** Multi-processing phase took %f ' + str((t1 - t0)/60.0) +' minutes to finish ***** ')

    
    




