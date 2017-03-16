import os
import sys
import pysam
import pybedtools
from re import sub
import shutil
import subprocess
from itertools import izip
from uuid import uuid4

def splitAndRePairHET(inbamfn, outbamfn, chr):
    print(" calling splitAndRePair het version" )
    splitfn1 = '/'.join([splittmpbams_hap,chr+'_sp1_h.bam'])
    splitfn2 = '/'.join([splittmpbams_hap,chr+'_sp2_h.bam'])
    
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

    num_reads_to_write = 0
    for row in spltcount:
        num_reads_to_write+=1
    
    itr1 = splt1.fetch(until_eof=True)
    itr2 = splt2.fetch(until_eof=True)
    
    num_reads_to_write = 0
    for row in spltcount:
        num_reads_to_write+=1
    
    writtencount = 0
    samplerate = 0
    start = True
    for read1, read2 in izip(itr1, itr2):
        
        try:
           
            if(read2.qname != read1.qname and start and read2.pos < read1.pos ):
                read2 = itr2.next()
                start = False
                continue
            read1next = itr1.next()
            read2next = itr2.next()     
        
        except StopIteration:
            pass   
            
        if(read2.rnext == read2.tid and read1.rnext == read1.tid and read1.qname != read2.qname and read2.tid == read1.tid and
           read2next.rnext == read2next.tid and read1next.rnext == read1next.tid and read1next.qname != read2next.qname and read2next.tid == read1next.tid ): 
            
            if ( abs(read1.pos - read2next.pos) < 200 and abs(read1.pos - read2next.pos) > 1  and
                 read1.mapping_quality >= 30 and read2next.mapping_quality >= 30 ):
                
                read2next.qname = read1.qname
                outbam.write(read1)
                outbam.write(read2next)
                writtencount = writtencount + 1
                
            if(abs(read2.pos - read1next.pos)  < 200 and  abs(read2.pos - read1next.pos)  > 1 and
                 read2.mapping_quality >= 30 and read1next.mapping_quality >= 30):
                
                read2.qname = read1next.qname
                outbam.write(read1next)
                outbam.write(read2)
                writtencount = writtencount + 1
    
    inbam.close()         
    inbam.close()         
    splt1.close()
    splt2.close()
    outbam.close()
    if(num_reads_to_write > 0):
        percentkept = float(writtencount)/float(num_reads_to_write)
        print("adjusted sampling rate for Non_Het: " + str(0.5/percentkept))
        return(0.5/percentkept)  
   
    
    
#################
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


def subsample(bamfn1, bamfn2, samplingrate = 0.5):
    command = " ".join(["samtools view -s", samplingrate ,"-b", bamfn1, ">", bamfn2])
    subprocess.check_output(command, shell = True)

def splitPairs(inbamfn,pair1fn, pair2fn):
    command1 = " ".join(["samtools view -u -h -f 0x0040", inbamfn, ">", pair1fn])
    command2 = " ".join(["samtools view -u -h -f 0x0080", inbamfn, ">", pair2fn])
    subprocess.check_output(command1, shell = True)
    subprocess.check_output(command2, shell = True)
    
    
def splitAndRePairNONHET(inbamfn, outbamfn, chr):
    print(" calling splitAndRePair non_het version" )
    splitfn1 = '/'.join([splittmpbams_hap,chr+'_sp1_nh.bam'])
    splitfn2 = '/'.join([splittmpbams_hap,chr+'_sp2_nh.bam'])
    
    splitPairs(inbamfn, splitfn1, splitfn2 )
    inbam = pysam.Samfile(inbamfn, 'rb')
    splt1sortedfn = sub('.bam$', '.sorted', splitfn1)
    splt2sortedfn = sub('.bam$', '.sorted', splitfn2)
    
    pysam.sort(splitfn1 , splt1sortedfn )
    pysam.sort(splitfn2 , splt2sortedfn)
    #os.remove(splitfn1)
    #os.remove(splitfn2)
    
    
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
               
            if(abs(read1.pos - read2next.pos) < 500 and abs(read1.pos - read2next.pos) > 10 and read1.mapping_quality >= 30 and read2next.mapping_quality >= 30):
               
                read2.qname = read1next.qname
                outbam.write(read1)
                outbam.write(read2next)
                writtencount = writtencount + 1
        
            if(abs(read1next.pos - read2.pos) < 500 and abs(read1next.pos - read2.pos) > 10 and read2.mapping_quality >= 30 and read1next.mapping_quality >= 30):
                read2next.qname = read1.qname 
                outbam.write(read1next)
                outbam.write(read2)
                writtencount = writtencount + 1
   
    inbam.close()         
    splt1.close()
    splt2.close()
    outbam.close()
    
    if(num_reads_to_write > 0):
        percentkept = float(writtencount)/float(num_reads_to_write)
        print("adjusted sampling rate for HET  " + str(0.5/percentkept))
        return(0.5/percentkept)
    

chr = sys.argv[1]
event = sys.argv[2]
splittmpbams_hap = sys.argv[3]
finalbams_hap = sys.argv[4]
NONHET = "/".join([splittmpbams_hap,  (chr + '_'+event).upper()  +'_NONHET.bam'])
HET_ALT =  "/".join([splittmpbams_hap, (chr +'_'+ event).upper() +'_HET_ALT.bam'])
            
gain_HET_ALT_RE_PAIR =  "/".join([splittmpbams_hap,  chr + '_GAIN_HET_ALT_RE_PAIR.bam'])
gain_HET_ALT_RE_PAIR_SAMPLED =  "/".join([splittmpbams_hap, chr + event +'_HET_ALT_RE_PAIR_SAMPLED.bam'])
gain_HET_ALT_RE_PAIR_SAMPLED_RENAMED =  "/".join([splittmpbams_hap ,  chr +'_GAIN_HET_ALT_RE_PAIR_SAMPLED_RENAMED.bam'])

gain_NONHET_RE_PAIR = "/".join([splittmpbams_hap, chr + event +'_NONHET_REPAIR.bam'])
gain_NONHET_RE_PAIR_SAMPLED = "/".join([splittmpbams_hap, chr + event +'_NONHET_SAMPLED.bam'])
gain_NONHET_RE_PAIR_SAMPLED_RENAMED = "/".join([splittmpbams_hap, chr + '_GAIN_NONHET_RE_PAIR_SAMPLED_RENAMED.bam'])

gain_NONHET_FINAL = "/".join([finalbams_hap,  chr.upper() +'_GAIN_NH'])
gain_HET_FINAL = "/".join([finalbams_hap,  chr.upper() +'_GAIN_H'])

if(os.path.isfile(HET_ALT)):
   
    samapleratehet = splitAndRePairHET(HET_ALT,  gain_HET_ALT_RE_PAIR, chr )
    subsample(gain_HET_ALT_RE_PAIR,gain_HET_ALT_RE_PAIR_SAMPLED, str(samapleratehet)) # we need to keep a bit more (by 15-20%)
    renamereads(gain_HET_ALT_RE_PAIR_SAMPLED, gain_HET_ALT_RE_PAIR_SAMPLED_RENAMED )
    pysam.sort(gain_HET_ALT_RE_PAIR_SAMPLED_RENAMED,  gain_HET_FINAL)
    #os.remove(gain_HET_ALT_RE_PAIR_SAMPLED_RENAMED)
    #logger.debug('sampling rate for HET'+chr +' = '+ str(samapleratehet))
else:
    print('HET_ALT FILE NOT EXISTING')
       
if(os.path.isfile(NONHET)):
    samapleratenonhet = splitAndRePairNONHET(NONHET, gain_NONHET_RE_PAIR, chr)
    subsample(gain_NONHET_RE_PAIR, gain_NONHET_RE_PAIR_SAMPLED,  str(samapleratenonhet))
    renamereads(gain_NONHET_RE_PAIR_SAMPLED, gain_NONHET_RE_PAIR_SAMPLED_RENAMED)
    pysam.sort(gain_NONHET_RE_PAIR_SAMPLED_RENAMED, gain_NONHET_FINAL)
    #os.remove(gain_NONHET_RE_PAIR_SAMPLED_RENAMED)
else:
    print('NON_HET FILE NOT EXISTING')
    
    