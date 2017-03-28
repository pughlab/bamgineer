import sys
import pysam
import subprocess
from re import sub
import os
from itertools import izip

bamsortfn = sys.argv[1]

bamrepairedfn = sub('.sorted.bam$',  ".repaired.bam", bamsortfn)
bamrepairedsortfn = sub('.sorted.bam$', ".repaired.sorted.bam", bamsortfn)

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
        
        read1_strand1sortfn =  sub('.bam$', '.read1_pos.bam', bamsortfn)
        read1_strand2sortfn =  sub('.bam$', '.read1_neg.bam', bamsortfn)
        read2_strand1sortfn =  sub('.bam$', '.read2_pos.bam', bamsortfn)
        read2_strand2sortfn =  sub('.bam$', '.read2_neg.bam', bamsortfn)
        command1 = " ".join(["samtools view -u -h -f 0x0063", bamsortfn, ">", read1_strand1sortfn])
        command2 = " ".join(["samtools view -u -h -f 0x0053", bamsortfn, ">", read1_strand2sortfn])
        command3 = " ".join(["samtools view -u -h -f 0x0093", bamsortfn, ">", read2_strand1sortfn])
        command4 = " ".join(["samtools view -u -h -f 0x00A3", bamsortfn, ">", read2_strand2sortfn])
        subprocess.check_output(command1, shell = True)
        subprocess.check_output(command2, shell = True)
        subprocess.check_output(command3, shell = True)
        subprocess.check_output(command4, shell = True)
  
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
      
      
    command = " ".join(["sambamba sort", bamrepairedfn, "-o", bamrepairedsortfn])  
    subprocess.check_output(command, shell = True)
    os.remove(bamrepairedfn) 