import sys
import pysam
import subprocess
from re import sub
import os
import numpy as np

bamrepairedsortfn = sys.argv[1]
bedfn = sys.argv[2]
haplotype_path= sys.argv[3]

outhetfn = sub('.sorted.bam$',".mutated.het.bam", bamrepairedsortfn)

if(os.path.isfile(bamrepairedsortfn) and os.path.isfile(bedfn) ):
    samfile = pysam.Samfile(bamrepairedsortfn, "rb" )
    alignmentfile = pysam.AlignmentFile(bamrepairedsortfn, "rb" )
    outbam = pysam.Samfile(outhetfn, 'wb', template=samfile) 
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
    
 
    