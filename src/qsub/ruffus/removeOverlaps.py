import os
import sys
import pysam
import pybedtools
from re import sub
import shutil
import subprocess
#from helpers import pipeline

global bases

nonhetfn = sys.argv[1]
hetaltfn=  sys.argv[2]
nhunique = sys.argv[3]
hetuniq = sys.argv[4]
path =  sys.argv[5]

def bamDiff(bamfn1, bamfn2,path):
    command = " ".join(["bam diff", "--in1", bamfn1, "--in2", bamfn2, "--out" ,"/".join([path,"diff.bam"])]) # ("roi.bam - vcf.bam"; seperate reads that do not overlap SNP regions from the ones that do)
    subprocess.check_output(command, shell = True)

path1, filename1 = os.path.split(nonhetfn)
path2, filename2 = os.path.split(hetaltfn)

#alignmentfile = pysam.AlignmentFile(bam1sorted+'.bam', "rb" )
#ovbam = pysam.Samfile(overlapbamfn, 'wb', template=alignmentfile )

bamDiff(nonhetfn, hetaltfn , path )
if(os.path.isfile("/".join([path, 'diff_only1_'+ filename1]))):
    os.rename("/".join([path, 'diff_only1_'+ filename1]), nhunique )
else:
    os.rename("/".join([path, filename1]), nhunique )
if(os.path.isfile("/".join([path, 'diff_only2_'+ filename2]))):
    os.rename("/".join([path, 'diff_only2_'+ filename2]), hetuniq )
else:
    os.rename("/".join([path, filename2]), hetuniq )



#bedfile = open(bedfn, 'r')
#line=0
#for bedline in bedfile:
#    c = bedline.strip().split()
#    
#    if (len(c) == 3 ):
#        chr2 = c[0]
#        chr = c[0].strip("chr")
#        start = int(c[1])
#        end   = int(c[2])
#    else :
#        #print(len(c))
#        continue
#    
#    try:
#        readmappings = alignmentfile.fetch(chr2, start, end)
#    except  ValueError as e:
#        print("problem fetching the read ")
#    
#    
#    for shortread in readmappings:
#        try:
#            ovbam.write(shortread)
#            line+=1
#        except ValueError as e:
#            print ("problem removing read :" + shortread.qname)
#
#ovbam.close()
         

#if line > 0:
#    overlapsorted =  sub('.bam$','.sorted',overlapbamfn)   
#    pysam.sort(overlapbamfn, overlapsorted)
#    
#    bamDiff(bam1sorted+'.bam', overlapsorted +'.bam', path )
#    os.rename("/".join([path, 'diff_only1_'+ bam1sorted+'.bam']), outputfn )
#    print('renaming ' + "/".join([path, 'diff_only1_'+ bam1sorted+'.bam']) + '  to  ' + outputfn)
#    os.remove(overlapbamfn)
#else:
#    os.rename("/".join([path, bam1sorted+'.bam']), outputfn )
#    os.remove(overlapbamfn)
#    print ('removing '+overlapbamfn )



