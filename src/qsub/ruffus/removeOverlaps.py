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

bamDiff(nonhetfn, hetaltfn , path )
if(os.path.isfile("/".join([path, 'diff_only1_'+ filename1]))):
    os.rename("/".join([path, 'diff_only1_'+ filename1]), nhunique )
else:
    os.rename("/".join([path, filename1]), nhunique )
if(os.path.isfile("/".join([path, 'diff_only2_'+ filename2]))):
    os.rename("/".join([path, 'diff_only2_'+ filename2]), hetuniq )
else:
    os.rename("/".join([path, filename2]), hetuniq )

