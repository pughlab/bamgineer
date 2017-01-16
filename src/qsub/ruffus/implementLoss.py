import os
import sys
import pysam
import subprocess
from re import sub

chr = sys.argv[1]
event = sys.argv[2]
splittmpbams_hap = sys.argv[3]
finalbams_hap = sys.argv[4]
split_path = sys.argv[5]
sortbyCoord = split_path + '/' + chr.lower() + '.bam'

refbam = "/".join([splittmpbams_hap, chr +'_'+ event.lower() +"_REFREADS.bam"])
altbam = "/".join([splittmpbams_hap, chr +'_'+event.lower() +"_ALTREADS.bam"])
        
refbamsorted = sub('.bam$','.sorted', refbam)
altbamsorted = sub('.bam$','.sorted', altbam)

if(os.path.isfile(refbam)):
    pysam.sort(refbam,refbamsorted)
    #os.remove(refbam)

inbam_deletion = "/".join([finalbams_hap, str(chr).upper() + 'LOSS_FINAL.bam']) 
NONHET = "/".join([splittmpbams_hap,  (chr + '_'+event).upper()  +'_NONHET.bam'])
HET_ALT =  "/".join([splittmpbams_hap, (chr +'_'+ event).upper() +'_HET_ALT.bam'])

def subsample(bamfn1, bamfn2, samplingrate = 0.5):
    command = " ".join(["samtools view -s", samplingrate ,"-b", bamfn1, ">", bamfn2])
    subprocess.check_output(command, shell = True)
def bamDiff(bamfn1, bamfn2,path):
    command = " ".join(["bam diff", "--in1", bamfn1, "--in2", bamfn2, "--out" ,"/".join([path,"diff.bam"])]) # ("roi.bam - vcf.bam"; seperate reads that do not overlap SNP regions from the ones that do)
    subprocess.check_output(command, shell = True)


if(not os.path.isfile(NONHET) and not os.path.isfile(HET_ALT)): 
    os.symlink(sortbyCoord, inbam_deletion)
    
else:
    
    if(os.path.isfile(NONHET) and os.path.isfile(HET_ALT)):
        loss_NONHET_SAMPLED =   "/".join([splittmpbams_hap,  chr + event +'_NONHET_SAMPLED.bam'])
        subsample(NONHET ,loss_NONHET_SAMPLED , str(0.5))
        inbam_minus_loss_NONHET_ALT = "/".join([splittmpbams_hap, chr + event+ '_MINUS_NONHET_ALT.bam'])
        
        bamDiff(sortbyCoord , loss_NONHET_SAMPLED, splittmpbams_hap)
        os.rename("/".join([splittmpbams_hap,  'diff_only1_' + chr + '.bam']), inbam_minus_loss_NONHET_ALT)
        
        
        bamDiff(inbam_minus_loss_NONHET_ALT, refbamsorted+'.bam', splittmpbams_hap)
        os.rename("/".join([splittmpbams_hap,'diff_only1_'+ chr + event+ '_MINUS_NONHET_ALT.bam']),inbam_deletion)
