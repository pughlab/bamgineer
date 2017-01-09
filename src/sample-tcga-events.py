#!/usr/bin/env python

import sys
import os
from random import random
import random
import argparse
from re import sub
from itertools import izip
import subprocess
from re import sub
import itertools
import numpy as np
from numpy.random import choice
import scipy.stats as stats
import ntpath


def median(lst):
    return np.median(np.array(lst))


def createDirectory(path):
    if not os.path.isdir(path):
        os.makedirs(path)
        
def runCommand(cmd):
    
    try:
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
       
    except Exception as e:
        print(e)

def findGISTICOverlap(segfilepath, gisticgainbedpath, gisticlossbedpath ):
    cmd = """awk '{if (($1!="Sample")&&($2 != "23") &&  ($6 < 0)){print $2"\t"$3"\t"$4"\t"$6}}' """ + segfilepath + ' > ' + 'loss.tmp' 
    cmd2 = """awk '{if (($1 != "Sample")&&($2 != "23") &&  ($6 > 0)){print $2"\t"$3"\t"$4"\t"$6}}' """ + segfilepath + ' > ' + 'gain.tmp' 
    runCommand(cmd)
    runCommand(cmd2)

    cmd3 = "bedtools intersect -a " + gisticgainbedpath + " -b gain.tmp -wa -u > gain.tmp2" 
    cmd4 = "bedtools intersect -a " + gisticlossbedpath + " -b loss.tmp  -wa -u > loss.tmp2" 
    runCommand(cmd3)
    runCommand(cmd4)
    
    nogainoverlaps= sum(1 for line in open('gain.tmp2' ))
    nolossoverlaps= sum(1 for line in open('loss.tmp2' ))
    
    return nogainoverlaps,nolossoverlaps


def mergeSegments(cnvsegfile,numgains, numlosses):
    path, filename = os.path.split(cnvsegfile)
    
    intermediatePath = "/".join([path, "intermediate"])
    cnvsegfn = "/".join([intermediatePath,sub('.seg$', '',filename)])
    finalCNVs = "/".join([path, "CNVs"])
    gisticpath = "/".join([path, "GISTIC"])
    fn = sub('.seg$', '',filename)
    gisticgainfn = "/".join([gisticpath,".".join([ntpath.basename(cnvsegfn),"GAIN.bed"])])
    gisticlossfn = "/".join([gisticpath,".".join([ntpath.basename(cnvsegfn),"LOSS.bed"])])
    sampe_cnt=0  
   
    gisticgainoverlapfn =  cnvsegfn+'_'+str(sampe_cnt) +'.gistic_gain.bed'
    gisticlossoverlapfn =  cnvsegfn+'_'+str(sampe_cnt) +'.gistic_loss.bed'
    
    finalbedfile =  "/".join([finalCNVs, ntpath.basename(cnvsegfn)+str(sampe_cnt) +'_final.bed'])
   
    tumorID  = "/".join([intermediatePath, sub('.seg$', '_IDs', ntpath.basename(cnvsegfile))])
    command= "cut -f 1 "+ cnvsegfile + " | uniq > " + tumorID
    runCommand(command)

    sampled_seg_path="/".join([intermediatePath,ntpath.basename(cnvsegfn)+"_sample"+ str(sampe_cnt)])+'.seg'
    seglength=0
    found= False
    gistic_gain=0
    gistic_loss=0
    
    tumor_scores=[]
    
    n_gains_written=0
    n_loss_written=0
    cntt=0
    tumorIDh=open(tumorID, 'r')
    for tumid in tumorIDh:
        #if(cntt < 10):
            #cntt +=1
        score=0
        
        sampled_seg_files = open(sampled_seg_path, 'w')
        idtmp = ''.join(str(e) for e in tumid.strip('\n').split("\t")  )
        sampled_seg_files.writelines([ line for line in open(cnvsegfile) if idtmp in line.strip('\n').split("\t")])
        sampled_seg_files.close()
        
        seglength= sum(1 for line in open(sampled_seg_path))
        gistic_gain,gistic_loss= findGISTICOverlap(sampled_seg_path,gisticgainfn,gisticlossfn)
        
        gainfn = "/".join([intermediatePath,".".join([ntpath.basename(cnvsegfn),str(sampe_cnt),"GAIN.bed"])]) 
        lossfn = "/".join([intermediatePath,".".join([ntpath.basename(cnvsegfn),str(sampe_cnt),"LOSS.bed"])])
        gainoverlapfn_gain_minus_loss = cnvsegfn+'_'+str(sampe_cnt) +'.gistic_gain_minus_loss.bed'
        lossoverlapfn_loss_minus_gain = cnvsegfn+'_'+str(sampe_cnt) +'.gistic_loss_minus_gain.bed'
              
        
        cmd = """awk '{if (($1!="Sample")&&($2 != "23") && ($4 -$3 > 2000000) && ($6 < -0.01)){print $2"\t"$3"\t"$4"\t"$6}}' """ + sampled_seg_path + ' > ' + lossfn 
        cmd2 = """awk '{if (($1 != "Sample")&&($2 != "23") && ($4 -$3 > 2000000) &&  ($6 > 0.01)){print $2"\t"$3"\t"$4"\t"$6}}' """ + sampled_seg_path + ' > ' + gainfn 
        
        cmd3= "bedtools merge -i " + lossfn + ' -d 3000000 > ' + lossfn+'.merged.bed' 
        cmd4= "bedtools merge -i " + gainfn +' -d 3000000 > '  + gainfn+'.merged.bed' 
        cmdsx = """awk '{if (($2 != "23") &&  ($6 < -0.01 || $6 > 0.01)){print $1"\t"$2"\t"$3"\t"$4"\t"$6}}' """ + sampled_seg_path + ' > ' + cnvsegfn+'.filtered.seg' 
        
        runCommand(cmd)
        runCommand(cmd2)
        runCommand(cmd3)
        runCommand(cmd4)
        runCommand(cmdsx)
        
        nl1 = sum(1 for line in open(gainfn+'.merged.bed'))
        nl2 = sum(1 for line in open(lossfn+'.merged.bed'))   
         
        cmd6 = "bedtools subtract -a " + gainfn+'.merged.bed' + " -b " + lossfn+'.merged.bed' + "  -A > " + lossoverlapfn_loss_minus_gain
        cmd7 = "bedtools subtract -a " + lossfn+'.merged.bed'+ " -b " + gainfn+'.merged.bed' + "  -A > " + gainoverlapfn_gain_minus_loss
                   
        runCommand(cmd6)
        runCommand(cmd7)
        n_gains_written = sum(1 for line in open(gainoverlapfn_gain_minus_loss))
        n_loss_written = sum(1 for line in open(lossoverlapfn_loss_minus_gain))
        
        score = abs(n_gains_written-numgains)+ abs(n_loss_written - numlosses) + abs(gistic_gain-numgains) +abs(gistic_loss - numlosses)
        tumor_scores.append((tumid,score) )
       
    tumor_scores=sorted(tumor_scores, key=lambda tup: tup[1])
    
    finalscores =  open("/".join([finalCNVs, ntpath.basename(cnvsegfn)+str(sampe_cnt) +'_final_scores.txt']), "w")
    for item in tumor_scores:
        finalscores.write("%s\n" % (item,))
    
    finalscores.close()
    
    cnt=0
    for score in tumor_scores[:3]:
        cnt +=1
        #tmid= score[0].strip('\n').split("\t")
    
        sampled_seg_path="/".join([finalCNVs,ntpath.basename(cnvsegfn)+"_sample"+ str(cnt)])+'.seg'
        sampled_seg_files = open(sampled_seg_path, 'w')
    
        idtmp = ''.join(str(e) for e in  score[0].strip('\n').split("\t")  )
        sampled_seg_files.writelines([ line for line in open(cnvsegfile) if idtmp in line.strip('\n').split("\t")])
        sampled_seg_files.close()
        
        seglength= sum(1 for line in open(sampled_seg_path))
        gistic_gain,gistic_loss= findGISTICOverlap(sampled_seg_path,gisticgainfn,gisticlossfn)
        
        gainfn = "/".join([finalCNVs,".".join([ntpath.basename(cnvsegfn),str(cnt),"GAIN.bed"])]) 
        lossfn = "/".join([finalCNVs,".".join([ntpath.basename(cnvsegfn),str(cnt),"LOSS.bed"])])
        gainoverlapfn_gain_minus_loss = cnvsegfn+'_'+str(cnt) +'.gain.bed'
        lossoverlapfn_loss_minus_gain = cnvsegfn+'_'+str(cnt) +'.loss.bed'
              
        
        cmd = """awk '{if (($1!="Sample")&&($2 != "23") && ($4 -$3 > 1000000) && ($6 < -0.01)){print $2"\t"$3"\t"$4"\t"$6}}' """ + sampled_seg_path + ' > ' + lossfn 
        cmd2 = """awk '{if (($1 != "Sample")&&($2 != "23") && ($4 -$3 > 1000000) &&  ($6 > 0.01)){print $2"\t"$3"\t"$4"\t"$6}}' """ + sampled_seg_path + ' > ' + gainfn 
        
        cmd3= "bedtools merge -i " + lossfn + ' -d 3000000 > ' + lossfn+'.merged.bed' 
        cmd4= "bedtools merge -i " + gainfn +' -d 3000000 > '  + gainfn+'.merged.bed' 
        cmdsx = """awk '{if (($2 != "23") &&  ($6 < -0.01 || $6 > 0.01)){print $1"\t"$2"\t"$3"\t"$4"\t"$6}}' """ + sampled_seg_path + ' > ' + cnvsegfn+'.filtered.seg' 
        
        runCommand(cmd)
        runCommand(cmd2)
        runCommand(cmd3)
        runCommand(cmd4)
        runCommand(cmdsx)
        
        nl1 = sum(1 for line in open(gainfn+'.merged.bed'))
        nl2 = sum(1 for line in open(lossfn+'.merged.bed'))   
         
        cmd6 = "bedtools subtract -a " + gainfn+'.merged.bed' + " -b " + lossfn+'.merged.bed' + "  -A > " + lossoverlapfn_loss_minus_gain
        cmd7 = "bedtools subtract -a " + lossfn+'.merged.bed'+ " -b " + gainfn+'.merged.bed' + "  -A > " + gainoverlapfn_gain_minus_loss    

        runCommand(cmd6)
        runCommand(cmd7)
    
    return tumor_scores
            
            

 
def findOverlapWithGistic(tumorgainbed, tumorlossbed, gisticgainbed, gisticlossbed, exp_gain_num, expected_loss_num):
    
    gistgain_overlp=sub('.gain$', 'gistic_overlap_gain', tumorgainbed)
    gistloss_overlp=sub('.gain$', 'gistic_overlap_gain', tumorgainbed)
    command1 = "bedtools intersect -a " + tumorgainbed + " -b " + gisticgainbed + " > " + gistgain_overlp + ' -wa'
    command2 = "bedtools intersect -a " + tumorlossbed + " -b " + gisticlossbed + " > " +  gistloss_overlp + ' -wa'
    runCommand( command1 )
    runCommand( command2 )

 
def random_line(afile):
    line = next(afile)
    for num, aline in enumerate(afile):
      if random.randrange(num + 2): continue
      line = aline
    return line 
        

def removeOverlap(bedfn1, bedfn2):
   
    bf1 = sub('.bed$', '.unique',bedfn1)
    bf2 = sub('.bed$', '.unique',bedfn2)
    
    command3 = "bedtools intersect -a " +bedfn1 + " -b "+ bedfn2 +" -v >" + bf1+'.bed'
    command4 = "bedtools intersect -a " +bedfn2  + " -b "+ bedfn1 +" -v >" + bf2+'.bed'
  
    runCommand( command3)
    runCommand( command4)
    
 
            
def mergeBedfiles(bedDir):
    for file in os.listdir(bedDir):
        
        if file.endswith("sorted2.bed"):
            path = os.path.abspath(bedDir)
            flp = "/".join([path, file])
            merged = "/".join([path, sub('.sorted2.bed$', '.merged', file)])
            
            command = "bedtools merge -i " +  flp + " -c 1 -o count > " +merged+'.bed'
            runCommand(command)

def main():
    
    
    Median_amp= 19600000
    Median_loss=900000
    
    BLCA_gain = 23
    BLCA_loss = 30
    BRCA_gain = 26
    BRCA_loss = 36
    CRC_gain = 24
    CRC_loss = 42
    GBM_gain = 26
    GBM_loss = 36
    HNSC_gain = 24
    HNSC_loss = 40
    KIRC_gain =  6
    KIRC_loss = 16
    LUAD_gain = 26
    LUAD_loss = 34
    LUSC_gain = 25
    LUSC_loss = 43
    OV_gain = 29
    OV_loss = 37
    UCEC_gain = 39
    UCEC_loss = 47
    
    for file in os.listdir('/mnt/work1/users/pughlab/projects/BAMgineer/inputs/TCGA-SEGS/'):
            
        if(file.endswith(".seg")):
            frawname = sub('.seg$', '', file)
            
            mergeSegments("/".join(['/mnt/work1/users/pughlab/projects/BAMgineer/inputs/TCGA-SEGS/',file]), eval(frawname+"_gain"), eval(frawname+'_loss'))
            
            
            
            #sample_events("/".join(['/mnt/work1/users/pughlab/projects/BAMgineer/TCGA-SEGS/',file]), eval(frawname+"_gain"), eval(frawname+'_loss'))
    
if __name__ == '__main__':
   main()
   
#def findUniqIDs(segfile, outseg):
#    ID =  bf1 = sub('.seg$', '_IDs',basename(segfile))
#    command= "cut -f 1 "+ segfile + " | uniq > " + "/".join([intermediatePath, ID])
#    runCommand(command)

