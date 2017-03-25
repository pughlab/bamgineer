import sys
import pysam
import subprocess
from re import sub
import os
import numpy as np
import utils 

mergedsortfn = sys.argv[1]
origroibamfn = sys.argv[2]
GAIN_FINAL= sys.argv[3]


ratio = float(utils.countReads(mergedsortfn))/float(utils.countReads(origroibamfn))
samplerate= round(0.5/(ratio*0.98),2)



success = False
if(samplerate < 1.0):
    utils.subsample(mergedsortfn, GAIN_FINAL,str(samplerate))
    success = True
elif(samplerate< 1.5):
    print('sample rate is larger than!  ' + str(samplerate))
    os.rename(mergedsortfn, GAIN_FINAL)
    success = True
else:
    print("not enough number of reads found for "+ mergedsortfn)