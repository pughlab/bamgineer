#!/usr/bin/env python
import sys
import argparse
import time
from helpers import parameters as params
from methods import *

def main(args):
    configReader = params.GetConfigReader()
    params.InitConfigReader(args.configfile)
    outbamfn = args.outBamFile
    params.SetGainCNV(args.cnvAmpFile)
    params.SetCopyNumber(args.copyNumber)
    params.SetSplitBamsPath(args.splitbams)
    params.SetCancerType(args.cancerType)

    if(args.copyNumber <= 3):
        print('please use simulate module for copy number below 4')

    if (args.inbamFile):
        params.SetInputBam(args.inbamFile)

    results_path = configReader.get('RESULTS', 'results_path')

    # set software paths
    java_path, beagle_path, samtools_path, bedtools_path, vcftools_path, sambamba_path = params.GetSoftwarePath()
    params.SetSoftwarePath(java_path, beagle_path, samtools_path, bedtools_path, vcftools_path, sambamba_path)


    run_amp_pipeline(results_path)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='adds CN spikes to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-outbam', dest='outBamFile', required=True,
                        help='.bam file name for output')
    parser.add_argument('-cnv_amp', dest='cnvAmpFile', required=False,
                        help='CNV amplification .bed file name')
    parser.add_argument('-copynumber', dest='copyNumber', required=False,
                        help='The desired Copy number for CN >=4')
    parser.add_argument('-inbam', dest='inbamFile', required=False,
                        help='sam/bam file from which to obtain reads')
    parser.add_argument('-c', '--configFile', action='store', required=True, dest='configfile',
                        help='/path/to/config_file.cfg')
    parser.add_argument('-splitbamdir', dest='splitbams', required=False,
                        help='input bam split by chromosomes')
    parser.add_argument('-cancertype', dest='cancerType', required=False,
                        help='acronyms for cancer type')

    args = parser.parse_args()

    t0 = time.time()
    main(args)
    t1 = time.time()