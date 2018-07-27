#!/bin/bash
#
#$ -cwd
module load samtools/1.2
module load bedtools
module load sambamba/0.5.4
module load vcftools

python ../../../../src/simulate.py -vcf ../../../../inputs/vcfs/normal_varscan.vcf -cnv_amp ../../../../inputs/beds/tcga/GBM_3.gain.bed -cnv_del ../../../../inputs/beds/tcga/GBM_3.loss.bed  -inbam ../../../../inputs/bams/Normal.bam -outbam gbm3.bam -exons ../../../../inputs/beds/exons.bed -r ../../../../ref/hg19.fasta -cancertype gbm -phase

