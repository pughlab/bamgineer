#!/bin/bash
#
#$ -cwd
module load samtools/1.2
module load bedtools
module load sambamba/0.5.4
module load vcftools

python /mnt/work1/users/pughlab/projects/BAMgineer/src/simulate.py -vcf /mnt/work1/users/pughlab/projects/BAMgineer/inputs/vcfs/normal_varscan.vcf -cnv_amp /mnt/work1/users/pughlab/projects/BAMgineer/inputs/beds/tcga/OV_3.gain.bed -cnv_del /mnt/work1/users/pughlab/projects/BAMgineer/inputs/beds/tcga/OV_3.loss.bed  -inbam /mnt/work1/users/pughlab/projects/BAMgineer/inputs/bams/Normal.bam -outbam ov3.bam -exons /mnt/work1/users/pughlab/projects/BAMgineer/inputs/beds/exons.bed -r /mnt/work1/users/pughlab/projects/BAMgineer/ref/hg19.fasta -cancertype ov -phase

