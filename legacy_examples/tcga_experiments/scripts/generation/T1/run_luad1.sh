#!/bin/bash
#
#$ -cwd
module load samtools/1.2
module load bedtools
module load sambamba/0.5.4
module load vcftools

python ../../../../src/simulate.py -vcf ../../../../inputs/vcfs/normal_varscan.vcf -cnv_amp ../../../../inputs/beds/tcga/LUAD_1.gain.bed -cnv_del ../../../../inputs/beds/tcga/LUAD_1.loss.bed  -inbam ../../../../inputs/bams/Normal.bam -outbam luad1.bam -exons ../../../../inputs/beds/exons.bed -r ../../../../ref/hg19.fasta -cancertype luad -phase

