#!/bin/bash
#
#$ -cwd

PDIR='/bamgineer/docker-example'

python /bamgineer/src/simulate.py \
-cnv_bed ${PDIR}/inputs/cnv.bed \
-config ${PDIR}/inputs/config.cfg \
-splitbamdir ${PDIR}/splitbams \
-outbam tumor.bam \
-cancertype LUAC1 \
-results ${PDIR}/outputs \
-vcf ${PDIR}/inputs/normal_het.vcf \
-exons ${PDIR}/inputs/exons.bed \
