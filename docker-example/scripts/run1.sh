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
-phase
