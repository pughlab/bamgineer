#!/bin/bash
#
#$ -cwd

PDIR='/mnt/work1/users/home2/tpcoop1/pughlab/projects/BAMgineer/bamgineer-v2/test'

python /mnt/work1/users/home2/tpcoop1/git/bamgineer/src/simulate.py \
-cnv_bed ${PDIR}/inputs/cnv1.bed \
-config ${PDIR}/inputs/config1.cfg \
-splitbamdir ${PDIR}/splitbams \
-outbam CN1.bam \
-cancertype LUAC1 \
-phase
