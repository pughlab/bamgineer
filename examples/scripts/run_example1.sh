#!/bin/bash
#
#$ -cwd

python ../../scr/simulate.py \
-cnv_bed ../inputs/cnv.bed \
-c ../inputs/config.bed \
-splitbamdir ../splitbams \
-outbam test.bam \
-phase
