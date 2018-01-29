#!/bin/bash
#
#$ -cwd

python ../../src/simulate.py \
-cnv_bed ../inputs/cnv.bed \
-c ../inputs/config.cfg \
-splitbamdir ../splitbams \
-outbam test.bam \
-phase
