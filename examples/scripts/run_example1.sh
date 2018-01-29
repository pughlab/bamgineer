#!/bin/bash
#
#$ -cwd

python ../../src/simulate.py \
-cnv_bed ../inputs/cnv.bed \
-config ../inputs/config.cfg \
-splitbamdir ../splitbams \
-outbam test.bam \
-phase
