
        usage: python simulate.py [-h] [--inbam INPUT_BAM] [--outbam OUTPUT_BAM]
                         [--cnvgain CNV_GAIN_FILE] [--cnvloss CNV_LOSS_FILE] [--splitdir BAM_SPLIT_DIR] 
                         [--config CONFIG_FILE] [--cancer CANCER_TYPE][--phase PHASE] [--ctdna CT_DNA]
                         
        
        arguments:
           --inbam  INPUT_BAM ,         bam file from for input
           --outbam OUTPUT_BAM,         bam file name for output
           --config CONFIG_FILE,        configuration file including paths to executables and references
           --cnvgain CNV_GAIN_FILE,     bed file name containing non-overlapping gain regions
           --cnvloss CNV_LOSS_FILE,     bed file name containing non-overlapping loss regions
        
        optional arguments:
           --h, --help                  show this help message and exit
           --cancer CANCER_TYPE,        cancer type/acronym
           --splitdir BAM_SPLIT_DIR,    input bam split by chromosomes
           --p PHASE,                   whether SNP phasing should be applied
           --ctdna CT_DNA,              whether simulation is on reads obtained from ctDNA sequencing data   
