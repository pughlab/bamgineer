 
## Using Bamgineer

### Prerequisites

***General NGS tools*** 

*Samtools (version 1.2): [samtools](http://samtools.sourceforge.net)* \
*Bedtools:[bedtools](http://bedtools.readthedocs.io/en/latest/)*\
*VCFtools:[vcftools](http://vcftools.sourceforge.net/index.html)*

***Python packages***

*pysam (version 0.8.4): [pysam](https://pypi.python.org/pypi/pysam)* \
Note: the latest version of pysam (0.9.0) is not backward compatible with Samtools1.2 \
*pyVCF [pyvcf](https://pypi.python.org/pypi/PyVCF)* \
*pyBedTools [pybedtools](https://pypi.python.org/pypi/pybedtools)*


***Input parameters***

-inbam: input sorted and indexed normal bam file \
-cnv_amp: bed file for allele-specific and cancer-specific CNV amplifications (null is not specified) \
-cnv_del: bed file for allele-specific and cancer-specific CNV deletions (null is not specified) \
-vcf: normal heterozygous vcf file (could be from HaplotypeCaller output, not including indels and homozygous loci) \
-target_region: bed file containing the target regions (exons or any user-specified region) \
-splitdir: input bam split by chromosomes
-phased(optional): Binary flag to perform phasing (BEAGLE) of SNPs prior to spiking CNV's

It is recommend

***Output***

-outbam: output engineered, sorted bam file

## Running Bamgineer

It is recommended that users split the input bam file by chromosome (e.g. chr1.bam, chr2.bam) sorted by coordinates along with index 
files(chr.bam.bai, chr1.bam.bai ,etc) and also have them sorted by name using "by.name.bam" extension(chr1.byname.bam, chr2.byname.bam, 
etc) prior to running the program. The directory containing these bam is given to bamgineer as an argument ( --splitdir). If no such 
directory is given Bamgineer performs the splitting and sorting by coordinates and name as a preprocessing step. Please note that it is 
not necessary to have the data for all the chromosomes. For instance, a user can run bamgineer on a bam file that is only comprised of 
chromosomes 1,3 and 22.

###  Editting the Config file

### Running on single node



### Running on HPC (Sun Grid Engine)


### Running options:
    
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
