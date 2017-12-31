# Bamgineer
Introduces simulated allele-specific copy number variants into exome and targeted sequence data sets

## Author
Soroush Samadian


## Maintainer
Soroush Samadian <soroush.samadian@uhnresearch.ca>

## Description
Bamgineer is a tool that can be used to introduce user-defined haplotype-phased allele-specific copy number variations (CNV) into an 
existing Binary Alignment Mapping (BAM) file and demonstrated applicability to simulate somatic cancer CNVs in exome and targeted cell-
free DNA sequencing data sets. This is done by introducing new read pairs sampled from existing reads, thereby retaining biases of the 
original data such as local coverage, strand bias, and insert size. As input, Bamgineer requires a BAM file and two lists of non-
overlapping genomic coordinates to introduce allele-specific gains and losses. The user may explicitly provide known haplotypes or chose 
to use the BEAGLE phasing module that is already incorporated within Bamgineer. We implemented parallelization of the Bamgineer 
algorithm for both standalone and high performance computing cluster environments, significantly improving the scalability of the 
algorithm.Bamgineer has been extensively tested on whole exome sequencing and targeted gene panel applied to cell-free DNA sequencing 
data.

*Quick start: [Quick start](https://github.com/pughlab/bamgineer/blob/master/docs/quick_start.md)* 

## Google User Group (Q&A)
If you have any questions with the package, please feel free to post in our Google user group 
https://groups.google.com/d/forum/bamgineer or email us at bamgineer@googlegroups.com. We will try our best to reply as soon as 
possible.

## Using Bamgineer
Below is a general description of the file formats and prerequizites. For details on usage and example, please follow our detailed guide in [Quick start](https://github.com/pughlab/bamgineer/blob/master/docs/quick_start.md)

### Running options:

        usage: python simulate.py [-h] [--inbam INPUT_BAM] [--outbam OUTPUT_BAM]
                         [--cnvgain CNV_LIST] [--splitdir BAM_SPLIT_DIR] 
                         [--config CONFIG_FILE] [--cancer CANCER_TYPE][--phase PHASE] [--ctdna CT_DNA]
                         
        
        arguments:
           --inbam  INPUT_BAM ,         bam file from for input
           --outbam OUTPUT_BAM,         bam file name for output
           --config CONFIG_FILE,        configuration file including paths to executables and references
           --cnvgain CNV_LIST,     bed file name containing non-overlapping CNV regions
        
        optional arguments:
           --h, --help                  show this help message and exit
           --cancer CANCER_TYPE,        cancer type/acronym
           --splitdir BAM_SPLIT_DIR,    input bam split by chromosomes
           --p PHASE,                   whether SNP phasing should be applied
           --chr_list                   list of chromosomes to process (default: all)
           --ctdna CT_DNA,              whether simulation is on reads obtained from ctDNA sequencing data   


### Prerequisites

***General NGS tools*** 

*Samtools (version 1.2): [samtools](http://samtools.sourceforge.net)* \
*Bedtools:[bedtools](http://bedtools.readthedocs.io/en/latest/)*\
*VCFtools:[vcftools](http://vcftools.sourceforge.net/index.html)*\
*BamUtil:[bamutil](https://genome.sph.umich.edu/wiki/BamUtil)*

***Python packages***

*pysam (version 0.8.4): [pysam](https://pypi.python.org/pypi/pysam)* \
Note: the latest version of pysam (0.9.0) is not backward compatible with Samtools1.2 \
*pyVCF [pyvcf](https://pypi.python.org/pypi/PyVCF)* \
*pyBedTools [pybedtools](https://pypi.python.org/pypi/pybedtools)*


***Input parameters***

-inbam: input sorted and indexed normal bam file \
-cnv_list: path to the directory contaiing bed files for allele-specific and cancer-specific CNV events. Each bed file is a tab seperated file (no header) containing the following: 

chr     start   stop    absolute_copy_number 

-vcf: normal heterozygous vcf file (could be from HaplotypeCaller output, not including indels and homozygous loci) \
-target_region: bed file containing the target regions (exons or any user-specified region) \ 
-phased(optional): Binary flag to perform phasing (BEAGLE) of SNPs prior to spiking CNV's


***Output***

-outbam: output engineered, sorted bam file 
