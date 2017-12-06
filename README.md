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

*Quick start: [Quick start](https://github.com/pughlab/bamgineer/blob/master/docs/quick_start)* 

## Google User Group (Q&A)
If you have any questions with the package, please feel free to post in our Google user group 
https://groups.google.com/d/forum/bamgineer or email us at bamgineer@googlegroups.com. We will try our best to reply as soon as 
possible.


## Using Bamgineer
Below is a general description of the file formats and prerequizites. For details on usage and example, please follow our detailed guide in [Quick start](https://github.com/pughlab/bamgineer/blob/master/docs/quick_start)

### Prerequisites

***General NGS tools*** 

Samtools1.2 \
Bedtools \
VCFtools

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
-r: reference hg19 fasta file (should be indexed: .fai, .amb, .ann, .pac, .awb) 
-phased(optional): Binary flag to perform phasing (BEAGLE) of SNPs prior to spiking CNV's


***Output***

-outbam: output engineered, sorted bam file 
