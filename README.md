# bamgineer
Introduces simulated allele-specific copy number variants into exome and targeted sequence data sets

## Author
Soroush Samadian


## Maintainer
Soroush Samadian <soroush.samadian@uhnresearch.ca>

## Description
Bamgineer is a tool to introduce user-defined haplotype-phased allele-specific copy number events into an existing Binary Alignment Mapping (BAM) file, with a focus on targeted and exome sequencing experiments. As input, this tool requires a read alignment file (BAM format), lists of non-overlapping genome coordinates for introduction of gains and losses (bed file), and an optional file defining known haplotypes (vcf format). 


## Prerequisites



## General NGS tools 

Samtools1.2

Bedtools

VCFtools


***Python packages***

*pysam (version 0.8.4): [pysam](https://pypi.python.org/pypi/pysam)*

Note: the latest version of pysam (0.9.0) is not backward compatible with Samtools1.2

*pyVCF [pyvcf](https://pypi.python.org/pypi/PyVCF)*

*pyBedTools [pybedtools](https://pypi.python.org/pypi/pybedtools)*


*** Parameters ***

-inbam: input sorted and indexed normal bam file 

-cnv_amp: bed file for allele-specific and cancer-specific CNV amplifications (null is not specified)

-cnv_del: bed file for allele-specific and cancer-specific CNV deletions (null is not specified)

-vcf: normal heterozygous vcf file (could be from HaplotypeCaller output, remember to filter indels and homozygous loci)

-exons: exon bed files (SureSelect V5 + UTR)

-r: reference hg19 fasta file (should be indexed: .fai, .amb, .ann, .pac, .awb)

-outbam: output engineered, sorted bam file 

-phased(optional): Binary flag to perform phasing (BEAGLE) of SNPs prior to spiking CNV's
