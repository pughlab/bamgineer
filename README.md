# bamgineer
Introduces simulated allele-specific copy number variants into exome and targeted sequence data sets

Prerequisites
General NGS tools
Samtools1.2
Bedtools
VCFtools
Python packages
pysam (version 0.8.4): pysam
Note: the latest version of pysam (0.9.0) is not backward compatible with Samtools1.2
pyVCF pyvcf
pyBedTools pybedtools
Parameters
-inbam: input sorted and indexed normal bam file
-cnv_amp: bed file for allele-specific and cancer-specific CNV amplifications (null is not specified)
-cnv_del: bed file for allele-specific and cancer-specific CNV deletions (null is not specified)
-vcf: normal heterozygous vcf file (could be from HaplotypeCaller output, remember to filter indels and homozygous loci)
-exons: exon bed files (SureSelect V5 + UTR)
-r: reference hg19 fasta file (should be indexed: .fai, .amb, .ann, .pac, .awb)
-outbam: output engineered, sorted bam file
-phased(optional): Binary flag to perform phasing (BEAGLE) of SNPs prior to spiking CNV's
