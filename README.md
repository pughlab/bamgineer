# Bamgineer
Introduces simulated allele-specific copy number variants into exome and targeted sequence data sets

## Author
Soroush Samadian

## Maintainer
Suluxan Mohanraj <suluxan.mohanraj@uhnresearch.ca>

## Description
Bamgineer is a tool that can be used to introduce user-defined haplotype-phased allele-specific copy number variations (CNV) into an existing Binary Alignment Mapping (BAM) file with demonstrated applicability to simulate somatic cancer CNVs in phased whole-genome sequencing datsets. This is done by introducing new read pairs sampled from existing reads, thereby retaining biases of the original data such as local coverage, strand bias, and insert size. As input, Bamgineer requires a BAM file and a list of non-overlapping genomic coordinates to introduce allele-specific gains and losses. We implemented parallelization of the Bamgineer algorithm for both standalone and high performance computing cluster environments, significantly improving the scalability of the algorithm. Bamgineer has been extensively tested on phased, whole-genome sequencing samples.

## Contact
If you have any questions with the package, please feel free to email Suluxan at <suluxan.mohanraj@uhnresearch.ca>.

## Running example Bamgineer workflow with Docker
### Please see bamgineer/docs/input_preparation for preparing your own files 

#### CHR21 bam files for NA12878 10X can be found here:
https://drive.google.com/file/d/1km9gupGi7W6aUE9XqBsiamGnwTpDrGpZ/view?usp=sharing

Tested with Docker version 17.05.0-ce, build 89658be

```sh
docker pull suluxan/bamgineer-v2
git clone https://github.com/pughlab/bamgineer.git
cd bamgineer/docker-example
# download google drive file (link above) and move into this directory
tar xjf splitbams.tar.bz2 
# start of bamgineer command
docker run --rm \
-v $(pwd):/src \
-it suluxan/bamgineer-v2 \
-config /src/inputs/config.cfg \
-splitbamdir src/splitbams \
-cnv_bed /src/inputs/cnv.bed \
-vcf src/inputs/normal_het.vcf \
-exons src/inputs/exons.bed \
-outbam tumour.bam \
-results src/outputs \
-cancertype LUAC1 
```

### Running without Docker:

        usage: python bamgineer/src/simulate.py 
      
        arguments:
           -config CONFIG_FILE,        configuration file including paths to tools/executables
           -splitdir SPLIT_BAMS_DIR,   input bam split by chromosomes
           -cnv_bed CNV_BED_FILE,      bed file containing non-overlapping CNVs following template
           -vcf FILTERED_VCF_FILE,     phased vcf file with indels removed 
           -exons COORDINATES_BED      bed file of whole genome coordinates or exons
           -outbam OUTPUT_BAM,         bam file name for output
           -results RESULTS_DIR,       output directory for final simulated bam results
           -cancertype CANCER_TYPE,    cancer type/acronym (OPTIONAL)


### Prerequisites - NOTE: many dependencies are outdated but will be continuously updated throughout the bamgineer development, docker image and docker workflow outlined above is ideal for ease-of-use

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
