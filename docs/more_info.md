 
## More Information for running Bamgineer

NOTE: many dependencies are outdated but will be continuously updated throughout the bamgineer development, docker image and docker workflow outlined above is ideal for ease-of-use

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

## Running Bamgineer

Users may split the input bam file by chromosome (e.g. chr1.bam, chr2.bam) sorted by coordinates along with index files (chr.bam.bai, chr1.bam.bai, etc.) and also have them sorted by name using "by.name.bam" extension(chr1.byname.bam, chr2.byname.bam, etc) prior to running the program. Please see bamgineer/docs/input_preparation for more information. The directory containing these bam is given to bamgineer as an argument ( -splitbamdir). If no such directory is given Bamgineer performs the splitting and sorting by coordinates and name as a preprocessing step (use -inbam argument). Please note that it is not necessary to have the data for all the chromosomes. For example, a user can run bamgineer on a bam file that is only comprised of chromosomes 1, 3 and 22.

###  Editing the Config file

Config file includes path to executable software and references. An example configuration file can be found in "docker-example/inputs/config.cfg". Template:

**[SOFTWARE]**

This category includes the path to external tools(such as Java, Beagle, etc). If you need to use a different version of and existing file (e.g. Java 1.8 instead of Java 1.7), you can do so by editing the path. In HPC cluster environment these path can be loaded using "module load" command (e.g. module load java) and then finding path to executable using "which" command (e.g. which java).