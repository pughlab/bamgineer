## Input data preparation

### Note: example files are already prepared - use "Running example Bamgineer workflow with Docker" in the README; this is for preparing your own files

In this document, we will walk users through an example run of Bamgineer using phased data downloaded from 10X Genomics. Specifically, sample NA12878 (https://support.10xgenomics.com/genome-exome/datasets/2.1.0/NA12878_WGS_210). 

### 1. Create a directory called splitbams and enter the directory

mkdir splitbams
cd splitbams

### 2. Download bam file 

wget http://s3-us-west-2.amazonaws.com/10x.files/samples/genome/NA12878_WGS_210/NA12878_WGS_210_phased_possorted_bam.bam

### 4. Download bam index file 

wget http://cf.10xgenomics.com/samples/genome/NA12878_WGS_210/NA12878_WGS_210_phased_possorted_bam.bam.bai

### 3. Filter for specific chromosome (example - smallest chromosome)

samtools view -b NA12878_WGS_210_phased_possorted_bam.bam chr21 > chr21.bam

### 4. Index filtered bam

samtools index chr21.bam

### 5. Sort reads by name

samtools sort -n chr21.bam chr21.byname   

### 4. Create a directory called inputs

cd ..
cd inputs

### 5. Create "exons.bed" file using whole genome coordinates of chr21 for hg19. The format should be tab-separated as follows: 

chromosome     chr_start_position     chr_end_position

For instance, navigate to inputs directory and create the following file and name it "exons.bed":

chr21     1     48129895

***Note: the current algorithm uses phased whole genome but phased exome data can be run and the exons.bed file in that case will be the bed file used to generate the bam. Further updates will be made regarding phasing and bamgineering exomes.

### 6. Download the VCF file

wget http://cf.10xgenomics.com/samples/genome/NA12878_WGS_210/NA12878_WGS_210_phased_variants.vcf.gz

### 7. Download the VCF index file

wget http://cf.10xgenomics.com/samples/genome/NA12878_WGS_210/NA12878_WGS_210_phased_variants.vcf.gz.tbi

### 8. Filter VCF for chr21 (example - smallest chromosome)

tabix NA12878_WGS_210_phased_variants.vcf.gz chr21 > chr21.vcf

### 9. Remove indels

zcat NA12878_WGS_210_phased_variants.vcf.gz | grep "^#" > header.txt  
cat header.txt chr21.vcf > chr21_header.vcf

vcftools --vcf chr21_header.vcf --remove-indels --recode --recode-INFO-all --out snps_only  
mv snps_only.recode.vcf normal_het.vcf

### 10. Create a bed file for desired CNVs. The format should be tab-separated as follows:

chromosome     start_cnv_position     stop_cnv_position     allelic_ratio     absolute_copy_number

For instance, navigate to inputs directory and create the following tab-separated file and name it "cnv.bed" (can use vim or other text editor):

chr21     30227447     30827447     AAB     3

### 12. Create or edit config.cfg file

The config.cfg file includes paths to executables and references (see quickstart). A config file corresponding to the above generated data is created in inputs directory as an example.

### 13. You are now ready to run Bamgineer (see run1.sh in /docker-example/scripts folder)