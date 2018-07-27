## Input data preparation

In this document, we will walk users through an example run of Bamgineer on the data downloaded from Personal Genome Project
(https://my.pgp-hms.org/profile_public?hex=huA2692E) website. For convinience the VCF and Bed files generated from this step are copied
to examples/inputs. Howvever, due to space limitation the user is required to download the Bam files for chromosomes 21 and 22 from the
website (see Steps 1-3 below).


### 1. Create a folder called splitbams; cd splitbams

### 2. Download chr21 and chr22 bam files 

wget https://my.pgp-hms.org/user_file/download/2312 \
wget https://my.pgp-hms.org/user_file/download/2313

### 3. Rename , index , sort by name:

mv 2312 chr21.bam \
mv 2313 chr22.bam

samtools index chr21.bam \
samtools index chr22.bam

samtools sort -n chr21.bam chr21.byname \
samtools sort -n chr22.bam chr22.byname

### 4. Create a directory called inputs; cd inputs

### 5. Download vcf file (for whole genome)

*Note: The reason for downloading whole genome VCF instead of whole exome 
is due to the fact that different vcf formats exist and the genome version should be compatible

wget https://my.pgp-hms.org/user_file/download/2816 \
mv huA2692-veritas-gVCF-4.2.vcf.bz2 genome.vcf.bz2 \
bzip2 -d genome.vcf.bz2


### 6. Filter chr21 and chr22
bgzip genome.vcf \
tabix genome.vcf.gz

tabix genome.vcf.gz chr21 > chr21.vcf \
tabix genome.vcf.gz chr22 > chr22.vcf


### 7. Extract heterozygous positions (used in step 9)


grep "0/1" chr21.vcf > chr21_het.vcf \
grep "0/1" chr22.vcf > chr22_het.vcf


### 8. Create "exons.bed" file 

Since the bed file used to generate the Bam file was not available, we can create a hypothetical 
exon file by making a bed file from all contiguous regions (>50bp) above a coverage threshold (10X) as
explained here : https://www.biostars.org/p/86027

bedtools genomecov -ibam chr21.bam -bg > chr21.bedgraph \
awk '$4 > 10' chr21.bedgraph > chr21.gt50.bedgraph

bedtools genomecov -ibam chr22.bam -bg > chr22.bedgraph \
awk '$4 > 10' chr22.bedgraph > chr22.gt50.bedgraph

bedtools merge -i chr21.gt50.bedgraph  > exons_chr21.bed \
bedtools merge -i chr22.gt50.bedgraph > exons_chr22.bed

cat exons_chr21.bed exons_chr22.bed > exons.bed

### 9. Extract heterozygous SNPs (in step 7) that reside in the exons bed region, remove InDels

zcat genome.vcf.gz | grep "^#" > header.txt \
cat header.txt chr21_het.vcf chr22_het.vcf > chr21_22_het_genome.vcf

vcftools --vcf chr21_22_het_genome.vcf --bed exons.bed --out normal_het --recode \
vcftools --vcf normal_het.recode.vcf --remove-indels --recode --recode-INFO-all --out snps_only \
mv snps_only.recode.vcf normal_het.vcf

### 10. move normal_het.vcf and exons.bed to inputs directory

mv normal_het.vcf ../inputs \
mv exons.bed ../inputs

### 11. Create an arbitrary bed file for desired CNVs. The format for this file should be a tab separated(no header) file as follows:

chromosome	start_cnv_position	stop_cnv_pos	absolute_copy_number

For instance navigate to input directory and create the following file and name it cnv.bed

chr21	30227447	47076809	3

### 12. Create/edit config.cfg file

The config.cfg file includes configuration file including paths to executables and references (see quickstart). A config file corresponding to the above generated data is created in inputs directory for convenience.

**Please note that depending the specific runtime environment the paths may vary. Specifically if you are using HPC environment (such Sun Grid or Slurm), the paths can be set differently ("module load" command). A seperate documentation will be created for HPC cluster. 

### 13. You are now ready to run Bamgineer (Please see run_example1.sh in /scripts folder)
