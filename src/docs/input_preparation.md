## Input data preparation

In this document, we will walk users through an example run of Bamgineer using data downloaded from the Personal Genome Project. Specifically, PGP participant huA2692E (https://my.pgp-hms.org/profile_public?hex=huA2692E). For your convenience, the VCF and bed files generated from this step are copied to examples/inputs. However, due to space limitations, the user is required to download the bam files for chromosome 21 from the PGP website (Steps 1-3 below).


### 1. Create a directory called splitbams; cd splitbams

### 2. Download chr21 bam file 

wget https://my.pgp-hms.org/user_file/download/2312 

### 3. Rename, index, and sort by name

module load samtools/1.2

mv 2312 chr21.bam \
samtools index chr21.bam \
samtools sort -n chr21.bam chr21.byname \

### 4. Create a directory called inputs; cd inputs

### 5. Create "exons.bed" file using whole genome coordinates of chr21 for hg19. The format should be tab-separated as follows: 

chromosome	chr_start_position	chr_end_position

For instance, navigate to inputs directory and create the following file and name it "exons.bed":

chr21   1       48129895

***Note: the current algorithm performs uses whole genome for phasing purposes but phased exome data can be run and the exons.bed file will be the bed file used to generate the bam. Further updates will be made regarding phasing and bamgineering exomes.

### 6. Download the VCF file

wget https://my.pgp-hms.org/user_file/download/2291 \
mv AE7EZ7WG3KG-EXT.vcf.gz genome.vcf.gz \
gunzip genome.vcf.gz

### 7. Filter VCF for chr21

bgzip genome.vcf 
tabix genome.vcf.gz
tabix genome.vcf.gz chr21 > chr21.vcf

### 8. Extract heterozygous SNPs from VCF

grep "0/1" chr21.vcf > chr21_het.vcf 

### 9. Remove indels

zcat genome.vcf.gz | grep "^#" > header.txt \
cat header.txt chr21_het.vcf > chr21_het_genome.vcf

vcftools --vcf chr21_het_genome.vcf --remove-indels --recode --recode-INFO-all --out snps_only \
mv snps_only.recode.vcf normal_het.vcf

### 10. Create an arbitrary bed file for desired CNVs. The format should be a tab-separated as follows:

chromosome	start_cnv_position	stop_cnv_position	allelic_ratio	absolute_copy_number

For instance, navigate to inputs directory and create the following file and name it "cnv.bed":

chr21	30227447	47076809	AAB	3

### 12. Create or edit config.cfg file

The config.cfg file includes paths to executables and references (see quickstart). A config file corresponding to the above generated data is created in inputs directory as an example.

### 13. You are now ready to run Bamgineer (see run_example1.sh in /scripts folder)
