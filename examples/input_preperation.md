## Input data preperation

In this document, we will walk users through an example run of Bamgineer on the data downlowded from Personal Genome Project
(https://my.pgp-hms.org/profile_public?hex=huA2692E) website. For convinience the VCF and Bed files generated from this step are copied
to examples/inputs. Howvever, due to space limitation the user is required to download the Bam files for chromosomes 21 and 22 from the
website (see Steps 1-3 below).


### 1. Create a foldr called splitbams; cd splitbams

### 2. Download chr21 and chr22 bam files 

wget https://my.pgp-hms.org/user_file/download/2312
wget https://my.pgp-hms.org/user_file/download/2313

### 3. Rename , index , sort by name:

mv 2312 chr21.bam
mv 2313 chr22.bam

samtools index chr21.bam
samtools index chr22.bam

samtools sort -n chr21.bam chr21.byname
samtools sort -n chr22.bam chr22.byname

### 4. Create a directory called inputs; cd inputs

### 5. Download vcf file (for whole genome)

*Note: The reason for downloading whole genome VCF instead of whole exome 
is due to the fact that different vcf formats exist and the genome version is compatible

wget https://my.pgp-hms.org/user_file/download/2816
mv huA2692-veritas-gVCF-4.2.vcf.bz2 genome.vcf.bz2
bzip2 -d genome.vcf.bz2


### 6. Filter chr21 and chr22
bgzip genome.vcf
tabix genome.vcf.gz

tabix genome.vcf.gz chr21 > chr21.vcf
tabix genome.vcf.gz chr21 > chr22.vcf


### 7. Extract heterozygous positions


grep "0/1" chr21.vcf > chr21_het.vcf
grep "0/1" chr22.vcf > chr22_het.vcf

cat chr21_het.vcf chr22_het.vcf > chr21_22_het_genome.vcf

### 8. Create "exons.bed" file

Since the bed file used to generate the Bam file was not available, we can create a hypothetical 
exon file by making a bed file from all contiguous regions (>50bp) above a coverage threshold (50X) as
explained here : https://www.biostars.org/p/86027

bedtools genomecov -ibam chr21.bam -bg > chr21.bedgraph
awk '$4 > 50' chr21.bedgraph > chr21.gt50.bedgraph

bedtools genomecov -ibam chr22.bam -bg > chr22.bedgraph
awk '$4 > 50' chr22.bedgraph > chr22.gt50.bedgraph

bedtools merge -i chr21.gt50.bedgraph  > exons_chr21.bed
bedtools merge -i chr22.gt50.bedgraphh > exons_chr22.bed

cat exons_chr21.bed exons_chr22.bed > exons.bed

### 9. Extract heterozygous SNPs (in step 7) that reside in the exons bed region


vcftools --vcf chr21_22_het_genome.vcf --bed exons.bed --out normal_het --recode

mv normal_het.recode.vcf normal_het.vcf

### 10. move normal_het.vcf and exons.bed to inputs dir

mv normal_het.vcf ../inputs
mv exons.bed ../inputs

### 11. Create an arbitrary bed file for desired CNVs. The format for this file should be a tab separated(no header) file as follows:

chromosome	start_cnv_position	stop_cnv_pos	absolute_copy_number

For instance navigate to input directory and create the following file and name it cnv.bed

chr21	30227447	47076809	3

### 12. You are now ready to run Bamgineer (Please see run_example)
