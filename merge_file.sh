#!/bin/sh
conda init bash
conda activate env_tools

# Region of interest
CHR=22
REGION=chr${CHR}:16600000-16800000

# Concatenate the panel files

PANEL_CHR21=./data/panel/21/1000GP.chr21.noNA12878.s
PANEL_CHR22=./data/panel/22/1000GP.chr22.noNA12878.s
PANEL_S=./data/panel/21_22/1000GP.chr21_22.noNA12878.s

# Concatenate the panel files
echo 'Concatenate the panel files'
mkdir -p data/panel/21_22
bcftools concat -Ob -o ${PANEL_S}.bcf ${PANEL_CHR21}.bcf ${PANEL_CHR22}.bcf
bcftools index -f ${PANEL_S}.bcf --threads 4

# Select only the SNPS and drop Genotypes
echo 'Select only SNPs and drop Genotypes'
bcftools view -G -m 2 -M 2 -v snps ${PANEL_S}.bcf -Oz -o ${PANEL_S}.sites.vcf.gz
bcftools index -f ${PANEL_S}.sites.vcf.gz

# Convert to TSV
echo 'Convert to TSV'
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${PANEL_S}.sites.vcf.gz | bgzip -c > ${PANEL_S}.sites.tsv.gz
tabix -s1 -b2 -e2 ${PANEL_S}.sites.tsv.gz

# Merge fasta files
echo 'Merge fasta files'
REF_FASTA=./data/reference_genome/hs38DH.chr
REF_FASTA_21=${REF_FASTA}21.s.fa
REF_FASTA_22=${REF_FASTA}22.s.fa

cat ${REF_FASTA_21} ${REF_FASTA_22} > ${REF_FASTA}21_22.s.fa
samtools faidx ${REF_FASTA}21_22.s.fa

# Merge SNP files
echo 'Merge SNP files'
SNP_FILE=./data/affi/snp6
SNP_FILE_21=${SNP_FILE}.chr21.s.map
SNP_FILE_22=${SNP_FILE}.chr22.s.map

cat ${SNP_FILE_21} ${SNP_FILE_22} > ${SNP_FILE}.chr21_22.s.map

# Merge individuals file
echo 'Merge individuals file'
IND_FILE=./data/NA12878/
IND_FILE_21=${IND_FILE}/21/NA12878.chr21.s
IND_FILE_22=${IND_FILE}/22/NA12878.chr22.s

mkdir -p ${IND_FILE}/21_22
# Merge BAM file
samtools merge -r ${IND_FILE}/21_22/NA12878.21_22.s.bam ${IND_FILE_21}.1x.bam ${IND_FILE_22}.1x.bam
samtools index ${IND_FILE}/21_22/NA12878.21_22.s.bam

# Merge GL file
bcftools concat -Oz -o ${IND_FILE}/21_22/NA12878.21_22.s.1x.vcf.gz ${IND_FILE_21}.1x.vcf.gz ${IND_FILE_22}.1x.vcf.gz
bcftools index -f ${IND_FILE}/21_22/NA12878.21_22.s.1x.vcf.gz --threads 4

# Merge truth bcf
bcftools concat -Oz -o ${IND_FILE}/21_22/NA12878.21_22.s.bcf ${IND_FILE_21}.bcf ${IND_FILE_22}.bcf
bcftools index -f ${IND_FILE}/21_22/NA12878.21_22.s.bcf --threads 4

