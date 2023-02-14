#!/bin/sh
#Variables
REF_GENOME=./data/reference_genome/hs38DH.chr21.fa.gz

IND_LOC=./data/NA12878/
IND_NAME=${IND_LOC}NA12878
IND_S=${IND_NAME}.chr21.s
IND_S_1X=${IND_S}.1x

REGION=chr21:16600000-16800000

PANEL_NAME=./data/panel/1000GP.chr21.noNA12878.s
VCF=${PANEL_NAME}.sites.vcf.gz
TSV=${PANEL_NAME}.sites.tsv.gz

MAP_SNP=./data/affi/snp6.map

# Compute genotype likelihood based on the panel
bcftools mpileup -f ${REF_GENOME} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${IND_S_1X}.bam -Ou |
bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${IND_S_1X}.vcf.gz
bcftools index -f ${IND_S_1X}.vcf.gz

# Get individual SNP
plink2 --chr-set 22 --bcf ${IND_S}.bcf \
    --set-missing-var-ids @:# \
    --max-alleles 2 \
    --output-chr 'chr26' \
    --snps-only --allow-extra-chr \
    --extract ${MAP_SNP} \
    --chr 1-38,X --recode vcf bgz \
    --out ${IND_S}.snp

bcftools view -e 'GT="./."||GT="."' ${IND_S}.snp.vcf.gz -Oz -o ${IND_S}.snp.filtered.vcf.gz
bcftools index -f ${IND_S}.snp.filtered.vcf.gz

# Phase the SNP data
SHAPEIT5_phase_common -I ${IND_S}.snp.filtered.vcf.gz -H ${PANEL_NAME}.phased.vcf.gz -O ${IND_S}.phased.vcf -R ${REGION}