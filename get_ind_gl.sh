#!/bin/sh
conda init bash
conda activate env_tools
#Variables
CHR=21
REF_GENOME=./data/reference_genome/hs38DH.chr${CHR}.fa.gz

IND_LOC=./data/NA12878/
IND_NAME=${IND_LOC}${CHR}/NA12878
IND_S=${IND_NAME}.chr${CHR}.s
IND_S_1X=${IND_S}.1x

REGION=chr${CHR}:16600000-16800000

PANEL_NAME=./data/panel/${CHR}/1000GP.chr${CHR}.noNA12878.s
VCF=${PANEL_NAME}.sites.vcf.gz
TSV=${PANEL_NAME}.sites.tsv.gz

# Compute genotype likelihood based on the panel
bcftools mpileup -f ${REF_GENOME} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${IND_S_1X}.bam -Ou |
bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${IND_S_1X}.vcf.gz
bcftools index -f ${IND_S_1X}.vcf.gz

# Phase the SNP data
SHAPEIT5_phase_common -I ${IND_S}.snp.filtered.vcf.gz -H ${PANEL_NAME}.phased.vcf.gz -O ${IND_S}.phased.vcf -R ${REGION}