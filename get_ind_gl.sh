#!/bin/sh
#Variables
REF_GENOME=./data/reference_genome/hs38DH.chr21.fa.gz

IND_LOC=./data/NA12878/
IND_S_1X=${IND_LOC}NA12878.chr21.s.1x

REGION=chr21:16600000-16800000

PANEL_NAME=./data/panel/1000GP.chr21.noNA12878.s
VCF=${PANEL_NAME}.sites.vcf.gz
TSV=${PANEL_NAME}.sites.tsv.gz

# Compute genotype likelihood based on the panel
bcftools mpileup -f ${REF_GENOME} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${IND_S_1X}.bam -Ou |
bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${IND_S_1X}.vcf.gz
bcftools index -f ${IND_S_1X}.vcf.gz