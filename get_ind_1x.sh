#!/bin/sh
conda init bash
conda activate env_tools

CHR=21
REF_GENOME=./data/reference_genome/hs38DH.${CHR}.fa.gz
IND_LOC=./data/NA12878/
IND_NAME=${IND_LOC}NA12878
IND_S=${IND_LOC}/${CHR}/NA12878.chr${CHR}.s
IND_S_1X=${IND_S}.1x
REGION=chr${CHR}:16600000-16800000

# Filter out the region of interest and format to BAM
samtools view -T ${REF_GENOME} -bo ${IND_S}.bam ${IND_NAME}.final.cram ${REGION}
samtools index ${IND_S}.bam

# Get the genotype likelihood based on the panel for the validation file
PANEL_NAME=./data/panel/1000GP.chr21.noNA12878.s
VCF=${PANEL_NAME}.sites.vcf.gz
TSV=${PANEL_NAME}.sites.tsv.gz

bcftools mpileup -f ${REF_GENOME} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${IND_S}.bam -Ou |
bcftools call -Aim -C alleles -T ${TSV} -Ob -o ${IND_S}.bcf
bcftools index -f ${IND_S}.bcf

# Downsampling the individual data to 1X
samtools view -T ${REF_GENOME} -s 1.03045 -bo ${IND_S_1X}.bam ${IND_NAME}.final.cram ${REGION}
samtools index ${IND_S_1X}.bam