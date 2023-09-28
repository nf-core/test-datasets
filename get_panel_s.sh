#!/bin/sh
conda init bash
conda activate env_tools

CHR=22
PANEL_LOC=./data/panel/${CHR}/
PANEL_ORIGIN=${PANEL_LOC}panel_2020-08-05_chr${CHR}.phased.vcf.gz
PANEL_NAME=${PANEL_LOC}1000GP.chr${CHR}
PANEL_S=${PANEL_NAME}.noNA12878.s

REGION=chr${CHR}:16600000-16800000

# Filter the region of interest of the panel file
echo 'Filter region of panel'
bcftools view ${PANEL_ORIGIN} -r ${REGION} -Oz -o ${PANEL_S}.vcf.gz
bcftools index -f ${PANEL_S}.vcf.gz --threads 4

# Normalise the panel and filter out related individual to NA12878
echo 'Normalise panel and take out related individual'
bcftools norm -m -any ${PANEL_S}.vcf.gz -Ou --threads 4 |
bcftools view -m 2 -M 2 -v snps -s ^NA12878,NA12891,NA12892 --threads 4 -Ob -o ${PANEL_S}.bcf
bcftools index -f ${PANEL_S}.bcf --threads 4

# Select only the SNPS
echo 'Select only SNPs'
bcftools view -G -m 2 -M 2 -v snps ${PANEL_S}.bcf -Oz -o ${PANEL_S}.sites.vcf.gz
bcftools index -f ${PANEL_S}.sites.vcf.gz

# Convert to TSV
echo 'Convert to TSV'
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${PANEL_S}.sites.vcf.gz | bgzip -c > ${PANEL_S}.sites.tsv.gz
tabix -s1 -b2 -e2 ${PANEL_S}.sites.tsv.gz

# Get phased haplotypes
#echo 'Get phased haplotypes'
#SHAPEIT5_phase_common -I ${PANEL_S}.bcf -T 4 -O ${PANEL_S}.phased.vcf.gz -R ${REGION}
#bcftools index -f ${PANEL_S}.phased.vcf.gz
