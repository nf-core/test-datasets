#!/bin/sh
conda init bash
conda activate env_tools

CHR=22
REF_GENOME=./data/reference_genome/hs38DH.chr${CHR}.fa.gz
IND_LOC=./data/NA12878/
IND_NAME=${IND_LOC}NA12878
IND_S=${IND_LOC}${CHR}/NA12878.chr${CHR}.s
IND_S_1X=${IND_S}.1x
REGION=chr${CHR}:16600000-16800000

# Filter out the region of interest and format to BAM
echo 'Filter out the region of interest and format to BAM'
samtools view -T ${REF_GENOME} -bo ${IND_S}.bam ${IND_NAME}.final.cram ${REGION}
samtools index ${IND_S}.bam

# Get the genotype likelihood based on the panel for the validation file
echo 'Get the genotype likelihood based on the panel for the validation file'
PANEL_NAME=./data/panel/${CHR}/1000GP.chr${CHR}.noNA12878.s
VCF=${PANEL_NAME}.sites.vcf.gz
TSV=${PANEL_NAME}.sites.tsv.gz

bcftools mpileup -f ${REF_GENOME} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${IND_S}.bam -Ou |
bcftools call -Aim -C alleles -T ${TSV} -Ob -o ${IND_S}.bcf
bcftools index -f ${IND_S}.bcf

# Downsampling the individual data to 1X
echo 'Downsampling the individual data to 1X'
MEAN_DEPTH=$(samtools coverage ${IND_S}.bam -r ${REGION} | \
    awk -F'\t' '(NR==2){ print $7}')
FRAC_DEPTH=$(echo "scale=5; 1/$MEAN_DEPTH" | bc)
samtools view -T ${REF_GENOME} -s 1${FRAC_DEPTH} -bo ${IND_S_1X}.bam ${IND_NAME}.final.cram ${REGION}
samtools index ${IND_S_1X}.bam

# Compute genotype likelihood based on the panel
echo 'Compute genotype likelihood based on the panel'
bcftools mpileup -f ${REF_GENOME} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${IND_S_1X}.bam -Ou |
bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${IND_S_1X}.vcf.gz
bcftools index -f ${IND_S_1X}.vcf.gz

# Get individual SNP
echo 'Get individual SNP'
MAP_SNP=./data/affi/snp6.chr${CHR}.s.map
plink2 --chr-set 22 --bcf ${IND_S}.bcf \
    --set-missing-var-ids @:# \
    --max-alleles 2 \
    --output-chr '26' \
    --snps-only --allow-extra-chr \
    --extract ${MAP_SNP} \
    --chr 1-22,X --recode vcf bgz\
    --out ${IND_S}.snp

bcftools index -f ${IND_S}.snp.vcf.gz

bcftools view -e 'GT="./."||GT="."' ${IND_S}.snp.vcf.gz -Oz -o ${IND_S}.snp.filtered.vcf.gz
bcftools index -f ${IND_S}.snp.filtered.vcf.gz