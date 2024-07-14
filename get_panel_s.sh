#!/bin/sh

#SBATCH --job-name=get_panel
#SBATCH --output=get_panel.out
#SBATCH --error=get_panel.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

conda init bash
conda activate env_tools

PANEL_NAME=$1
REF_FASTA=$2
REGION_LST=$3
PANEL_NOREL=${PANEL_NAME}.s.norel

# Extract only necessary region from fasta
echo 'Extract region from fasta'
samtools faidx ${REF_FASTA}.fa \
    --region-file ${REGION_LST} --output ${REF_FASTA}.s.fa
samtools faidx ${REF_FASTA}.s.fa

# Filter the region of interest of the panel file
## Convert REGION_LST to tab delimited format
REGION_CSV=${REGION_LST}.csv
awk -v OFS='\t' -F':' '{split($2, coords, "-"); print $1, coords[1], coords[2]}' ${REGION_LST} > ${REGION_CSV}

echo 'Filter region of panel'
bcftools view ${PANEL_NAME}.vcf.gz \
    --regions-file ${REGION_CSV} | \
    bcftools annotate \
        --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' \
        -Oz -o ${PANEL_NAME}.s.vcf.gz
bcftools index -f ${PANEL_NAME}.s.vcf.gz --threads 4

# Normalise the panel and filter out related individual to selected individuals
# Read the all the related individuals
REL_IND=$(cat ./analysis/all_rel_individuals.txt | tr '\n' ',' | sed 's/,$//')

echo 'Normalise panel and take out related individual'
bcftools norm -m -any ${PANEL_NAME}.s.vcf.gz -Ou --threads 4 |
bcftools view -m 2 -M 2 -v snps -s ^${REL_IND} --threads 4 -Ob -o ${PANEL_NOREL}.bcf
bcftools index -f ${PANEL_NOREL}.bcf --threads 4

# Select only the SNPS and drop Genotypes
echo 'Select only SNPs and drop Genotypes'
bcftools view -G -m 2 -M 2 -v snps ${PANEL_NOREL}.bcf -Oz -o ${PANEL_NOREL}.sites.vcf.gz
bcftools index -f ${PANEL_NOREL}.sites.vcf.gz

# Convert to TSV
echo 'Convert to TSV'
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${PANEL_NOREL}.sites.vcf.gz | bgzip -c > ${PANEL_NOREL}.tsv.gz
tabix -s1 -b2 -e2 ${PANEL_NOREL}.tsv.gz

# Get phased haplotypes
#echo 'Get phased haplotypes'
#SHAPEIT5_phase_common -I ${PANEL_NOREL}.bcf -T 4 -O ${PANEL_NOREL}.phased.vcf.gz -R ${REGION}
#bcftools index -f ${PANEL_NOREL}.phased.vcf.gz
