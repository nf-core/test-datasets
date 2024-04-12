#!/bin/sh

#SBATCH --job-name=get_ind
#SBATCH --output=get_ind.out
#SBATCH --error=get_ind.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

#conda init bash
#conda activate env_tools

PANEL_NAME=$1
REF_FASTA=$2
SNP_FILE=$3
REGION_LST=$4

# Selected individuals
INDS=($(cat ./analysis/selected_individuals.txt))

# Set region file to bcftools format
REGION_CSV=${REGION_LST}.csv
awk -v OFS='\t' -F':' '{split($2, coords, "-"); print $1, coords[1], coords[2]}' ${REGION_LST} > ${REGION_CSV}

while IFS="," read IND; do
    echo "Individuals: ${IND}"
    IND_LOC=./data/individuals/${IND}/
    IND_CRAM=${IND_LOC}${IND}.final.cram
    IND_S=${IND_LOC}${IND}.s

    # Filter out the region of interest and format to BAM
    echo 'Filter out the region of interest and format to BAM'
    samtools view \
        --regions-file ${REGION_CSV} \
        -bo ${IND_S}.bam ${IND_CRAM}
    samtools index ${IND_S}.bam

    # Get the genotype likelihood based on the panel for the validation file
    echo 'Get the genotype likelihood based on the panel for the validation file'
    VCF=${PANEL_NAME}.sites.vcf.gz
    TSV=${PANEL_NAME}.tsv.gz

    bcftools mpileup -f ${REF_FASTA} \
        -I -E -a 'FORMAT/DP' -T ${VCF} \
        ${IND_S}.bam -Ou | \
        bcftools call -Aim -C alleles -T ${TSV} | \
        bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' -Ob -o ${IND_S}.bcf
    bcftools index -f ${IND_S}.bcf

    # Downsampling the individual data to 1X
    echo 'Downsampling the individual data to 1X'
    MEAN_DEPTH=$(samtools coverage ${IND_S}.bam | \
        awk -F'\t' '(NR==2){ print $7}')
    FRAC_DEPTH=$(echo "scale=5; 1/$MEAN_DEPTH" | bc)

    # Downsample to 1X and convert to CRAM or BAM
    samtools view -T ${REF_FASTA} \
        --regions-file ${REGION_CSV} -s 1${FRAC_DEPTH} \
        -bo ${IND_S}.1x.bam ${IND_S}.bam
    samtools index ${IND_S}.1x.bam

    # Compute genotype likelihood based on the panel
    echo 'Compute genotype likelihood based on the panel'
    bcftools mpileup -f ${REF_FASTA} \
        -I -E -a 'FORMAT/DP' -T ${VCF} \
        --regions-file ${REGION_CSV} ${IND_S}.1x.bam -Ou |
        bcftools call -Aim -C alleles -T ${TSV} | \
        bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' -Ob -o ${IND_S}.1x.bcf
    bcftools index -f ${IND_S}.1x.bcf

    # Get individual SNP
    echo 'Get individual SNP'
    plink --bcf ${IND_S}.bcf \
        --allow-extra-chr \
        --geno 0 \
        --make-bed --out ${IND_S}
    # Replace all variants id by .
    awk '{$2 = "."; print}' ${IND_S}.bim > ${IND_S}_NID.bim
    plink --bfile ${IND_S} \
        --bim ${IND_S}_NID.bim \
        --allow-extra-chr \
        --set-missing-var-ids chr@:# \
        --extract ${SNP_FILE} \
        --recode vcf-iid bgz\
        --output-chr chrM \
        --out ${IND_S}.snp

    bcftools index -f ${IND_S}.snp.vcf.gz

    #bcftools view -e 'GT="./."||GT="."' ${IND_S}.snp.vcf.gz -Oz -o ${IND_S}.snp.filtered.vcf.gz
    #bcftools index -f ${IND_S}.snp.filtered.vcf.gz
done < ./analysis/selected_individuals.txt