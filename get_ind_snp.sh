#!/bin/sh

#SBATCH --job-name=get_ind
#SBATCH --output=get_ind.out
#SBATCH --error=get_ind.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

#conda init bash
#conda activate env_tools

DIR_IND=$1
DIR_PANEL=$2
PANEL_NAME=$3
REF_FASTA=$4
SNP_FILE=$5
IND_SEL=$6
REGION_LST=$7

# Create full tsv file from legend file
> ${DIR_PANEL}/${PANEL_NAME}.tsv
for REGION in $( cat $REGION_LST );
do
    CHR=$(echo $REGION | cut -d':' -f1)
    zcat ${DIR_PANEL}/$CHR/${PANEL_NAME}.${CHR}.legend.gz | \
        awk -v OFS='\t' 'NR>1 { split($1, a, "[:-_]"); print a[1], $2, $3 "," $4 }' >> ${DIR_PANEL}/${PANEL_NAME}.tsv
done

bgzip ${DIR_PANEL}/${PANEL_NAME}.tsv
tabix -s1 -b2 -e2 ${DIR_PANEL}/${PANEL_NAME}.tsv.gz

TSV=${DIR_PANEL}/${PANEL_NAME}.tsv.gz

while IFS="," read IND; do
    # Get the genotype likelihood based on the panel for the validation file
    echo 'Get the genotype likelihood based on the panel for the validation file'

    IND_DIR=${DIR_IND}/${IND}
    IND_FILE=${IND_DIR}/${IND}.s

    bcftools mpileup -f ${REF_FASTA} \
        -I -E -a 'FORMAT/DP' -T ${TSV} \
        ${IND_S}.bam | \
        bcftools call -Aim -C alleles -T ${TSV} | \
        bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' -Ob -o ${IND_S}.bcf

    bcftools index -f ${IND_S}.bcf

    # Get individual SNP
    echo 'Get individual SNP'
    mkdir -p ${IND_DIR}/tmp
    plink --bcf ${IND_S}.bcf \
        --allow-extra-chr \
        --geno 0 \
        --make-bed --out ${IND_DIR}/tmp/${IND}

    # Replace all variants id by . to filter them with snp array
    awk '{$2 = "."; print}' ${IND_DIR}/tmp/${IND}.bim > ${IND_DIR}/tmp/${IND}_NID.bim
    plink --bfile ${IND_DIR}/tmp/${IND} \
        --bim ${IND_DIR}/tmp/${IND}_NID.bim \
        --allow-extra-chr \
        --set-missing-var-ids chr@:# \
        --extract ${SNP_FILE} \
        --recode vcf-iid bgz\
        --output-chr chrM \
        --out ${IND_S}.snp

    bcftools index -f ${IND_S}.snp.vcf.gz

    #bcftools view -e 'GT="./."||GT="."' ${IND_S}.snp.vcf.gz -Oz -o ${IND_S}.snp.filtered.vcf.gz
    #bcftools index -f ${IND_S}.snp.filtered.vcf.gz
done < $IND_SEL