#!/bin/sh

#SBATCH --job-name=get_ind_imputed
#SBATCH --output=get_ind_imputed.out
#SBATCH --error=get_ind_imputed.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

#conda init bash
#conda activate env_tools

echo 'Get imputed individual data'

PANEL_NAME=$1
REF_FASTA=$2
REGION_LST=$3

# Selected individuals
INDS=($(cat ./analysis/selected_individuals.txt))

# Set region file to bcftools format
REGION_CSV=${REGION_LST}.csv
awk -v OFS='\t' -F':' '{split($2, coords, "-"); print $1, coords[1], coords[2]}' ${REGION_LST} > ${REGION_CSV}

VCF=${PANEL_NAME}.sites.vcf.gz
TSV=${PANEL_NAME}.tsv.gz

# Chunk panel file
echo 'Chunk panel file'
while IFS="\t" read REGION; do
    echo "Chunk: ${REGION}"
    GLIMPSE2_split_reference \
        --reference ${PANEL_NAME}.bcf --input-region ${REGION} --output-region ${REGION} \
        --output ${PANEL_NAME}_split
    REGIONN=$(echo ${REGION} | sed 's/:/_/' | sed 's/-/_/')
    while IFS="," read IND; do
        IND_S="./data/individuals/${IND}/${IND}.s"
        BAM="${IND_S}.1x.bam"
        # Impute with glimpse
        GLIMPSE2_phase --bam-file ${BAM} --ind-name ${IND} --reference ${PANEL_NAME}_split_${REGIONN}.bin --output ${IND_S}_${REGIONN}_imputed.bcf
    done < ./analysis/selected_individuals.txt

done < ${REGION_LST}

while IFS="," read IND; do
        IND_S="./data/individuals/${IND}/${IND}.s"
        # Impute with glimpse
        ls -1v ${IND_S}_*_imputed.bcf >  ${IND_S}_list.txt
        GLIMPSE2_ligate --input ${IND_S}_list.txt --output ${IND_S}_imputed.bcf
done < ./analysis/selected_individuals.txt

printf "chr21:16650000-16700000 ${PANEL_NAME}.bcf ${IND_S}.bcf ${IND_S}_imputed.bcf" > input.txt

GLIMPSE2_concordance \
        --bins 0 0.01 0.05 0.1 0.2 0.5 \
        --min-val-gl 0.9 \
        --min-val-dp 5 \
        --input input.txt \
        --thread 2 \
        --output NA12878