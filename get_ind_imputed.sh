#!/bin/sh

#SBATCH --job-name=get_ind_imputed
#SBATCH --output=get_ind_imputed.out
#SBATCH --error=get_ind_imputed.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

#conda init bash
#conda activate env_tools

echo 'Get imputed individual data'

PANEL_DIR=$1
PANEL_NAME=$2
DIR_IND=$3
REF_FASTA=$4
REGION_LST=$5
IND_LST=$6

# Chunk panel file
echo 'Chunk panel file'

while IFS="," read IND; do
    echo "Individual: ${IND}"
    IND_S="${DIR_IND}/${IND}/${IND}.s"
    IND_TMP="${DIR_IND}/${IND}/tmp"
    BAM="${IND_S}.1x.bam"
    while IFS="\t" read REGION; do
        echo "Region: ${REGION}"
        CHR=$(echo $REGION | cut -d':' -f1)
        echo "Chromosome: ${CHR}"
        mkdir -p ${IND_TMP}/${CHR}
        PANEL_FILE=${PANEL_DIR}/${CHR}/${PANEL_NAME}.${CHR}.s.norel.vcf.gz
        awk -v OFS='\t' '{print $1,$3,$4}' ${PANEL_DIR}/${CHR}/${PANEL_NAME}.${CHR}_chunks.txt > ${IND_TMP}/${IND}_${CHR}_chunks.txt
        while IFS=$'\t' read REG INPUT OUTPUT; do
            echo "Chunk: reg:${REG} in:${INPUT} out:${OUTPUT}"
            # Impute with glimpse
            GLIMPSE2_phase --bam-file ${BAM} --ind-name ${IND} \
                --reference ${PANEL_FILE} --input-region ${INPUT} --output-region ${OUTPUT} \
                --keep-monomorphic-ref-sites \
                --output ${IND_TMP}/${CHR}/${IND}_${REG}_imputed.bcf
        done < ${IND_TMP}/${IND}_${CHR}_chunks.txt
         # Impute with glimpse
        ls -1v ${IND_TMP}/${CHR}/${IND}_*_imputed.bcf >  ${IND_TMP}/${IND}_${CHR}_list.txt
        GLIMPSE2_ligate --input ${IND_TMP}/${IND}_${CHR}_list.txt --output ${IND_TMP}/${IND}_${CHR}_imputed.bcf
    done < ${REGION_LST}
    # Merge all regions
    ls -1v ${IND_TMP}/${IND}_chr*_imputed.bcf >  ${IND_TMP}/${IND}_list.txt
    GLIMPSE2_ligate --input ${IND_TMP}/${IND}_list.txt --output ${IND_S}_imputed.bcf
done < ${IND_LST}
