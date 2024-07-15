#!/bin/sh

#SBATCH --job-name=get_panel
#SBATCH --output=get_panel.out
#SBATCH --error=get_panel.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

conda init bash
conda activate env_tools

REF_MAP=$1
REF_GEN=$2
SNP_FILE=$3
REGION_LST=$4

# Get SNP file
echo 'Extracting region from SNP array file'

# For each line in the region file, extract the SNPs and map file

# Initialize the SNP and map file
> ${SNP_FILE}.s.map

while IFS=':' read -r CHR REGION; do
    # Use col1 and col2 as arguments
    START=$(echo $REGION | awk -F'-' '{print $1}')
    END=$(echo $REGION | awk -F'-' '{print $2}')
    CHR_NUM=$(echo $CHR | sed 's/chr//g')
    echo "$CHR: $START - $END"
    # Extract the SNPs
    zcat ${SNP_FILE}.txt.gz | \
        awk -F'\t' '$5 == "SNP" && $2 == '"$CHR_NUM"' { print $2":"$3}' \
        >> ${SNP_FILE}.s.map
    # Unzip the map file and keep only the chromosome file
    unzip -p ${REF_MAP}${REF_GEN}.map.zip plink.${CHR}.${REF_GEN}.map | \
        awk -v OFS='\t' -F' ' '{ print $1, $3, $4 }' \
        >  ${REF_MAP}/${REF_GEN}_${CHR_NUM}.map
done < $REGION_LST
