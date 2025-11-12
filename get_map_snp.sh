#!/bin/sh

#SBATCH --job-name=get_panel
#SBATCH --output=get_panel.out
#SBATCH --error=get_panel.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

conda init bash
conda activate env_tools

FOLDER=$1
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
        awk -F'\t' '$5 == "SNP" && $2 == '"$CHR_NUM"' { print "chr"$2":"$3}' \
        >> ${SNP_FILE}.s.map
    # Unzip the map file and keep only the chromosome file
    echo -e "pos\tchr\tcM" > ${FOLDER}/${REF_GEN}_${CHR}.glimpse.map
    unzip -p ${FOLDER}${REF_GEN}.map.zip chr_in_chrom_field/plink.chr${CHR}.${REF_GEN}.map | \
        awk -v OFS='\t' -F' ' '{
            if ($1 !~ /^chr/) $1 = "chr" $1
            print $4, $1, $3
        }' \
        >>  ${FOLDER}/${REF_GEN}_${CHR}.glimpse.map
    # Unzip the map file and keep all 4 PLINK columns + add chr prefix
    unzip -p ${FOLDER}${REF_GEN}.map.zip chr_in_chrom_field/plink.chr${CHR}.${REF_GEN}.map | \
        awk -v OFS=' ' -F' ' '{ 
            if ($1 !~ /^chr/) $1 = "chr" $1
            print $1, $2, $3, $4 
        }' \
        >  ${FOLDER}/${REF_GEN}_${CHR}.plink.map
done < $REGION_LST
