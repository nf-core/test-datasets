#!/bin/sh

#SBATCH --job-name=get_panel
#SBATCH --output=get_panel.out
#SBATCH --error=get_panel.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

conda init bash
conda activate env_tools

PANEL_DIR=$1
PANEL_NAME=$2
REF_FASTA=$3
REGION_LST=$4
PREFIX=$5

# Extract only necessary region from fasta
echo 'Extract region from fasta'
# Keep only chr
cat ${REGION_LST} | cut -d':' -f1 | sort | uniq > ${REGION_LST}.chr
samtools faidx ${REF_FASTA}.fa.bgz \
    --region-file ${REGION_LST}.chr --output ${REF_FASTA}.s.fa
bgzip -f ${REF_FASTA}.s.fa
samtools faidx ${REF_FASTA}.s.fa.gz

# Filter the region of interest of the panel file
echo 'Filter region of panel'

for REGION in $( cat $REGION_LST );
do
    if [ "$PREFIX" == "nochr" ]; then
        REGION=$(echo $REGION | sed 's/chr//')
    fi
    echo $REGION
    CHR=$(echo $REGION | cut -d':' -f1)
    PANEL_FILE=${PANEL_DIR}/${CHR}/${PANEL_NAME}.${CHR}
    REL_IND=$(cat ./analysis/all_rel_individuals.txt | tr '\n' ',' | sed 's/,$//')

    # Filter the panel file
    # Select the region of interest, annotate the variants, normalise the panel and filter out related individual to selected individuals
    bcftools norm -m +any ${PANEL_FILE}.vcf.gz \
        --regions $REGION --threads 4 -Ov | \
    bcftools view \
        -m 2 -M 2 -v snps -s ^${REL_IND} --force-samples --threads 4 | \
    bcftools annotate --threads 4 -Ov \
        --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' | \
    vcffixup - | bgzip -c > ${PANEL_FILE}.s.norel.vcf.gz

    # Index the panel file
    bcftools index -f ${PANEL_FILE}.s.norel.vcf.gz --threads 4

    # Get phased haplotypes
    #echo 'Get phased haplotypes'
    #SHAPEIT5_phase_common -I ${PANEL_FILE}.norel.bcf -T 4 -O ${PANEL_FILE}.phased.vcf.gz -R ${REGION}
    #bcftools index -f ${PANEL_FILE}.phased.vcf.gz

    echo "Chunk: ${REGION}"
    GLIMPSE_chunk \
        --input ${PANEL_FILE}.s.norel.vcf.gz --region ${REGION} \
        --window-size 10000 --window-count 400 --buffer-size 5000 --buffer-count 30 \
        --output ${PANEL_FILE}_chunks.txt
        #--sequential --window-mb 0.01 --window-cm 0.01 --window-count 200 --buffer-mb 0.0005 --buffer-cm 0.0005 --buffer-count 5 \
        

    # Select only the SNPS and drop Genotypes
    echo 'Select only SNPs and drop Genotypes'
    bcftools view -G -m 2 -M 2 -v snps ${PANEL_FILE}.s.norel.vcf.gz -Oz -o ${PANEL_FILE}.sites.vcf.gz
    bcftools index -f ${PANEL_FILE}.sites.vcf.gz

    # Convert to hap legend format
    echo 'Convert to hap legend format'
    bcftools convert --haplegendsample ${PANEL_NAME}.s.norel ${PANEL_FILE}.s.norel.vcf.gz -Oz -o ${PANEL_FILE}
done
