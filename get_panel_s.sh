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
REF_MAP=$3
SNP_FILE=$4
REGION_LST=$5
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
    --regions-file ${REGION_CSV} \
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
#SHAPEIT5_phase_common -I ${PANEL_S}.bcf -T 4 -O ${PANEL_S}.phased.vcf.gz -R ${REGION}
#bcftools index -f ${PANEL_S}.phased.vcf.gz

# Get SNP file
echo 'Extracting region from SNP array file'

# For each line in the region file, extract the SNPs and map file

# Initialize the SNP and map file
> ${SNP_FILE}.s.map
> ${REF_MAP}.s.map

while IFS=':' read -r CHR REGION; do
    # Use col1 and col2 as arguments
    START=$(echo $REGION | awk -F'-' '{print $1}')
    END=$(echo $REGION | awk -F'-' '{print $2}')
    CHR_NUM=$(echo $CHR | sed 's/chr//g')
    echo "$CHR: $START - $END"
    # Extract the SNPs
    cat ${SNP_FILE}.txt | \
        awk -F'\t' '$5 == "SNP" && $2 == '"$CHR_NUM"' && $3 >='"$START"' && $3 <= '"$END"' { print "chr"$2":"$3}' \
        >> ${SNP_FILE}.s.map
    # Unzip the map file and keep only the chromosome file
    unzip -p ${REF_MAP}.map.zip plink.${CHR}.GRCh38.map | \
        awk -v OFS='\t' -F' ' '$4 >='"$START"' && $4 <= '"$END"' { print $1, $3, $4}' \
    >> ${REF_MAP}.s.map
done < $REGION_LST
unzip -p ${REF_MAP}.zip plink.${CHR}.GRCh38.map > test.map