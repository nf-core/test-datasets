#!/bin/sh

#SBATCH --job-name=get_ind
#SBATCH --output=get_ind.out
#SBATCH --error=get_ind.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

#conda init bash
#conda activate env_tools

DIR_IND=$1
REF_FASTA=$2
IND_SEL=$3
REGION_LST=$4

# Set region file to bcftools format
REGION_CSV=${REGION_LST}.csv
awk -v OFS='\t' -F':' '{split($2, coords, "-"); print $1, coords[1], coords[2]}' ${REGION_LST} > ${REGION_CSV}

while IFS="," read IND; do
    echo "Individuals: ${IND}"
    IND_LOC=${DIR_IND}/${IND}/
    IND_CRAM=${IND_LOC}${IND}.final.cram
    IND_S=${IND_LOC}${IND}.s

    # Filter out the region of interest and format to BAM
    echo 'Filter out the region of interest and format to BAM'
    samtools view \
        --regions-file ${REGION_CSV} \
        -bo ${IND_S}.bam ${IND_CRAM}
    samtools index ${IND_S}.bam

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

done < $IND_SEL
