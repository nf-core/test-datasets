#!/bin/bash

# Directory containing bam files
BAM_DIR="out/GSE182201/star_salmon/genome"

# Directory containing fastq files
FASTQ_DIR="out/GSE182201/fastq"

# Directory for storing output
OUT_DIR="test_data"
TMP_DIR="test_data/tmp"

mkdir -p $TMP_DIR

# Filter bam files for chromosome 20 and extract sequence names
for bam_file in $BAM_DIR/*.bam; do
    base_name=$(basename "$bam_file" _genome.bam)
    echo "Generating sequence lists for $base_name"
    if [ ! -e "$TMP_DIR/${base_name}_20.txt" ]; then
        samtools view -h "$bam_file" "20" | samtools view -b - > "$TMP_DIR/${base_name}_20.bam"
        samtools view "$TMP_DIR/${base_name}_20.bam" | cut -f 1 > "$TMP_DIR/${base_name}_20.txt"
    fi
done

# Loop through fastq files and filter based on sequence names
for fastq_file in $FASTQ_DIR/*_1.fastq.gz; do
    base_name=$(basename "$fastq_file" _1.fastq.gz)
    bam_name="${base_name%_SRR*}"
    echo "Subsetting $base_name"

    r1=${FASTQ_DIR}/${base_name}_1.fastq.gz
    r2=${FASTQ_DIR}/${base_name}_2.fastq.gz
    seq_selection=$TMP_DIR/${bam_name}_20.txt
    if [ ! -e "$seq_selection" ]; then
        echo "Can't find seq selection $seq_selection"
	exit 1
    fi

    if [ -e "$r1" ]; then
	echo "Subsetting $r1 by $seq_selection"
	zcat "$r1" | seqtk subseq - "$seq_selection" | head -n 2000000 | gzip > "$OUT_DIR/${base_name}_chr20_1.fastq.gz"
    else
	echo "$r1 does not exist"
    fi
    if [ -e "$r2" ]; then
	echo "Subsetting $r2 by $seq_selection"
	zcat "$r2" | seqtk subseq - "$seq_selection" | head -n 2000000 | gzip > "$OUT_DIR/${base_name}_chr20_2.fastq.gz"
    fi
done

# Remove temporary directory
#rm -r "$TMP_DIR"
