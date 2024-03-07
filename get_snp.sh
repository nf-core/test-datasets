#!/bin/sh

# Get the file of all SNP and CN
CHR=22
REGION=chr${CHR}:16600000-16800000
SNP_FILE=./data/affi/snp6

# Select only chr 21 in the region of interest, filter out only SNP and concatenate it as chrX:X
cat ${SNP_FILE}.txt | awk -F'\t' '$5 == "SNP" && $2 == '"$CHR"' && $3 >=16600000 && $3 <= 16800000 { print "chr"$2":"$3}' > ${SNP_FILE}.chr${CHR}.s.map