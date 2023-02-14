#!/bin/sh

# Get the file of all SNP and CN
wget https://api.gdc.cancer.gov/data/9bd7cbce-80f9-449e-8007-ddc9b1e89dfb snp6.txt.gz
gunzip snp6.txt.gz

# Filter out only SNP and concatenate it as chrX:X
cat snp6.txt | awk -F'\t' '$5 == "SNP" { print "chr"$2":"$3}' > snp6.map