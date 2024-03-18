#!/bin/bash

#SBATCH --job-name=dwld
#SBATCH --output=dwld.out
#SBATCH --error=dwld.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

# Selected individuals
FL_IND=$1

# Location in the listing file
FL_LST=$2

# FTP server
FTP="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323"

# Download the selected individuals
while read -r IND; do
    echo $IND
    mkdir -p data/individuals/$IND
    grep -w $IND $FL_LST | xargs -I {} wget -c ${FTP}/{} -O data/individuals/$IND/$IND.{cram,cram.crai}
done < $FL_IND
