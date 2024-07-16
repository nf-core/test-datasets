#!/bin/bash

#SBATCH --job-name=dwld
#SBATCH --output=dwld.out
#SBATCH --error=dwld.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

FL_IND=$1 # Selected individuals
DIR_IND=$2 # Directory to save the individuals
FL_LST=$3 # Location of the listing file
FTP=$4 # FTP server

# Download the selected individuals
while read -r IND; do
    echo $IND
    mkdir -p $DIR_IND/$IND
    grep -w $IND $FL_LST | xargs -I {} wget -c ${FTP}/{} -O $DIR_IND/$IND/$IND.{cram,cram.crai}
done < $FL_IND
