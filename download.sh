#!/bin/bash

#SBATCH --job-name=dwld
#SBATCH --output=dwld.out
#SBATCH --error=dwld.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

mkdir -p data/NA19401 data/NA20359

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239749/NA19401.final.cram -O data/NA19401/NA19401.final.cram
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239996/NA20359.final.cram -O data/NA20359/NA20359.final.cram
