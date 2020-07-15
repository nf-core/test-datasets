# Use 1000 genomes to use as as shrinked set, in chunks and in a combined form. 
# A data reseource, meant to be used for development and testing.

# Download 22+X chromosomes
source_fold="data/data_source"
mkdir -p ${source_fold}
for chr in {1..22}; do
  echo ${chr}
  wget -P ${source_fold} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz*
done
wget -P ${source_fold} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz*

# Load modules (how software is installed on HPC)
module load tools
module load bcftools/1.9
module load tabix/1.2.1

# Shrink files to only contain 4500 lines each (to not take too much space)
chunk_fold="data/data_shrink_chunk_4500"
mkdir -p ${chunk_fold}
for chr in {1..22}; do
  fileToRead="data_source/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" 
  bcftools view -Ov ${fileToRead} | head -n4500 | bgzip -c > ${chunk_fold}/chr${chr}.vcf.bgz
  tabix -p vcf ${chunk_fold}/chr${chr}.vcf.bgz
done
fileToRead="data_source/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
bcftools view -Ov ${fileToRead} | head -n4500 | bgzip -c > ${chunk_fold}/chrX.vcf.bgz
tabix -p vcf ${chunk_fold}/chrX.vcf.bgz

# Use bcftools to combine the shrinked data
combine_fold="data/data_shrink_combined_4500"
mkdir -p ${combine_fold}
bcftools concat -Ov ${chunk_fold}/*.vcf.bgz | bgzip -c > ${combine_fold}/chr1_to_22_and_X.vcf.bgz
tabix -p vcf ${combine_fold}/chr1_to_22_and_X.vcf.bgz


