Reference panel
- Region: chr1:24300000-24340000. 
- Reference genome: CanFam3.1 
- Number of samples: 658
- Source: Animal-ImputeDB

Reference genome
- Region: chr1:23719985-24920014
- Reference genome: Broad CanFam3.1/canFam3
- Source: UCSC


```bash
mkdir -p data/dog/reference_genome data/dog/panel/21 data/dog/panel/22
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/285/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_genomic.fna.gz -O data/dog/reference_genome/canFam3.fa.gz
wget -c http://gong_lab.hzau.edu.cn/static/imputeDB/download/species/dog/panel/chr21_dog_impute.vcf.gz -O data/dog/panel/21/652_dog.chr21.vcf.gz
wget -c http://gong_lab.hzau.edu.cn/static/imputeDB/download/species/dog/panel/chr22_dog_impute.vcf.gz -O data/dog/panel/22/652_dog.chr22.vcf.gz

gunzip data/dog/reference_genome/canFam3.fa.gz | \
samtools faidx data/dog/reference_genome/canFam3.fa --region-file region.lst --output data/dog/reference_genome/canFam3.s.fa
samtools faidx data/dog/reference_genome/canFam3.s.fa
```
