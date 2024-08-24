Reference panel
- Region: 21:16570000-16610000 & 22:16570000-16610000
- Reference genome: CanFam3.1 
- Number of samples: 658
- Source: Animal-ImputeDB

Reference genome
- Region: chr21:16570000-16610000 & chr22:16570000-16610000
- Reference genome: Broad CanFam3.1/canFam3
- Source: UCSC


```bash
# Download the reference panel and reference genome
mkdir -p data/dog/reference_genome data/dog/panel/21 data/dog/panel/22
wget -c -O- http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz | gunzip | bgzip > dog_data/reference_genome/canFam3.fa.bgz
wget -c -O- http://gong_lab.hzau.edu.cn/static/imputeDB/download/species/dog/panel/chr21_dog_impute.vcf.gz | gunzip | bgzip > dog_data/panel/21/658_dog.21.vcf.gz
wget -c -O- http://gong_lab.hzau.edu.cn/static/imputeDB/download/species/dog/panel/chr22_dog_impute.vcf.gz | gunzip | bgzip > dog_data/panel/22/658_dog.22.vcf.gz

bcftools index dog_data/panel/21/658_dog.21.vcf.gz
bcftools index dog_data/panel/22/658_dog.22.vcf.gz


# Subset the files
. get_panel_s.sh \
    dog_data/panel \
    658_dog \
    dog_data/reference_genome/canFam3 \
    region.lst \
    nochr
```

Individual data
- Region: chr21:16570000-16610000 & chr22:16570000-16610000
- Source: iDog (https://ngdc.cncb.ac.cn/idog/)
- License: Free for academic use

```
mkdir -p dog_data/csv dog_data/individuals/1735
cd dog_data/individuals/1735
wget ftp://download.big.ac.cn/idog/dogsd/bam/1735.rmdup.realign.recalibration.bam
samtools index 1735.rmdup.realign.recalibration.bam
samtools view -b 1735.rmdup.realign.recalibration.bam "chr21:16570000-16610000" "chr22:16570000-16610000" > 1735_subset_chr21_chr22.bam
samtools index 1735_subset_chr21_chr22.bam
```