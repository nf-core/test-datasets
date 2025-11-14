# This file contains information on how the nallo testdata files are created.
# If new regions need to be added, try to keep them in the beginning of the chromosomes.

mkdir tmp

# 1. From somalier sites (https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz), extract sites to make sex inference work.
zcat data/sites.hg38.vcf.gz | grep -v "#" | awk -v OFS='\t' '{print $1,$2-1,$2+1}' > tmp/somalier_sites.bed
cat tmp/somalier_sites.bed | grep "chrX" | head -35 > tmp/chrX.bed
cat tmp/somalier_sites.bed | grep "chr16" | head -35 > tmp/chr16.bed
cat tmp/chr16.bed tmp/chrX.bed > tmp/somalier_sites_cut.bed

# 2. Make a new somalier_sites containing only the variants we cut out
bcftools view -R tmp/somalier_sites_cut.bed data/sites.hg38.vcf.gz -Oz > tmp/somalier_sites.vcf.gz

# 3. Add paraphase (paraphase.bed), and general (test_data.bed) regions

cat <<EOF > tmp/test_data.bed
chr16	160000	173710
chr16	176680	177522
chr20	2652400	2653100
EOF

# HBA1, HBA2
cat <<EOF > tmp/paraphase.bed
chr16	172876	173710
chr16	176680	177522
chr17	3053138	3073138
EOF

cat tmp/paraphase.bed tmp/test_data.bed tmp/chr16.bed tmp/chrX.bed > tmp/reference_regions.bed
cat tmp/paraphase.bed tmp/test_data.bed tmp/chr16.bed tmp/chrX.bed > tmp/test_somalier_small.bed
# 4. Use the BED file to cut out regions in the BAM files

# Merge the regions to create masking for reference
# Here we will chop the chromosomes of at the end of the regions
bedtools sort -i tmp/reference_regions.bed > tmp/reference_regions.sorted.bed
bedtools merge -i tmp/reference_regions.sorted.bed -d 999999999 |\
    awk -v OFS='\t' '{print $1,0,$3}' > tmp/chromosomes_with_sequence.bed

# 5. Make a reference with the cut out chromosomes,
# 5.1 mask everything which is not in the bed files
# 5.2 fix chr names,
# 5.3 add extra_chr.fa:

cat <<EOF > tmp/extra_chr.fa
>chr1
N
>chr2
N
>chr3
N
>chr4
N
>chr5
N
>chr6
N
>chr7
N
>chr8
N
>chr9
N
>chr10
N
>chr11
N
>chr12
N
>chr13
N
>chr14
N
>chr15
N
>chr18
N
>chr19
N
>chr21
N
>chr22
N
>chrY
N
>chrMT
N
EOF

seqtk subseq data/GRCh38_no_alt_analysis_set.fasta tmp/chromosomes_with_sequence.bed | seqtk seq -c -M tmp/test_somalier_small.bed | sed '/^>/ s/:.*//' | cat - tmp/extra_chr.fa | pigz -p 36 > tmp/hg38.test.fa.gz

# Subsample to reduce size
samtools view -@ 36 -s 0.8 -M -L tmp/test_somalier_small.bed data/HG002_haplotagged.bam -h -O BAM -u -x HP,PS \
  | samtools ampliconclip -b clipping.bed --both-ends --tolerance 1 -u \
  | samtools calmd -b - data/GRCh38_no_alt_analysis_set.fasta > tmp/hg002_somalier_small_revio.bam
samtools view -@ 36 -s 0.8 -M -L tmp/test_somalier_small.bed data/HG003_haplotagged.bam -h -O BAM -u -x HP,PS \
  | samtools ampliconclip -b clipping.bed --both-ends --tolerance 1 -u \
  | samtools calmd -b - data/GRCh38_no_alt_analysis_set.fasta > tmp/hg003_somalier_small_revio.bam
samtools view -@ 36 -s 0.6 -M -L tmp/test_somalier_small.bed data/HG004_haplotagged.bam -h -O BAM -u -x HP,PS \
  | samtools ampliconclip -b clipping.bed --both-ends --tolerance 1 -u \
  | samtools calmd -b - data/GRCh38_no_alt_analysis_set.fasta > tmp/hg004_somalier_small_revio.bam
samtools view -@ 36 -s 0.5 -M -L tmp/test_somalier_small.bed data/hg002.haplotagged.bam -h -O BAM -u -x HP,PS \
  | samtools ampliconclip -b clipping.bed --both-ends --tolerance 1 -u \
  | samtools calmd -b - data/GRCh38_no_alt_analysis_set.fasta > tmp/hg002_somalier_small_ont.bam

# Make fastq
samtools fastq -T \* -@ 36 tmp/hg002_somalier_small_revio.bam | pigz -p 36 > tmp/hg002_somalier_small_revio.fastq.gz
samtools fastq -T \* -@ 36 tmp/hg002_somalier_small_ont.bam |chopper --maxlength 200000| pigz -p 36 > tmp/hg002_somalier_small_ont.fastq.gz

echo "reads Revio HG002, HG003, HG004:"
samtools view tmp/hg002_somalier_small_revio.bam|wc -l
samtools view tmp/hg003_somalier_small_revio.bam|wc -l
samtools view tmp/hg004_somalier_small_revio.bam|wc -l
echo "reads ONT HG002"
samtools view tmp/hg002_somalier_small_ont.bam|wc -l

# 6. Create a SVDB file from https://zenodo.org/records/11511513/files/CoLoRSdb.GRCh38.v1.0.0.pbsv.jasmine.vcf.gz
bcftools view -R tmp/test_somalier_small.bed  data/CoLoRSdb.GRCh38.v1.0.0.pbsv.jasmine.vcf.gz --output-type z --output tmp/colorsdb.test_data.vcf.gz

# 7. Copy files, and make one copy smaller
cp tmp/colorsdb.test_data.vcf.gz reference/colorsdb.test_data.vcf.gz
cp tmp/hg38.test.fa.gz reference/hg38.test.fa.gz
cp tmp/test_somalier_small.bed reference/test_data.bed

cp tmp/hg002_somalier_small_revio.bam testdata/HG002_PacBio_Revio.bam
cp tmp/hg003_somalier_small_revio.bam testdata/HG003_PacBio_Revio.bam
cp tmp/hg004_somalier_small_revio.bam testdata/HG004_PacBio_Revio.bam
samtools index testdata/HG002_PacBio_Revio.bam
samtools index testdata/HG003_PacBio_Revio.bam
samtools index testdata/HG004_PacBio_Revio.bam
cp tmp/hg002_somalier_small_revio.fastq.gz testdata/HG002_PacBio_Revio.fastq.gz
samtools view -@ 36 -s 0.01 -h -O BAM tmp/hg002_somalier_small_revio.bam > testdata/HG002_PacBio_Revio_copy.bam
samtools index testdata/HG002_PacBio_Revio_copy.bam

cp tmp/hg002_somalier_small_ont.bam testdata/HG002_ONT.bam
samtools index testdata/HG002_ONT.bam
samtools view -@ 36 -s 0.01 -h -O BAM tmp/hg002_somalier_small_ont.bam > testdata/HG002_ONT_copy.bam
samtools index testdata/HG002_ONT_copy.bam
cp tmp/hg002_somalier_small_ont.fastq.gz testdata/HG002_ONT.fastq.gz

# 8. Unzip reference file
gunzip -c reference/hg38.test.fa.gz > reference/hg38.test.fa

# 9. Gnomad test data set
bcftools view --regions-file tmp/test_data.bed data/gnomad.genomes.v4.1.sites.chr16.vcf.bgz --output-type b --output tmp/gnomad_chr16.test_data.bcf.gz
bcftools view --regions-file tmp/test_data.bed data/gnomad.genomes.v4.1.sites.chrX.vcf.bgz --output-type b --output tmp/gnomad_chrX.test_data.bcf.gz
bcftools concat tmp/gnomad_chr16.test_data.bcf.gz tmp/gnomad_chrX.test_data.bcf.gz --output-type b --output tmp/gnomad_chr16_chrX.test_data.bcf.gz
docker run -v $(pwd):/data docker.io/fellen31/echtvar:0.2.2 echtvar encode /data/tmp/gnomad_af.zip /data/assets/gnomad.json /data/tmp/gnomad_chr16_chrX.test_data.bcf.gz
cp tmp/gnomad_af.zip reference/gnomad_af.zip