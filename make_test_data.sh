# This file contains information on how the nallo testdata files are created.
# If new regions need to be added, try to keep them in the beginning of the chromosomes. 

# 1. From somalier sites (https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz), extract sites to make sex inference work. 
zless sites.hg38.vcf.gz| grep -v "#" | awk -v OFS='\t' '{print $1,$2-1,$2+1}' > somalier_sites.bed
less somalier_sites.bed | grep "chrX" | head -35 > chrX.bed
less somalier_sites.bed | grep "chr16" | head -35 > chr16.bed
cat chr16.bed chrX.bed > somalier_sites_cut.bed

# 2. Make a new somalier_sites containing only the variants we cut out
bcftools view -R somalier_sites_cut.bed sites.hg38.vcf.gz -Oz > somalier_sites.vcf.gz

# 3. Add paraphase (paraphase.bed), and general (test_data.bed) regions
# 
# test_data.bed
# chr16	160000	173710
# chr16	176680	177522
# chr20	2652400	2653100
# 
# paraphase.bed
# chr16	172876	173710
# chr16	176680	177522
cat paraphase.bed test_data.bed chr16.bed chrX.bed > test_somalier_small.bed

# 4. Use the BED file to cut out regions in the BAM files (
# pacbio, this file aligned to GRCh38_no_alt_analysis_set.fasta: https://human-pangenomics.s3.amazonaws.com/submissions/80d00e88-7a92-46d8-88c7-48f1486e11ed--HG002_PACBIO_REVIO/m84039_230117_233243_s1.hifi_reads.default.bam
# ont: s3://ont-open-data/giab_2023.05/analysis/variant_calling/hg002_sup_60x/hg002.haplotagged.bam 
# )  

# Merge the regions to create masking for reference
# Here we will chop the chromosomes of at the end of the regions
bedtools sort -i test_somalier_small.bed > test_somalier_small.sorted.bed
bedtools merge -i test_somalier_small.sorted.bed -d 999999999| awk -v OFS='\t' '{print $1,0,$3}' > chr16X.bed

# 5. Make a reference with the cut out chromosomes,
# 5.1 mask everything which is not in the bed files
# 5.2 fix chr names,
# 5.3 add extra_chr.fa:
# >chr1
# N
# >chr2
# N
# >chr3
# N
# >chr4
# N
# >chr5
# N
# >chr6
# N
# >chr7
# N
# >chr8
# N
# >chr9
# N
# >chr10
# N
# >chr11
# N
# >chr12
# N
# >chr13
# N
# >chr14
# N
# >chr15
# N
# >chr17
# N
# >chr18
# N
# >chr19
# N
# >chr21
# N
# >chr22
# N
# >chrY
# N
# >chrMT
# N

seqtk subseq /data/shared/reference_genomes/GRCh38_no_alt_analysis_set.fasta chr16X.bed|seqtk seq -c -M test_somalier_small.bed | sed '/^>/ s/:.*//' | cat - extra_chr.fa | pigz -p 36 > hg38.test.fa.gz

# Subsample to reduce size
samtools view -@ 36 -s 0.8 -M -L test_somalier_small.bed HG002_PACBIO_REVIO.bam -h -O BAM > hg002_somalier_small_revio.bam
samtools view -@ 36 -s 0.8 -M -L test_somalier_small.bed /media/hgst_12tb_1/felix/HG003_phased.bam -h -O BAM > hg003_somalier_small_revio.bam
samtools view -@ 36 -s 0.8 -M -L test_somalier_small.bed /media/hgst_12tb_1/felix/HG004_phased.bam -h -O BAM > hg004_somalier_small_revio.bam

samtools view -@ 36 -s 0.5 -M -L test_somalier_small.bed hg002.haplotagged.bam -h -O BAM > hg002_somalier_small_ont.bam

# Make fastq
samtools fastq -T \* -@ 36 hg002_somalier_small_revio.bam | pigz -p 36 > hg002_somalier_small_revio.fastq.gz
samtools fastq -T \* -@ 36 hg002_somalier_small_ont.bam |chopper --maxlength 200000| pigz -p 36 > hg002_somalier_small_ont.fastq.gz

echo "reads Revio HG002, HG003, HG004:"
samtools view hg002_somalier_small_revio.bam|wc -l
samtools view hg003_somalier_small_revio.bam|wc -l
samtools view hg004_somalier_small_revio.bam|wc -l
echo "reads ONT HG002, HG003, HG004:"
samtools view hg002_somalier_small_ont.bam|wc -l

# 6. Create a SVDB file from https://zenodo.org/records/11511513/files/CoLoRSdb.GRCh38.v1.0.0.pbsv.jasmine.vcf.gz 
bcftools view -R test_somalier_small.bed  CoLoRSdb.GRCh38.v1.0.0.pbsv.jasmine.vcf.gz --output-type z --output colorsdb.test_data.vcf.gz

# 7. Copy files, and make one copy smaller
cp colorsdb.test_data.vcf.gz /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/reference/colorsdb.test_data.vcf.gz
cp hg38.test.fa.gz /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/reference/hg38.test.fa.gz
cp test_somalier_small.bed /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/reference/test_data.bed

cp hg002_somalier_small_revio.bam /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG002_PacBio_Revio.bam
cp hg003_somalier_small_revio.bam /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG003_PacBio_Revio.bam
cp hg004_somalier_small_revio.bam /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG004_PacBio_Revio.bam
samtools index /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG002_PacBio_Revio.bam
samtools index /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG003_PacBio_Revio.bam
samtools index /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG004_PacBio_Revio.bam
cp hg002_somalier_small_revio.fastq.gz /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG002_PacBio_Revio.fastq.gz
samtools view -@ 36 -s 0.01 -h -O BAM hg002_somalier_small_revio.bam > /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG002_PacBio_Revio_copy.bam
samtools index /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG002_PacBio_Revio_copy.bam

cp hg002_somalier_small_ont.bam /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG002_ONT.bam
samtools index /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG002_ONT.bam
samtools view -@ 36 -s 0.01 -h -O BAM hg002_somalier_small_ont.bam > /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG002_ONT_copy.bam
samtools index /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG002_ONT_copy.bam
cp hg002_somalier_small_ont.fastq.gz /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/testdata/HG002_ONT.fastq.gz

# Copy this file
cp somalier.sh /media/ssd_4tb/felix/projects/genomic-medicine-sweden/test-datasets/make_test_data.sh
