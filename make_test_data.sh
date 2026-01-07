#!/bin/bash
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
>chrM
N
>chrY
N
EOF

seqtk subseq data/GRCh38_no_alt_analysis_set.fasta tmp/chromosomes_with_sequence.bed | seqtk seq -c -M tmp/test_somalier_small.bed | sed '/^>/ s/:.*//' | cat - tmp/extra_chr.fa | pigz -p 36 > tmp/hg38.test.fa.gz

# Prepare small bam files by subsampling and remapping
function prepare_bam {
  local subsample_to_n_reads=$1
  local minimap_preset=$2
  local in_bam=$3
  local out_bam=$4
  local rg_file="tmp/$(basename $out_bam .bam).rg.txt"

  # Extract read group information to keep it after remapping
  samtools view -H ${in_bam} | grep "^@RG" > ${rg_file}

  # Subsample and remap. We must reheader to keep read group information
  samtools view -M -L tmp/test_somalier_small.bed ${in_bam} -h -O BAM -u -x HP,PS,AS,CC,CG,CP,H1,H2,HI,H0,IH,MC,MD,MQ,NM,SA,TS \
    | samtools reset \
    | samtools sort -n -O SAM \
    | awk -v n="$subsample_to_n_reads" '/^@/ {print; next} ++c <= n {print} c==n {exit}' \
    | samtools fastq -T '*' \
    | minimap2 -a -x ${minimap_preset} -y --secondary=no -Y --MD -t 36 tmp/hg38.test.fa.gz - \
    | samtools view -@ 36 -h -O BAM -u - \
    | samtools reheader  -c "cat - ${rg_file}" - \
    | samtools sort -o ${out_bam}

    samtools index ${out_bam}
}

prepare_bam 1050 map-hifi data/HG002_haplotagged.bam testdata/HG002_PacBio_Revio.bam
prepare_bam 690 map-hifi data/HG003_haplotagged.bam testdata/HG003_PacBio_Revio.bam
prepare_bam 1150 map-hifi data/HG004_haplotagged.bam testdata/HG004_PacBio_Revio.bam
prepare_bam 1200 lr:hq data/hg002.haplotagged.bam testdata/HG002_ONT.bam
prepare_bam 11 map-hifi data/HG002_haplotagged.bam data/HG002_PacBio_Revio_copy.bam
prepare_bam 35 lr:hq data/hg002.haplotagged.bam testdata/HG002_ONT_copy.bam

# Make fastq
samtools fastq -T \* -@ 36 testdata/HG002_PacBio_Revio.bam | pigz -p 36 > testdata/HG002_PacBio_Revio.fastq.gz
samtools fastq -T \* -@ 36 testdata/HG002_ONT.bam |chopper --maxlength 200000| paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | pigz -p 36 > testdata/HG002_ONT.fastq.gz

# 6. Create a SVDB file from https://zenodo.org/records/11511513/files/CoLoRSdb.GRCh38.v1.0.0.pbsv.jasmine.vcf.gz
bcftools view -R tmp/test_somalier_small.bed  data/CoLoRSdb.GRCh38.v1.0.0.pbsv.jasmine.vcf.gz --output-type z --no-version --output tmp/colorsdb.test_data.vcf.gz

# 7. Copy files
cp tmp/colorsdb.test_data.vcf.gz reference/colorsdb.test_data.vcf.gz
cp tmp/hg38.test.fa.gz reference/hg38.test.fa.gz
cp tmp/test_somalier_small.bed reference/test_data.bed

# Print read counts
echo "reads Revio HG002, HG003, HG004:"
samtools view tmp/hg002_somalier_small_revio.bam|wc -l
samtools view tmp/hg003_somalier_small_revio.bam|wc -l
samtools view tmp/hg004_somalier_small_revio.bam|wc -l
echo "reads ONT HG002"
samtools view tmp/hg002_somalier_small_ont.bam|wc -l
echo "Reads in HG002_PacBio_Revio_copy.bam:"
samtools view testdata/HG002_PacBio_Revio_copy.bam | wc -l
echo "Reads in HG002_ONT_copy.bam:"
samtools view testdata/HG002_ONT_copy.bam | wc -l

# 8. Unzip reference file
gunzip -c reference/hg38.test.fa.gz > reference/hg38.test.fa

# 9. Gnomad test data set
bcftools view --regions-file tmp/test_data.bed data/gnomad.genomes.v4.1.sites.chr16.vcf.bgz --output-type b --output tmp/gnomad_chr16.test_data.bcf.gz
bcftools view --regions-file tmp/test_data.bed data/gnomad.genomes.v4.1.sites.chrX.vcf.bgz --output-type b --output tmp/gnomad_chrX.test_data.bcf.gz
bcftools concat tmp/gnomad_chr16.test_data.bcf.gz tmp/gnomad_chrX.test_data.bcf.gz --output-type b --output tmp/gnomad_chr16_chrX.test_data.bcf.gz
docker run -v $(pwd):/data docker.io/fellen31/echtvar:0.2.2 echtvar encode /data/tmp/gnomad_af.zip /data/assets/gnomad.json /data/tmp/gnomad_chr16_chrX.test_data.bcf.gz
cp tmp/gnomad_af.zip reference/gnomad_af.zip