# 10x Genomics - 10k Human PBMCs, Multiome v1.0, Chromium X

Human peripheral blood mononuclear cells (PBMCs) of a healthy male donor aged 30-35 were obtained by 10x Genomics from AllCells. Nuclei were isolated using the demonstrated protocol Nuclei Isolation for Single Cell Multiome ATAC + Gene Expression Sequencing (CG000365).

Gene Expression and ATAC libraries were generated from ~16,100 transposed nuclei (10,974 nuclei recovered) as described in the Chromium Next GEM Single Cell Multiome ATAC + Gene Expression Reagent Kits User Guide (CG000338 Rev E) using the Chromium X and sequenced on an Illumina NovaSeq 6000.

The reads from lane 001 of both GEX and ATAC-seq modalities were aligned using cellranger-arc, then reads aligned to chr21 were extracted and subsampled, as described below.

## Download original input data

```bash
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_library.csv
```

## Extracting reads aligned to chr21 and subsampling 

```bash
CHR_NAME="chr21"
NREADS=40000

singularity run -B /scratch/ https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0 samtools view atac_possorted_bam.bam $CHR_NAME | cut -f1  > atac_reads_name
singularity run -B /scratch/ https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0 samtools view gex_possorted_bam.bam $CHR_NAME | cut -f1  > gex_reads_name

singularity run -B /scratch/ https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0 seqkit grep -f atac_reads_name 10k_PBMC_Multiome_nextgem_Chromium_X_atac_S2_L001_R1_001.fastq.gz -o 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_S2_L001_R1_001.fastq.gz
singularity run -B /scratch/ https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0 seqkit grep -f atac_reads_name 10k_PBMC_Multiome_nextgem_Chromium_X_atac_S2_L001_R2_001.fastq.gz -o 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_S2_L001_R2_001.fastq.gz
singularity run -B /scratch/ https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0 seqkit grep -f atac_reads_name 10k_PBMC_Multiome_nextgem_Chromium_X_atac_S2_L001_R3_001.fastq.gz -o 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_S2_L001_R3_001.fastq.gz

singularity run -B /scratch/ https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0 seqkit grep -f gex_reads_name 10k_PBMC_Multiome_nextgem_Chromium_X_gex_S2_L001_R1_001.fastq.gz -o 10k_PBMC_Multiome_nextgem_Chromium_X_gex_chr21_S2_L001_R1_001.fastq.gz
singularity run -B /scratch/ https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0 seqkit grep -f gex_reads_name 10k_PBMC_Multiome_nextgem_Chromium_X_gex_S2_L001_R2_001.fastq.gz -o 10k_PBMC_Multiome_nextgem_Chromium_X_gex_chr21_S2_L001_R2_001.fastq.gz

zcat 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_S2_L001_R1_001.fastq.gz | head -n $NREADS > 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_subsample_S2_L001_R1_001.fastq
gzip 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_subsample_S2_L001_R1_001.fastq
zcat 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_S2_L001_R2_001.fastq.gz | head -n $NREADS > 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_subsample_S2_L001_R2_001.fastq
gzip 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_subsample_S2_L001_R2_001.fastq
zcat 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_S2_L001_R3_001.fastq.gz | head -n $NREADS > 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_subsample_S2_L001_R3_001.fastq
gzip 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr21_subsample_S2_L001_R3_001.fastq

zcat 10k_PBMC_Multiome_nextgem_Chromium_X_gex_chr21_S2_L001_R1_001.fastq.gz | head -n $NREADS > 10k_PBMC_Multiome_nextgem_Chromium_X_gex_chr21_subsample_S2_L001_R1_001.fastq
gzip 10k_PBMC_Multiome_nextgem_Chromium_X_gex_chr21_subsample_S2_L001_R1_001.fastq
zcat 10k_PBMC_Multiome_nextgem_Chromium_X_gex_chr21_S2_L001_R2_001.fastq.gz | head -n $NREADS > 10k_PBMC_Multiome_nextgem_Chromium_X_gex_chr21_subsample_S2_L001_R2_001.fastq
gzip 10k_PBMC_Multiome_nextgem_Chromium_X_gex_chr21_subsample_S2_L001_R2_001.fastq
```
