# 60k Human PBMCs Multiplexed, 4 CMOs

Peripheral blood mononuclear cells (PBMCs). Samples originate from two healthy donors (Donor 1 & 2), and from two patients diagnosed with Acute Lymphocytic Leukemia (Donor 3 & 4). Samples from donor 1,2, and 4 consist of peripheral blood mononuclear cells (PBMCs), while the sample from donor 3 consists of Bone Marrow Mononuclear Cells (BMMCs).

Libraries were prepared following the Single Cell V(D)J Reagent Kits User Guide CG000424 (5' v2 HT) and sequenced on Illumina NovaSeq 6000.

## Download original input data
Data is available from 10x:
[Demultiplexing 5’ Immune Profiling Libraries Pooled with Hashtags](https://www.10xgenomics.com/datasets/5-hashing-example-with-tabs-2-standard)

```bash
# Input Files
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/7.0.0/5p_hashing_demux/5p_hashing_demux_gex_fastq.tar
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/7.0.0/5p_hashing_demux/5p_hashing_demux_ab_fastq.tar
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/7.0.0/5p_hashing_demux/5p_hashing_demux_tcr_fastq.tar
curl -O https://cf.10xgenomics.com/samples/cell-vdj/7.0.0/5p_hashing_demux/5p_hashing_demux_bcr_fastq.tar
curl -O https://cf.10xgenomics.com/samples/cell-vdj/7.0.0/5p_hashing_demux/5p_hashing_demux_cmo-set.csv
curl -O https://cf.10xgenomics.com/samples/cell-vdj/7.0.0/5p_hashing_demux/5p_hashing_demux_config.csv
```

## Original `cellranger multi` configs
Since this data consist of both multiplexed samples and immune profiling data, data must be processed as follows:

1. Demultiplex by sample
2. Convert BAM to FASTQ per sample
3. Run cellranger for each sample

See [this](https://www.10xgenomics.com/analysis-guides/demultiplexing-and-analyzing-5’-immune-profiling-libraries-pooled-with-hashtags) analysis guide.

### Sample demultiplexing config (step 1)
This config is to demultiplex by sample. This step requires GEX and CMO data.

```bash
[gene-expression]
reference,/path/to/cellranger/prebuilt/human/refdata-gex-GRCh38-2020-A/
cmo-set,/path/to/wrkdir/demultiplexing/5p_hashing_demux_cmo-set.csv

[libraries]
fastq_id,fastqs,feature_types
PBMC-ALL_60k_universal_HashAB1-4_BL_4tags_Rep1_gex,Path/to/GEX/FASTQs,Gene Expression
PBMC-ALL_60k_universal_HashAB1-4_BL_4tags_Rep1_ab,Path/to/AB/FASTQs,Multiplexing Capture

[samples]
sample_id,cmo_ids
Donor1_healthy,Hash-tag1
Donor2_healthy,Hash-tag2
Donor3_ALLpatient,Hash-tag3
Donor4_ALLpatient,Hash-tag4
```

### Convert BAM to FASTQ (step 2)
For each sample output convert the GEX BAM file to FASTQ.
This step requires the 10x provided bamtofastq binaries which is bundled with cellranger.
Add `bamtofastq`to your `$PATH` like so:

```bash
# Command to change to the directory where the Cell Ranger executable file lives and put it in your $PATH:
cd /path/to/cellranger-7.0.0/
export PATH=${PWD}:$PATH

# Command to put other tools bundled with Cell Ranger in your path:
source /path/to/cellranger-7.0.0/sourceme.bash
```

Run `bamtofastq` per sample:
```bash
mkdir bamtofastq

bamtofastq /path/to/wrkdir/demultiplexing/demultiplexed_samples/outs/per_sample_outs/Donor1_healthy/count/sample_alignments.bam /path/to/wrkdir/bamtofastq/Donor1_healthy
bamtofastq /path/to/wrkdir/demultiplexing/demultiplexed_samples/outs/per_sample_outs/Donor2_healthy/count/sample_alignments.bam /path/to/wrkdir/bamtofastq/Donor2_healthy
bamtofastq /path/to/wrkdir/demultiplexing/demultiplexed_samples/outs/per_sample_outs/Donor3_ALLpatient/count/sample_alignments.bam /path/to/wrkdir/bamtofastq/Donor3_ALLpatient
bamtofastq /path/to/wrkdir/demultiplexing/demultiplexed_samples/outs/per_sample_outs/Donor4_ALLpatient/count/sample_alignments.bam /path/to/wrkdir/bamtofastq/Donor4_ALLpatient
```

Expected output:
```bash
$ ls /path/to/wrkdir/bamtofastq/Donor1_healthy

demultiplexed_samples_0_1_HW7KMDSX3/ demultiplexed_samples_1_1_HW7WGDSX3/
```

Please note that the order of the libraries depends on the order that the libraries are listed in 5p_hashing_demux_config.csv (Step 1).
The GEX library is listed before Multiplexing Capture, and therefore GEX data is found in `demultiplexed_samples_0_1_HW7KMDSX3/`.

### Immune profiling per sample (step 3)
This step must be run for each sample.

```bash
mkdir final-analysis
```

Below are listed configs for each sample.

1. Donor
```bash
[gene-expression]
reference,wrkdir/refdata-gex-GRCh38-2020-A/
force-cells,11638
check-library-compatibility,false

[feature]
reference,wrkdir/final-analysis/feature_ref.csv

[vdj]
reference,wrkdir/final-analysis/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0

[libraries]
fastq_id,fastqs,feature_types
bamtofastq,wrkdir/bamtofastq/Donor1_healthy/demultiplexed_samples_0_1_HW7KMDSX3,Gene Expression
PBMC_60k_10x_4tags_Rep1_Ab,wrkdir/fastqs/ab_fastq,Antibody Capture
PBMC_60k_10x_4tags_Rep1_TCR,wrkdir/fastqs/tcr_fastq,VDJ-T
PBMC_60k_10x_4tags_Rep1_BCR,wrkdir/fastqs/bcr_fastq,VDJ-B
```

2. Donor 2
```bash
[gene-expression]
reference,wrkdir/refdata-gex-GRCh38-2020-A/
force-cells,8873
check-library-compatibility,false

[feature]
reference,wrkdir/final-analysis/feature_ref.csv

[vdj]
reference,wrkdir/final-analysis/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0

[libraries]
fastq_id,fastqs,feature_types
bamtofastq,wrkdir/bamtofastq/Donor2_healthy/demultiplexed_samples_0_1_HW7KMDSX3,Gene Expression
PBMC_60k_10x_4tags_Rep1_Ab,wrkdir/fastqs/ab_fastq,Antibody Capture
PBMC_60k_10x_4tags_Rep1_TCR,wrkdir/fastqs/tcr_fastq,VDJ-T
PBMC_60k_10x_4tags_Rep1_BCR,wrkdir/fastqs/bcr_fastq,VDJ-B
```

3. Donor 3
```bash
[gene-expression]
reference,wrkdir/refdata-gex-GRCh38-2020-A/
force-cells,9024
check-library-compatibility,false

[feature]
reference,wrkdir/final-analysis/feature_ref.csv

[vdj]
reference,wrkdir/final-analysis/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0

[libraries]
fastq_id,fastqs,feature_types
bamtofastq,wrkdir/bamtofastq/Donor1_healthy/demultiplexed_samples_0_1_HW7KMDSX3,Gene Expression
PBMC_60k_10x_4tags_Rep1_Ab,wrkdir/fastqs/ab_fastq,Antibody Capture
PBMC_60k_10x_4tags_Rep1_TCR,wrkdir/fastqs/tcr_fastq,VDJ-T
PBMC_60k_10x_4tags_Rep1_BCR,wrkdir/fastqs/bcr_fastq,VDJ-B
```

4. Donor 4
```bash
[gene-expression]
reference,path/to/cellranger/prebuilt/human/refdata-gex-GRCh38-2020-A/
force-cells,7032
check-library-compatibility,false

[feature]
reference,/path/to/feature_ref.csv

[vdj]
reference,/path/to/vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0

[libraries]
fastq_id,fastqs,feature_types
bamtofastq,/path/to/bamtofastq_outs_v1.4/Donor4_ALLpatient/PBMC_60k_10x_4tags_Rep1_step1_v7_0_1_HMW77DSX3,Gene Expression
PBMC_60k_10x_4tags_Rep1_Ab,/path/to/Ab-fastqs,Antibody Capture
PBMC_60k_10x_4tags_Rep1_TCR,/path/to/Ab-fastqs/TCR,VDJ-T
PBMC_60k_10x_4tags_Rep1_BCR,/path/to/Ab-fastqs/BCR/,VDJ-B
```

Run cellranger multi:
```bash
cellranger multi --id=Donor1-final --csv=Donor1_config.csv
cellranger multi --id=Donor2-final --csv=Donor2_config.csv
cellranger multi --id=Donor3-final --csv=Donor3_config.csv
cellranger multi --id=Donor4-final --csv=Donor4_config.csv
```

## Cell Multiplexing Data 

The cell multiplexing libraries require more generous subsampling to yield sufficient reads for successful processing. 
Test data use the first 3,500,000 lines of each of the lane 1 FASTQs.

```bash
$ zcat ../SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R1_001.fastq.gz | head -n 3500000 | gzip -c > subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R1_001.fastq.gz
$ zcat ../SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R2_001.fastq.gz | head -n 3500000 | gzip -c > subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R2_001.fastq.gz
```
