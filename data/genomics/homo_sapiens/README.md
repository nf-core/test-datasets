# Test Data

## Table of contents

- [Data Access](#data-access)
- [Reference files](#reference-files)
- [Index files](#index-files)
- [Output data generation](#output-data-generation)
- [Limitations](#limitations)

## Data Access

1. The raw data was retrieved from [this](https://www.ncbi.nlm.nih.gov/bioproject/?term=prjeb39899) project. The two used datasets are [ERR4467723](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=ERR4467723) (tumor) and [ERR4467726](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=ERR4467726) (normal)

    For RNAseq, the raw data was retrieved from 'Genome in a Bottle' sample GM12878 (SRA accession [SRX2900878](https://www.ncbi.nlm.nih.gov/sra/?term=SRX2900878), sequenced on NextSeq 500 with 150bpx2 library).

2. The data was downloaded using the SRA Toolkit with:

    ```bash
    prefetch -v <Acc number>
    sam-dump <Acc number> | samtools view -bS - > <Acc number>.bam
    ```

3. Reads mapping to chr 22 were extracted and converted to `fastq.gz` using [qbic-pipelines/bamtofastq](https://github.com/qbic-pipelines/bamtofastq)

## Determine region covered by reads

1. Visual inspection as to where the reads map to with IGV.

    ```bash
    chr22   16570000        16610000
    ```

    For RNAseq, all the reads and their pair-mates that overlap the SNP sites in the above region were extracted using [VariantBAM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4920121/) and converted to `fastq.gz` using the method described above. Files are available in illumina/fastq/

2. Save length in `genome.bed` 0-40001

## Reference files

### VCF files
Following 'reference' vcf files are generated. All found in igenomes at `s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/`:

- dbsnp_146.hg38
- Mills_and_1000G_gold_standard.indels.hg38
- gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only

1. Downsize all vcf files to only cover the region in chromosome 22 with:

    ```bash
    tabix -l dbsnp_146.hg38.vcf.gz | parallel -j 5 'tabix -h dbsnp_146.hg38.vcf.gz {} > {}.vcf'
    bgzip dbsnp_chr22.vcf
    tabix chr22.vcf
    bcftools filter dbsnp_146.hg38.chr22.vcf.gz -r chr22:16570000-16610000 > region_22/dbsnp_146.hg38.chr22_region.vcf
    bgzip dbsnp_146.hg38.chr22_region.vcf
    tabix dbsnp_146.hg38.chr22_region.vcf.gz

    mv dbsnp_146.hg38.chr22_region.vcf.gz dbsnp_146.hg38.vcf.gz
    mv dbsnp_146.hg38.chr22_region.vcf.gz.tbi dbsnp_146.hg38.vcf.gz.tbi
    ```

2. Manipulated mills & gnomAD file, by changing chr length for chr22 to 40001.

3. The `syntheticvcf_short.vcf.gz` set is a synthetically generated dataset containing ~2000 common variants made for testing polygenic risk scoring. It is expected to score particularly high for type 2 diabetes.

4. justhusky_minimal.vcf.gz and associated files justhusky_minimal.vcf.gz.tbi and justhusky.ped is a subsampled minimal example vcf/ped combination made for testing family-related modules. justhusky_minimal.vcf.gz.tbi was generated with tabix.

### Fasta

As base reference `s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/Chromosomes/chr22.fasta` was used.

```bash
samtools faidx chr22.fasta chr22:16570000-16610000  > genome.fasta
```

The corresponding transcriptome file was extracted:
```bash
gffread -F -w transcriptome.fasta -g genome.fasta genome.gtf
```

An interval list file was prepared from the genome.bed using GATK4:
```bash
gatk BedToIntervalList -I genome.bed -SD genome.dict -O genome.interval_list
```

### GTF/GFF

Downloaded the gtf and gff3 files from Ensembl:

1. Download

    ```bash
    wget http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz
    wget http://ftp.ensembl.org/pub/release-103/gff3/homo_sapiens/Homo_sapiens.GRCh38.103.chromosome.22.gff3.gz
    ```

2. Unzip both with

    ```bash
    gzip -d
    ```

3. Copy `test.bed`, change chromosome name to `22`

    ```bash
    bedtools intersect -a test.bed -b Homo_sapiens.GRCh38.103.chr.gtf -wa -wb > genome_bed.gtf
    ```

4. Remove the first three columns in both files:

    ```bash
    awk '{ $1=""; $2=""; $3=""; print}' genome_bed.gtf > genome.gtf
    ```

5. Change chromosome name to `chr22`
6. Replace spaces with tabs
7. The coordinates in `genome.gtf` were adapted to start from 1, and the last entries that ended in coordinates >40000 were adapted to end at coordinate 40000.

## Index files

### salmon index
The salmon index  (`homo_sapiens/genome/index/salmon`) was created with the following command:

```bash
salmon index -t transcriptome.fasta -k 31 -i salmon
```

## Output data generation

### Sarek pipeline generation

1. Used release 2.7.1 container:
2. Add `publishDir` to all UMI related steps
3. Add to mapping process:

    ```bash
    gatk --java-options -Xmx${task.memory.toGiga()}g SamToFastq --INPUT=${inputFile1} --FASTQ=/dev/stdout --INTERLEAVE=true     --NON_PF=true > ${inputFile1}.fq.gz
    ```

    and `publish` the reads. Un-interleave reads after sarek is run:

    ```bash
    paste - - - - - - - - < test2_umi-consensus.bam.fq.gz | tee >(cut -f 1-4 | tr "\t" "\n" > test2_1.fq) | cut -f 5-8 | tr "\t" "\n" > test2_2.fq
    ```
    Double-check the integrity of the read files with `seqkit sana` and `seqkir pair`

4. Add `publishDir` to HaplotypeCaller process to publish `.g.vcf` files
5. Run sarek with the following command:

```bash
nextflow run nf-core/sarek -r 2.7.1 -profile cfc -c sarek.config \
--input 'testdata_dsl2_chr22.tsv' \
--outdir 'results_fix_bams_again' \
--intervals false  \
--bwa false \
--igenomes_ignore  \
--save_reference \
--fasta 'data/genomics/homo_sapiens/genome/genome.fasta' \
--save_bam_mapped \
--genome custom \
--dict false \
--dbsnp 'data/genomics/homo_sapiens/genome/vcf/dbsnp_146.hg38.vcf.gz' \
--dbsnp_index 'data/genomics/homo_sapiens/genome/vcf/dbsnp_146.hg38.vcf.gz.tbi' \
--known_indels 'data/genomics/homo_sapiens/genome/vcf/mills_and_1000G.indels.vcf.gz' \
--known_indels_index 'data/genomics/homo_sapiens/genome/vcf/mills_and_1000G.indels.vcf.gz.tbi' \
--germline_resource 'data/genomics/homo_sapiens/genome/vcf/gnomAD.r2.1.1.vcf.gz' \
--germline_resource_index 'data/genomics/homo_sapiens/genome/vcf/gnomAD.r2.1.1.vcf.gz.tbi' \
--pon 'data/genomics/homo_sapiens/genome/vcf/mills_and_1000G.indels.vcf.gz' \
--pon_index 'data/genomics/homo_sapiens/genome/vcf/mills_and_1000G.indels.vcf.gz.tbi' \
--tools 'CNVkit,mutect2,freebayes,mpileup,msisensor,cnvkit,strelka,HaplotypeCaller,Manta,tiddit' \
--max_memory 59.GB \
--max_cpus 19 \
-resume \
--generate_gvcf \
```

with the following TSV:

```bash
test	XY	0	testN	1	test_umi_1.fq.gz	test_umi_2.fq.gz
test	XY	1	testT	2	test2_umi_1.fq.gz	test2_umi_2.fq.gz
```

#### GVCF files

1. Set up conda environment:

    ```bash
    conda install -c bioconda gatk4=4.2.0.0
    conda install -c bioconda tabix=1.11
    ```

2. Take the vcf files generated with haplotypecaller within sarek and run commands on both `test.genome.vcf` and `test2.genome.vcf` files:

    ```bash
    gatk IndexFeatureFile -I test.genome.vcf
    bgzip test.genome.vcf
    tabix test.genome.vcf.gz
    ```



#### CRAM files

The cram files were generated with

    ```bash
    samtools view -C -T genome.fasta -o test2.paired_end.recalibrated.sorted.cram test2.paired_end.recalibrated.sorted.bam
    samtools index test2.paired_end.recalibrated.sorted.cram
    ```

#### GATK

`test_test2_paired_mutect2_calls.artifact-prior.tar.gz`:

```
gatk LearnReadOrientationModel -I ..illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.f1r2.tar.gz -O test_test2_paired_mutect2_calls.artifact-prior.tar.gz
```
#### Paired Mutect files

The unfiltered calls are used from the Mutect2 output directory in sarek.



#### GenomicsDB

```
gatk GenomicsDBImport -V ../gvcf/test.genome.vcf --genomicsdb-workspace-path test_genomicsdb -L ../../genome/genome.interval_list
```

The test_genomicsdb has been placed in a tar archive in order to make it easier to download for tests, please remember to untar the directory when using it for tests.
The directory name: /test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/
and the file:
/test_genomicsdb/chr22$1$40001/genomicsdb_meta_dir/genomicsdb_meta_2b25a6c2-cb94-4a4a-9005-acb7c595d322.json change with each run, but the contents of the file and directory will remain the same. Rename them to the above values to keep tests passing.

### 10X genomics scRNA-seq data
10X Genomics (v3) FastQ files covering chr22 are contained in `illumina/10xgenomics`
Data generation:

1. Data was downloaded from https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3
2. A STAR index was generated from `genome.fasta` and `genome.gtf`
3. Reads were aligned with STAR using the generated index:

    ```bash
    STAR --genomeDir star --readFilesIn pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz --runThreadN 6 --outFileNamePrefix test --soloCBwhitelist ../10x_V3_barcode_whitelist.txt --sjdbGTFfile genome/genome.gtf --soloType CB_UMI_Simple --readFilesCommand zcat --soloBarcodeReadLength 28
    ```

4. Extract the readnames: `samtools view testAligned.out.sam | cut -f 1 > readnames.txt`
5. Extract reads mapping to chr22 from fastq files:

    ```bash
    seqtk subseq pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq  readnames.txt  > pbmc_R1.fastq`
    seqtk subseq pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq  readnames.txt  > pbmc_R2.fastq`

    ```

6. Subsample 100 reads

    ```bash
    seqtk sample -s100 pbmc_R1.fastq 100 > test_1.fastq
    seqtk sample -s100 pbmc_R2.fastq 100 > test_2.fastq
    ```

### cooler test dataset

The raw data were downloaded from https://github.com/open2c/cooler/tree/master/tests/data

### PACBIO test dataset
The first 1000 raw reads were extracted from the [public Alzheimer dataset](https://downloads.pacbcloud.com/public/dataset/IsoSeq_sandbox/2020_Alzheimer8M_subset/alz.1perc.subreads.bam) using samtools.
```
samtools view -h alz.1perc.subreads.bam|head -n 1006|samtools view -bh > alz.bam
```
The lima, refine, clustered, singletons and gene models datasets were generated using the isoseq3 framework and TANA collapse.

## Limitations

1. Reads do not cover chromosome 6

### Missing files

1. Single-end reads
2. Methylated bams
3. Unaligned bams
4. Panel of Normals
5. Ploidy files for ASCAT
6. Mappability files for CONTROLFREEC
