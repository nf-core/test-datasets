# test-datasets: `raredisease`

This branch contains test data to be used for automated testing with the [nf-core/raredisease](https://github.com/nf-core/raredisease) pipeline.

## Content of this repository

`reference/`: background resources needed by tools of raredisease pipeline

`testdata/`: chr20 test resources

`reference/grch38_gnomad_reformated_-r3.1.1-.vcf.gz`: Gnomad vcf file containing entries for the region chr20:90000-92000

`reference/grch38_vcfanno_config_-v0.2-_chr20.toml`: TOML file for small test

`reference/vcfanno_grch38_small_test.tar.gz`: the archived files of grch38_*.{vcf.gz, vcf.gz.tbi} for small test

`reference/genome.ploidy_priors.tsv`: Contains contig ploidy priors for gatk4's DetermineGermlineContigPloidy

`reference/genome.ploidy_model.tar.gz`: tar gzipped directory containing the ploidy model files

`reference/genome.germline_cnv_model.tar.gz`: tar gzipped directory containing the cnv model files

`reference/mobile_elemement_references.tsv`: tsv file with paths to the mobile element locations on chromosome 21

### For Mitochondrial subworkflow

`reference/Homo_sapiens_assembly38_chr20_chrM.fasta`: chr20 and chrM hg38 reference fasta file

`reference/Homo_sapiens_assembly38_chr20_chrM.fasta.fai`: chr20 and chrM hg38 reference index fasta file

`reference/hg38.chrM.fa`: chrM hg38 reference fasta file

`reference/hg38.chrM.fa.fai`: chrM hg38 reference index fasta file

`reference/hg38.chrMshifted8000.fa`: chrM hg38 reference fasta file shifted by 8000 bp

`reference/hg38.chrMshifted8000.fa.fai`: chrM hg38 reference index fasta file shifted by 8000 bp

`reference/control_region_shifted.chrM.interval_list`: Subset of mitochondrial control regions shifted by 8000 bp

`reference/non_control_region.chrM.interval_list`: Subset of mitochondrial non control regions

`reference/ShiftBack.chain`: Used to liftover from shifted regions

`testdata/NA12878_mito_1.fq.gz`: Test fastq file 1 with chr20 and chrM reads

`testdata/NA12878_mito_2.fq.gz`: Test fastq file 2 with chr20 and chrM reads

`testdata/NA12878_sorted_chrM_chr20_rehead_60pdown.cram`: Test alignment file with chrM and chr20 downsampled and compressed 

`testdata/samplesheet_MT.csv`: samplesheet containing fastq from chr20 and chrM

`testdata/samplesheet_2_samples.csv`: samplesheet containing fastq from chr20 and chrM for 2 patients
