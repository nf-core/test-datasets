# Modules Test Data

This branch of the `nf-core/test-datasets` repository contains all data used for the individual module tests.
There are two main directories: `genomics` and `delete_me`. The former contains all datasets for genomics tools while the latter contains temporary datasets that will be deleted as better data gets available.

### delete_me

The `delete_me` folder does not adhere to a defined structure as data in this folder should be directly as fast as possible, whenever a more suitable dataset is found that can be added to the `genomics` folder.

### genomics

The genomics folder contains subfolders for all organisms for which test data is available. At the moment, there are two organisms available:
    *homo_sapiens
    * sarscov2

The two folders are structured in a similar way, with any genome-specific files in `genome` (e.g. fasta, gtf, ...) and technology specific raw-data files
in the `illumina` and `nanopore` subfolders.
It is currently organised in `genomics` and `generic`. The former contains all typical data required for genomics modules, such as fasta, fastq and bam files. Every folder in `genomics` corresponds to a single organisms. Any other data is stored in `generic`. This contains files that currently cannot be associated to a genomics category, but also depreciated files which will be removed in the future and exchanged by files in `genomics`. For every data file, a short description about how this file was generated is available either in this description or in the respective subfolder.

If you cannot find suitable test data on this repository, please contact us on the [nf-core Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)). The goal is to create a suitable, small test dataset that should be usable with the available test data and if possible be generated with modules available from `nf-core/modules`. If creating that test data is difficult but you want to add the module first, it is also possible to add a small subset to the `delete_me` folder to get your module tests running, and then add proper test data afterwards. This should be discussed on slack. In order to add test data. For a short description of the workflow for adding new data, have a look at [here](docs/ADD_NEW_DATA.md)

## Data Description

### genomics

* sarscov2
  * alignment
    * 'all_sites.fas'
    * 'informative_sites.fas'
  * bed
    * 'test.bed': exemplary bed file for the MT192765.1 genome (fasta/test_genomic.fasta)
    * 'test2.bed': slightly modified copy of the above file
    * 'test.bed.gz': gzipped version
    * 'baits.bed'
    * 'test.bed12'
  * db
    * 'kraken2': kraken2 DB
    * 'kraken2.tar.gz': kraken2 DB archive
  * genome
    * 'genome.fasta': MT192765.1 genomem including (GCA_011545545.1_ASM1154554v1)
    * 'genome.dict': GATK dict for 'genome.fasta'
    * 'genome.fasta.fai': fasta index for 'genome.fasta'
    * 'transcriptome.fasta': coding sequencing from MT192765.1 genome (transcripts)
    * 'transcriptome.paf': PAF file for MT192765.1  genome
    * 'genome.gtf': GTF for MT192765.1 genome
    * 'genome.gff3': GFF for MT192765.1 genome
    * 'genome.gff3.gz': bgzipped-version
    * 'genome.sizes': genome size for the MT192765.1 genome

  * illumina
    * bam
      * 'test.paired_end.{,methylated}.bam': sarscov2 sequencing reads aligned against test_genomic.fasta using minimap2
      * 'test.paired_end.{,methylated}.sorted.bam': sorted version of the above bam file
      * 'test.paired_end.{,methylated}.bam.sorted.bam.bai': bam index for the sorted bam file
      * 'test.single_end.bam': alignment (unsorted) of the 'test_1.fastq.gz' reads against test_genomic.fasta using minimap2
      * 'test.unaligned.bam': unmapped BAM file created from 'test_1.fastq.gz' using GATK4 SamToFastq
    * bedgraph
      * 'test.bedgraph'
    * bigwig
      * 'test.bw'
    * deeptools
      * 'test.computeMatrix.mat.gz': matrix generated with deeptools computeMatrix using 'test.bw'
    * fasta
      * 'contigs.fasta': sarscov2 contigs obtained running SPAdes `--rnaviral`on sample1 of the [nf-core/viralrecon tests-dataset](https://github.com/nf-core/test-datasets/tree/viralrecon/illumina/amplicon)
      * 'scaffolds.fasta': sarscov2 scaffolds obtained running SPAdes `--rnaviral`on sample1 of the [nf-core/viralrecon test-dataset](https://github.com/nf-core/test-datasets/tree/viralrecon/illumina/amplicon)
    * fastq
      * 'test_{1,2}.fastq.gz' sarscov2 paired-end sequencing reads
      * 'test_{1,2}.2.fastq.gz‘: copies of the above reads
      * 'test.methylated_{1,2}.fastq.gz' sarscov2 paired-end bisulfite sequencing reads (generated with [Sherman](https://github.com/FelixKrueger/Sherman))
    * gatk
      * 'test.baserecalibrator.table': Recalibration table generated with gatk4 BaseRecalibrator from 'test_paired_end.sorted.bam', using 'test.vcf.gz' as known sites.
    * vcf
      * 'test.vcf', 'test2.vcf': generated from 'test_paired_end.sorted.bam' using bcftools mpileup, call and filter
      * 'test3.vcf': generated from 'test_single_end.sorted.bam' using bcftools mpileup, call and filter
      * 'test2.targets.tsv.gz' from 'test2.vcf.gz' using bcftools query  and bgzip
      * '*.gz': generated from VCF files using bgzip
      * '.tbi': generated from '.vcf.gz' files using `tabix -p vcf -f <file>`
    * wig
      * 'test.wig.gz'
  * nanopore
    * bam
      * 'test.sorted.bam'
      * 'test.sorted.bam.bai'
    * fast5
      * 'fast5.tar.gz': compressed `fast5` folder with the following structure: `fast5/100read_fast5/*.fast5`
    * fastq
      * 'test.fastq.gz'
    * sequencing_summary
      * 'test.sequencing_summary.txt'

* homo_sapiens
  * genome
    * vcf
      * dbsnp: DBSnp file downsampled based on reference position
      * gnomAD: gnomAD file downsampled based on reference position
      * mills_and_1000G: Indels file downsampled based on reference position
    * dict: Sequence dictionary corresponding to `fasta`
    * genome.fasta: Reference fasta based on chr22:16570000-16610000
    * transcriptome.fasta: Reference transcriptome based on `genome.fasta`
    * gff3: Encode GFF3 file downsampled based on reference position
    * gtf: Encode GTF file downsampled based on reference position
    * sizes
    * .bed
    * index
      * salmon: salmon index created with `transcriptome.fasta`
  * illumina
    * 10xgenomics
      * test_{1,2}.fastq.gz: 10X Genomics V3 fastq files from human PBMCs
    * bam:
      * test.paired_end.sorted: Mapped, and sorted reads based on `test{,.umi}_{1,2}` (normal)
      * test.paired_end.markduplicates.sorted: Mapped, sorted, and duplicate marked reads based on `test{,.umi}_{1,2}` (normal)
      * test.paired_end.recalibrated.sorted: Mapped, sorted, duplicate marked, and recalibrated reads based on `test{,.umi}_{1,2}` (normal)
      * test2.paired_end.sorted: Mapped, and sorted reads based on `test2{,.umi}_{1,2}` (tumor)
      * test2.paired_end.markduplicates.sorted: Mapped, sorted, and duplicate marked reads based on `test2{,.umi}_{1,2}` (tumor)
      * test2.paired_end.recalibrated.sorted: Mapped, sorted, duplicate marked, and recalibrated reads based on `test2{,.umi}_{1,2}` (tumor)
      * umi:
        * test.paired_end.umi_*: Files base on  `test.umi_{1,2}` (normal)
        * test2.paired_end.umi_*: Files base on  `test2.umi_{1,2}` (tumor)
    * fastq:
      * test_{1,2}: reads corresponding to normal sample
      * test.umi_{1,2}: UMI tagged reads corresponding to normal sample
      * test2_{1,2}: reads corresponding to tumor sample
      * test2.umi_{1,2}: UMI tagged reads corresponding to tumor sample
    * gatk:
      * test: Recalibration table corresponding to `test{,.umi}_{1,2}` (normal) reads
      * test2: Recalibration table corresponding to `test2{,.umi}_{1,2}` (tumor) reads
    * gvcf:
      * test.genome.vcf: Genome vcf corresponding to `test{,.umi}_{1,2}` (normal) reads
      * test2.genome.vcf: Genome vcf corresponding to `test2{,.umi}_{1,2}` (tumor) reads
      * test{,2}.genome.vcf.gz: Bgzipped file based on `test{,2}.genome.vcf` file
      * test{,2}.genome.vcf.gz.tbi: Tbi index based on `test{,2}.genome.vcf.gz` file
      * test{,2}.genome.vcf.idx: Index feature file based on `test{,2}.genome.vcf` file
    * vcf:
      * test.rnaseq.vcf: RNAseq vcf corresponding to `test.rnaseq_{1,2}` reads
    * yak:
      * test.yak: Yak kmer index of 1000 of paternal paired-end reads from the GIAB Ashkenazim trio [RM8392](https://www-s.nist.gov/srmors/view_detail.cfm?srm=8392). These reads were selected from D2_S1_L001_R{1,2}_001.fastq.gz and D2_S1_L001_R{1,2}_002.fastq.gz so that they map to `pacbio/fastq/test_hifi.fastq.gz`.
      * test2.yak: Yak kmer index of 1000 of maternal reads from the GIAB Ashkenazim trio [RM8392](https://www-s.nist.gov/srmors/view_detail.cfm?srm=8392). These reads were selected from D3_S1_L001_R{1,2}_001.fastq.gz and D3_S1_L001_R{1,2}_001.fastq.gz so that they map to `pacbio/fastq/test_hifi.fastq.gz`.
  * pacbio:
    * fastq:
      * test_hifi.fastq.gz: Reads mapping to a randomly selected contig from the whole genome assembly by [Cheng et al., 2021](https://www.nature.com/articles/s41592-020-01056-5) of the child of the GIAB Ashkenazim trio [RM8392](https://www-s.nist.gov/srmors/view_detail.cfm?srm=8392). The reads were taken from [SRR10382244](https://www.ncbi.nlm.nih.gov/sra/?term=SRR10382244).
  * cooler:
    * cload:
      * hg19:
        * hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz, hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz.px2: hg19 pairix test file and its index file.
        * hg19.GM12878-MboI.pairs.subsample.sorted.possrt.txt.gz, hg19.GM12878-MboI.pairs.subsample.sorted.possrt.txt.gz.tbi: hg19 tabix test file and its index file.
        * hg19.sample1.pairs: hg19 pair test file.
        * hg19.chrom.sizes: hg19 chromosome sizes. Downloaded from [goldenpath](http://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes)
    * merge:
      * toy:
        * toy.symm.upper.2.cool, toy.symm.upper.2.cp2.cool: test file for cooler_merge. Downloaded from [open2c/cooler](https://github.com/open2c/cooler/master/tests/data/toy.symm.upper.2.cool)

### Uncategorized

* e_coli_k12_16s.fna: E. coli K-12 16S rRNA
* bac.16S_rRNA.hmm: Bacterial 16S HMM file
