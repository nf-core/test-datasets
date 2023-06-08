# Modules Test Data

This branch of the `nf-core/test-datasets` repository contains all data used for the individual module tests.
There are three main directories: `generic`, `genomics` and `delete_me`. The first contains generic files, the second contains all datasets for genomics tools while the latter contains temporary datasets that will be deleted as better data gets available.

## Adding New Data

If you cannot find suitable test data on this repository, please contact us on the [nf-core Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)). The goal is to create a suitable, small test dataset that should be usable with the available test data and if possible be generated with modules available from `nf-core/modules`. If creating that test data is difficult but you want to add the module first, it is also possible to add a small subset to the `delete_me` folder to get your module tests running, and then add proper test data afterwards. This should be discussed on slack. In order to add test data. For a short description of the workflow for adding new data, have a look at [here](docs/ADD_NEW_DATA.md).

### delete_me

The `delete_me` folder does not adhere to a defined structure as data in this folder should be delete as fast as possible, whenever a more suitable dataset is found that can be added to any other folder.

### generic

The `generic` folder contains generic files that currently cannot be associated to a genomics category. They are organised by their respective file extension. Also, it contains depreciated files which will be removed in the future and exchanged by files in `genomics`.

### genomics

The genomics folder contains subfolders for all organisms for which test data is available. At the moment, there are these organisms available:

- bacteroides_fragilis
- homo_sapiens
- sarscov2
- galaxea_fascicularis

The folders are structured in a similar way, with any genome-specific files in `genome` (e.g. fasta, gtf, ...) and technology specific raw-data files in the `10xgenomics`, `illumina`, `nanopore`, `pacbio`, `hic` and `cooler` subfolders whenever available.
`Genomics` contains all typical data required for genomics modules, such as fasta, fastq and bam files. Every folder in `genomics` corresponds to a single organism. For every data file, a short description about how this file was generated is available either in this description or in the respective subfolder.

### imaging

The imaging folder contains data related to the analysis of highly-multiplexed imaging data, such as image preprocessing, cell segmentation, signal quantification and others. The files are organized by their respective data type.

- h5: image files in HDF5 format
- ilp: ilastik specific project files
- ome-tiff: OME-TIFF image files
- tiff: TIFF image files

### pangenomics

The pangenomics folder contains subfolders for all organisms for which test data is available. At the moment, there is one organism available:

- homo_sapiens

The folder is structured in the following way: Any nonspecific-pangenome file is located in `pangenome` (e.g. PAF, GFA, ...) and software specific binary files in the `odgi` subfolder. `Pangenomics` contains all typical data required for pangenomics modules, such as PAF, GFA files including the binary formats ODGI, and LAY. Every folder in `pangenomics` corresponds to a single organism. For every data file, a short description about how this file was generated is available either in this description or in the respective subfolder. All files in the `pangenomics` folder originates from a [PGGB](https://github.com/pangenome/pggb) run using the [HLA V-352962 gene FASTA](https://github.com/pangenome/pggb/blob/master/data/HLA/V-352962.fa.gz).

## Data Description

### genomics

- sarscov2

  - alignment
    - last
      - 'contigs.genome.maf.gz': alignment of 'contigs.fasta' to 'scaffolds.fasta', in MAF format.
      - 'contigs.genome.par': alignment parameters for comparing 'contigs.fasta' to 'scaffolds.fasta' with LAST.
      - 'lastdb.tar.gz ': 'scaffolds.fasta' index archive for the LAST aligner.
    - 'all_sites.fas'
    - 'informative_sites.fas'
  - bed
    - 'test.bed': exemplary bed file for the MT192765.1 genome (fasta/test_genomic.fasta)
    - 'test2.bed': slightly modified copy of the above file
    - 'test.bed.gz': gzipped version
    - 'baits.bed'
    - 'test.bed12'
    - 'bed6alt.as': AutoSQL file to describe an alternative bed 3+3 format
  - cnn
    - 'reference.cnn': exemplary copy-number reference file for MT192765.1 genome
  - db
    - 'kaiju': Kaiju DB. Database created from ORF1ab polyprotein UNJ12943.1 and taxonomic ID 2697049
    - 'kraken2': kraken2 DB
    - 'kraken2.tar.gz': kraken2 DB archive
    - 'kraken2_bracken': kraken2 & Bracken DB
    - 'kraken2_bracken.tar.gz': kraken2 & Bracken DB archive
    - kofamscan: kofamscan DB files
      - 'ko_list.gz': compressed text file list of KO terms
      - 'profiles.tar.gz': directory archive with HMMER profiles
    - 'metamaps.tar.gz': metamaps DB archive
    - maltextract
      - 'taxon_list.txt': text file of list NCBI sarscov2 species IDs primarily used for MaltExtract
      - 'ncbi_taxmap.zip': mini-NCBI taxonomy map prmiarily used for MaltExtract
    - 'mmseqs.tar.gz': mmseqs DB archive
  - genome
    - 'genome.fasta': MT192765.1 genome including (GCA_011545545.1_ASM1154554v1)
    - 'genome.fasta.gz': bgzipped version of 'genome.fasta'
    - 'genome.dict': GATK dict for 'genome.fasta'
    - 'genome_strtablefile.zip': An StrTableFile zip folder for 'genome.fasta'
    - 'genome.fasta.fai': fasta index for 'genome.fasta'
    - 'transcriptome.fasta': coding sequencing from MT192765.1 genome (transcripts)
    - 'transcriptome.paf': PAF file for MT192765.1 genome
    - 'genome.gtf': GTF for MT192765.1 genome
    - 'genome.gff3': GFF for MT192765.1 genome
    - 'genome.gff3.gz': bgzipped-version
    - 'genome.sizes': genome size for the MT192765.1 genome
    - graphtyper: files for testing graphtyper‚
      - regions.txt: chromosome names and positions for MT192765.1 genome
    - PRG_test: zipped directory to build a test Population Reference Graph‚
  - illumina
    - bam
      - 'test.paired_end.{,methylated}.bam': sarscov2 sequencing reads aligned against test_genomic.fasta using minimap2
      - 'test.paired_end.{,methylated}.sorted.bam': sorted version of the above bam file
      - 'test.paired_end.{,methylated}.bam.sorted.bam.bai': bam index for the sorted bam file
      - 'test.paired_end.name.sorted.bam': Paired-end bam file sorted by name
      - 'test.single_end.bam': alignment (unsorted) of the 'test_1.fastq.gz' reads against test_genomic.fasta using minimap2
      - 'test.unaligned.bam': unmapped BAM file created from 'test_1.fastq.gz' using GATK4 SamToFastq
      - 'read_group_settings.txt': a simple textfile containing the readgroup settings '1 paired' primarily used for the atlas/splitmerge module
      - 'purecn_ex1.bam': Example BAM file used to test PureCN/coverage
      - 'purecn_ex1.bam.bai': BAM index for 'purecn_ex1.bam'
      - 'purecn_ex1_intervals.txt': File containing genomic intervals to be used for testing PureCN/coverage
      - 'purecn_ex1_normal.txt.gz': Example normal coverage file used to test PureCN/normaldb
      - 'purecn_ex2_normal.txt.gz': Example normal coverage file used to test PureCN/normaldb
      - 'purecn_normalpanel.vcf.gz': Example normal VCF to be used for testing PureCN/normaldb
      - 'purecn_normalpanel.vcf.gz.tbi': Index file for 'purecn_normalpanel.vcf.gz'
    - bedgraph
      - 'test.bedgraph'
    - bigwig
      - 'test.bw'
    - deeptools
      - 'test.computeMatrix.mat.gz': matrix generated with deeptools computeMatrix using 'test.bw'
    - fasta
      - 'contigs.fasta': sarscov2 contigs obtained running SPAdes `--rnaviral`on sample1 of the [nf-core/viralrecon tests-dataset](https://github.com/nf-core/test-datasets/tree/viralrecon/illumina/amplicon)
      - 'scaffolds.fasta': sarscov2 scaffolds obtained running SPAdes `--rnaviral`on sample1 of the [nf-core/viralrecon test-dataset](https://github.com/nf-core/test-datasets/tree/viralrecon/illumina/amplicon)
    - fastq
      - 'test\_{1,2}.fastq.gz' sarscov2 paired-end sequencing reads
      - 'test_interleaved.fastq.gz': Interleaved version of the above
      - 'test\_{1,2}.2.fastq.gz‘: copies of the above reads
      - 'test.methylated\_{1,2}.fastq.gz' sarscov2 paired-end bisulfite sequencing reads (generated with [Sherman](https://github.com/FelixKrueger/Sherman))
    - gatk
      - 'test.baserecalibrator.table': Recalibration table generated with gatk4 BaseRecalibrator from 'test_paired_end.sorted.bam', using 'test.vcf.gz' as known sites.
      - 'test_paired_end_sorted_dragstrmodel.txt': The DRAGEN STR model of 'test_paired_end.sorted.bam'.
    - gfa
      - 'assembly.gfa': assembly in Graphical Fragment Assembly (GFA) 1.0 format
      - 'assembly.gfa.bgz': compressed with block-compressed GZIP (BGZF)
      - 'assembly.gfa.gz': compressed with GZIP
      - 'assembly.gfa.zst': compressed with Zstandard (zstd)
    - sra
      - `SRR13255544.tar.gz`: Tar archive containing SRA file obtained from SRR13255544.
      - `SRR11140744.tar.gz`: Tar archive containing SRA file obtained from SRR11140744.
    - vcf
      - 'test.vcf', 'test2.vcf': generated from 'test_paired_end.sorted.bam' using bcftools mpileup, call and filter
      - 'test3.vcf': generated from 'test_single_end.sorted.bam' using bcftools mpileup, call and filter
      - 'test2.targets.tsv.gz' from 'test2.vcf.gz' using bcftools query and bgzip
      - 'sv_query.vcf.gz': a VCF file containing structural variants in chromosome 22
      - 'sv_query.vcf.gz.tbi': The index of the sv_query.vcf.gz file
      - '\*.gz': generated from VCF files using bgzip
      - '.tbi': generated from '.vcf.gz' files using `tabix -p vcf -f <file>`
      - ped
        - 'justhusky_minimal.vcf.gz': minimal combination example of VCF/PED file
        - 'justhusky.ped': minimal combination example of VCF/PED file
        - '.tbi': generated from '.vcf.gz' files using `tabix -p vcf -f <file>`
    - wig
      - 'test.wig.gz'
    - picard
      - 'test.single_end.bam.readlist.txt': text file of a list of two read IDs primarily for picard FilterSamReads
  - nanopore
    - bam
      - 'test.sorted.bam'
      - 'test.sorted.bam.bai'
    - fast5
      - 'fast5.tar.gz': compressed `fast5` folder with the following structure: `fast5/100read_fast5/*.fast5`
    - fastq
      - 'test.fastq.gz'
      - 'test_2.fastq.gz'
    - sequencing_summary
      - 'test.sequencing_summary.txt'
  - metagenome
    - 'test_1.kraken2.reads.txt': kraken classification of each input read of test file `test_1.fastq.gz`
    - 'test_1.kraken2.report.txt': kraken report after classification of test file `test_1.fastq.gz`
    - 'krona_taxonomy.tab': sars-cov-2 taxonomy tree only extracted from taxonomy.tab database created with ktUpdateTaxonomy
    - 'seqid2taxid.map': taxonomy mapping file of the SARS-CoV2 genome genbank ID with NCBI taxonomy ID, originally generated for KrakenUniq
    - 'nodes_dmp': file including sars-cov-2 taxonomy nodes, originally originated for Kaiju
    - 'names_dmp': file with sars-cov-2 taxonomy names, originally generated for Kaiju

- homo_sapiens
  - 10xgenomics
    - cellranger
      - test_10x_10k_pbmc_5fb_fastq_\{1,2\}_gz: 5' V2 Feature Barcode FASTQs from 10k PBMC data
      - test_10x_10k_pbmc_5gex_fastq_\{1,2\}_gz: 5' V2 gene expression FASTQs from 10k PBMC data
      - test_10x_10k_pbmc_b_fastq_\{1,2\}_gz: 5' V2 B-cell FASTQs from 10k PBMC data
      - test_10x_10k_pbmc_t_fastq_\{1,2\}_gz: 5' V2 T-cell FASTQs from 10K PBMC data
      - test_10x_10k_pbmc_feature_ref_csv: Feature Barcode reference for the 10k PBMC
      - test_10x_10k_pbmc_cmo_cmo_fastq_\{1,2\}_gz: 3' V3 Cell Multiplexing FASTQs from 10k PBMC data with Cell Multiplexing
      - test_10x_10k_pbmc_cmo_gex\{1,2\}_fastq_\{1,2\}_gz: 3' V3 gene expression FASTQs from 10k PBMC data with Cell Multiplexing
      - test_10x_10k_pbmc_cmo_feature_ref_csv: Feature Barcode reference for the 10k PBMC data with Cell Multiplexing
      - test_10x_5k_cmvpos_tcells_ab_fastq_\{1,2\}_gz: Antibody Capture FASTQs from the 5k CMV+ T-cell dataset
      - test_10x_5k_cmvpos_tcells_gex1_fastq_\{1,2\}_gz: Gene expression FASTQs from the 5k CMV+ T-cell dataset
      - test_10x_5k_cmvpos_tcells_vdj_fastq_\{1,2\}_gz: V(D)J FASTQs from the 5k CMV+ T-cell dataset
      - test_10x_5k_cmvpos_tcells_feature_ref_csv: Feature Barcode reference for the 5k CMV+ T-cell dataset
      - test_10x_vdj_ref_json: JSON from version 5 of the 10X-curated human V(D)J reference using GRCh38 and Ensembl 94
      - test_10x_vdj_ref_fasta: FASTA file of V(D)J sequences from version 5 of the 10X-curated human V(D)J reference using GRCh38 and Ensembl 94
      - test_10x_vdj_ref_suppfasta: Supplemental FASTA file from version 5 of the 10X-curated human V(D)J reference using GRCh38 and Ensembl 94
    - cellranger-atac
      - test_scATAC_S1_L001_I1_001.fastq.gz: Dual index i7 read (8bp) of a downsamples version of the cellranger-atac-tiny-bcl-simple-1.0.0.csv data (chr1).
      - test_scATAC_S1_L001_R\{1,3\}_001.fastq.gz: Read 1 and 2 of a downsamples version of the cellranger-atac-tiny-bcl-simple-1.0.0.csv data (chr1).
      - test_scATAC_S1_L001_R2_001.fastq.gz: Dual index i5 read (16 bp) of a downsamples version of the cellranger-atac-tiny-bcl-simple-1.0.0.csv data (chr1).
  - genome
    - BUSCO
      - 'chr22_odb10.tar.gz': BUSCO database 'primates_odb10.2021-02-19.tar.gz' purged of entries not matching 'genome.fasta'.
    - chr1: directory for reference files using chr1 rather than 22, used for cellranger-atac
      - genome.fasta.gz
      - genome.gtf
    - chr21: directory for reference files using chr21 rather than 22, used for most gatk4 testing
      - sequence: directory containing fasta, fai, dict and several other indexes for chr21 including:
        - 'genome_sdf.tar.gz': The SDF (RTG Sequence Data File) folder of the reference genome
        - .{1-4,rev.1-2}.bt2
        - .amb
        - .ann
        - .bwt
        - .pac
        - .sa
        - .interval_list
        - .bed
        - .gtf
        - .cnn: copy number reference file for chr_21
        - .snp: Eigenstrat snp file of 1240k snps on chr 21
      - germlineresources: directory containing several germline resource vcfs and tbis, including:
        - 1000G_omni2.5
        - 1000G_phase1.snps
        - axiom_exome_plus.genotypes.all_populations.poly
        - dbsnp_138.hg38
        - gnomAD.r2.1.1
        - haplotype_map
        - hapmap_3.3.hg38
        - hapmap_3.pop_stratified_chr21
        - mills_and_1000G.indels
        - dbNSFP4.1a.21: Created from dbNSFP database. Chromosome 21 extracted from main file and posterior subsampling of first 100K lines.
      - chromosomes.tar.gz: compressed directory containing the fasta genome file renamed to chr21 (needed for ControlFREEC)
    - vcf
      - dbsnp: DBSnp file downsampled based on reference position
      - gnomAD: gnomAD file downsampled based on reference position
      - mills_and_1000G: Indels file downsampled based on reference position
      - vcfanno
        - 'vcfanno_grch38_module_test.tar.gz': exac.vcf.gz + exac.vcf.gz.tbi and they're reference ExAC vcf used to query
        - 'vcfanno.toml': configuration file for vcfanno to operate
    - tsv
      - functional_genomics.counts.tsv : functional genomics count table for CNV correction
      - library_functional_genomics.tsv : functional genomics library for CNV correction
    - genome_config.json: json config file for cellranger-atac or cellranger-arc  
    - genome.ploidy_model.tar.gz: tar gzipped directory containing the ploidy model files
    - genome.ploidy_calls.tar.gz: tar gzipped directory containing the ploidy call files
    - genome.germline_cnv_model.tar.gz: tar gzipped directory containing the cnv model files
    - genome.germline_cnv_calls.tar.gz: tar gzipped directory containing the cnv calls files
    - dict: Sequence dictionary corresponding to `fasta`
    - genome.fasta: Reference fasta based on chr22:16570000-16610000
    - genome.fasta.gz: bgzipped version of 'genome.fasta'
    - genome.fasta.gz.fai: index file for 'genome.fasta.gz'
    - genome.fasta.gz.gzi: index file for 'genome.fasta.gz'
    - genome2.fasta: Reference fasta based on chr22:16600000-16800000
    - genome_motifs.txt: TF motifs used for cellranger-atac
    - transcriptome.fasta: Reference transcriptome based on `genome.fasta`
    - gff3: Encode GFF3 file downsampled based on reference position
    - gtf: Encode GTF file downsampled based on reference position
    - sizes
    - .bed
    - multi_intervals.bed: Contains the interval from `interval.list` split into two parts
    - blacklist_intervals.bed: Contains the intervals of problematic regions of the genome
    - annotated_intervals.tsv: Contains the intervals of the genome (annotated with gc-content) excluding problematic regions
    - filtered_intervals.interval_list: Contains the intervals of the genome that contain at least one read hit
    - ploidy_priors.tsv: Contains contig ploidy priors for gatk4's DetermineGermlineContigPloidy
    - preprocessed_intervals.counts.tsv: Contains the intervals of the genome excluding problematic regions and the respective read counts
    - preprocessed_intervals.interval_list: Contains the intervals of the genome excluding problematic regions
    - index
      - salmon: salmon index created with `transcriptome.fasta`
    - vep.tar.gz: Compressed VEP cache containing info.txt and synonyms of chr22 only. No annotations included.
  - illumina

    - bam:
      - test.paired*end.sorted: Mapped, and sorted reads based on `test{,.umi}*{1,2}` (normal)
      - test.paired*end.markduplicates.sorted: Mapped, sorted, and duplicate marked reads based on `test{,.umi}*{1,2}` (normal)
      - test.paired*end.recalibrated.sorted: Mapped, sorted, duplicate marked, and recalibrated reads based on `test{,.umi}*{1,2}` (normal)
      - test2.paired*end.sorted: Mapped, and sorted reads based on `test2{,.umi}*{1,2}` (tumor)
      - test2.paired*end.markduplicates.sorted: Mapped, sorted, and duplicate marked reads based on `test2{,.umi}*{1,2}` (tumor)
      - test2.paired*end.recalibrated.sorted: Mapped, sorted, duplicate marked, and recalibrated reads based on `test2{,.umi}*{1,2}` (tumor)
      - 'example_hla_pe.bam': Downsampled BAM file for HLATyping workflow / OptiType module. Using existing data did not work as it misses preparation steps.
      - 'example_hla_pe.sorted.bam': Sorted BAM file for HLATyping workflow / OptiType module.
      - 'example_hla_pe.sorted.bam.bai': Sorted BAM file index for HLATyping workflow / OptiType module.
      - mitochon_standin.recalibrated.sorted: copy of the old, smaller test2.paired_end.recalibrated.sorted, this is to be used to test mutect2's mitochondria mode, as the current recal bams are far too big. This should be replaced once rarediseases obtain an actual mitochondria sample.
      - 'test3.single_end.markduplicates.sorted.bam': Mapped, sorted, and duplicate removed reads from ancient DNA across all human chromosomes on the hs37d5 human reference. Data from [ERR2857053](https://www.ebi.ac.uk/ena/browser/view/ERR2857053) downsampled to 10% of original reads.
      - 'test.rna.paired_end.bam': STAR-aligned, unsorted, paired-end RNAseq bam file from the test*rnaseq*{1,2}.fastq.gz: chr22 of sample GM12878 (SRA accession: SRX2900878)
        - 'test.rna.paired_end.sorted.bam': STAR-aligned, sorted, paired-end RNAseq bam file based on test.rna.Aligned.unsorted.bam
        - 'test.rna.paired_end.sorted.bam.bai': STAR-aligned, sorted paired-end RNAseq bam index file for test.rna.Aligned.sorted.bam
      - 'test.rna.paired_end.sorted.chr6.bam': STAR-aligned, sorted, paired-end sampled RNAseq bam file of chromosome 6 of sample GM12878 (SRA accession: SRX2900878)
      - 'test.rna.paired_end.sorted.chr6.bam.bai': STAR-aligned, sorted, paired-end sampled RNAseq bam index file of chromosome 6 of sample GM12878 (SRA accession: SRX2900878)
      - umi:
        - test.paired*end.umi*\*: Files base on `test.umi_{1,2}` (normal)
        - test2.paired*end.umi*\*: Files base on `test2.umi_{1,2}` (tumor)
        - test.paired_end.duplex_umi_unmapped.bam: file originating from `test_duplex_umi\_{1,2}` (spiked)
        - test.paired_end.duplex_umi_mapped.bam: file originating from `test.paired_end.duplex_umi_unmapped.bam`
        - test.paired_end.duplex_umi_mapped_tagged.bam: file originating from `test.paired_end.duplex_umi_unmapped.bam` and `test.paired_end.duplex_umi_mapped.bam`
        - test.paired_end.duplex_umi_grouped.bam: file originating from `test.paired_end.duplex_umi_mapped_tagged.bam`
        - test.paired_end.duplex_umi_duplex_consensus.bam: file originating from `test.paired_end.duplex_umi_grouped.bam`
    - bcl:
      - flowcell.tar.gz: bcl data generated on a MiSeq sequencer. Contains only data for the first tile.
      - flowcell_samplesheet.csv: SampleSheet for data on flowcell
    - cram:
      - test.paired*end.sorted: Mapped, and sorted reads based on `test*{1,2}` (normal)
      - test.paired*end.markduplicates.sorted: Mapped, sorted, and duplicate marked reads based on `test*{1,2}` (normal)
      - test.paired*end.recalibrated.sorted: Mapped, sorted, duplicate marked, and recalibrated reads based on `test*{1,2}` (normal)
      - test2.paired*end.sorted: Mapped, and sorted reads based on `test2*{1,2}` (tumor)
      - test2.paired*end.markduplicates.sorted: Mapped, sorted, and duplicate marked reads based on `test2*{1,2}` (tumor)
      - test2.paired*end.recalibrated.sorted: Mapped, sorted, duplicate marked, and recalibrated reads based on `test2*{1,2}` (tumor)
      - test3.paired*end.recalibrated.sorted: Mapped, sorted, duplicate marked, and recalibrated reads based on `test2*{1,2}` (tumor) Sample is renamed to allow multi-sample testing
    - fastq:

      - test\_{1,2}: reads corresponding to normal sample
      - test.umi\_{1,2}: UMI tagged reads corresponding to normal sample
      - test_duplex_umi\_{1,2}.fastq.gz: duplex UMI tagged reads corresponding to spiked samples (SRA accession: SRR7041712)
      - test2\_{1,2}: reads corresponding to tumor sample
      - test2.umi\_{1,2}: UMI tagged reads corresponding to tumor sample
      - test\_{1,2}germline.fq.gz: Synthetic raw reads file used to generate normal test data for HaplotypeCaller, simulated from chr21
      - test2\_{1,2}germline.fq.gz: Synthetic raw reads file used to generate disease test data for HaplotypeCaller
      - test*rnaseq*{1,2}.fastq.gz: reads from chr22 of sample GM12878 (SRA accession: SRX2900878) for transcriptome analysis.
      - test*airrseq*{umi_R1,R2}.fastq.gz: reads from MiSEQ sequencing of BCR data.
      - rCRS_simulated_test.fq.gz: Synthetic raw mitochondrial reads from the rCRS mitochondrial reference genome for use in testing HaploCart.

    - gatk:
      - test: Recalibration table corresponding to `test{,.umi}_{1,2}` (normal) reads
      - test2: Recalibration table corresponding to `test2{,.umi}_{1,2}` (tumor) reads
      - contig_ploidy_priors_table.tsv: The contig ploidy priors table needed for GATK DetermineGermlineContigPloidy
      - test.pileups.table: Summarises normal allele counts, based on test.paired_end.recalibrated.sorted.bam, used by CalculateContamination
      - test2.pileups.table: Summarises tumor allele counts, based on test2.paired_end.recalibrated.sorted.bam, used by CalculateContamination
      - test_test2_paired_mutect2_calls.artifact-prior.tar.gz: Table of artifact priors, generated from test_test2_paired_mutect2_calls.f1r2.tar.gz, used by FilterMutectCalls
      - test_test2_paired.contamination.table: Table of contamination estimates, generated using test.pileups.table and test2.pileups.table, used by FilterMutectCalls
      - test_test2_paired.segmentation.table: Table of tumor segmentations, generated using test.pileups.table and test2.pileups.table, used by FilterMutectCalls
      - paired_mutect2_calls:
        - test_test2_paired_mutect2_calls.vcf.gz: Output vcf of mutect2 tumor_normal mode based on test.paired_end.recalibrated.sorted.bam (normal) and test2.paired_end.recalibrated.sorted.bam (tumor)
        - test_test2_paired_mutect2_calls.vcf.gz.tbi: Index file for test_test2_paired_mutect2_calls.vcf.gz
        - test_test2_paired_mutect2_calls.vcf.gz.stats: Stats table output along with test_test2_paired_mutect2_calls.vcf.gz
        - test_test2_paired_mutect2_calls.f1r2.tar.gz: Output file generated along with test_test2_paired_mutect2_calls.vcf.gz used by LearnReadOrientationModel to generate artifact_priors
        - test_test2_paired_filtered_mutect2_calls.vcf.gz: tumor normal vcf file after being passed through filtermutectcalls
        - test_test2_paired_filtered_mutect2_calls.vcf.gz.tbi: tumor normal tbi file after being passed through filtermutectcalls
        - test_test2_paired_filtered_mutect2_calls.vcf.gz.filteringStats.tsv: filtering stats file for the tumor normal vcf file after being passed through filtermutectcalls
      - pon_mutect2_calls:
        - test_pon.vcf.gz: variant calls of normal sample run through mutect2 in panel of normals mode
        - test_pon.vcf.gz.tbi: variant calls index of normal sample run through mutect2 in panel of normals mode
        - test_pon.vcf.gz.stats: variant calls stats file of normal sample run through mutect2 in panel of normals mode
        - test2_pon.vcf.gz: variant calls of tumour sample run through mutect2 in panel of normals mode tumour data used as standin would normally be another normal sample in a real pon
        - test2_pon.vcf.gz.tbi: variant calls index of tumour sample run through mutect2 in panel of normals mode tumour data used as standin would normally be another normal sample in a real pon
        - test2_pon.vcf.gz.stats: variant calls stats file of tumour sample run through mutect2 in panel of normals mode tumour data used as standin would normally be another normal sample in a real pon
      - haplotypecaller_calls:
        - test_haplotc.vcf.gz: vcf output from HaplotypeCaller using germline normal reads
        - test_haplotc.vcf.gz.tbi: vcf.tbi output from HaplotypeCaller using germline normal reads
        - test_haplotcaller.cnn.vcf.gz: vcf output from HaplotypeCaller using germline normal reads and run through CNNScoreVariants to add CNN_Key=1D
        - test_haplotcaller.cnn.vcf.gz.tbi: vcf.tbi output from HaplotypeCaller using germline normal reads and run through CNNScoreVariants to add CNN_Key=1D
        - test2_haplotc.vcf.gz: vcf output from HaplotypeCaller using germline disease reads
        - test2_haplotc.vcf.gz.tbi: vcf.tbi output from HaplotypeCaller using germline disease reads
        - test_haplotc.ann.vcf.gz: vcf output from HaplotypeCaller using germline normal reads annotated with snpeff
        - test_haplotc.ann.vcf.gz.tbi: vcf.tbi output from HaplotypeCaller using germline normal reads annotated with snpeff
        - test2_haplotc.ann.vcf.gz: vcf output from HaplotypeCaller using germline disease reads annotated with snpeff
        - test2_haplotc.ann.vcf.gz.tbi: vcf.tbi output from HaplotypeCaller using germline disease reads annotated with snpeff
        - test.g.vcf.gz: vcf output from HaplotypeCaller using germline normal reads in GVCF mode
        - test.g.vcf.gz.tbi: vcf.tbi output from HaplotypeCaller using germline normal reads in GVCF mode
        - test2.g.vcf.gz: vcf output from HaplotypeCaller using germline disease reads in GVCF mode
        - test2.g.vcf.gz.tbi: vcf.tbi output from HaplotypeCaller using germline disease reads in GVCF mode
      - variantrecalibrator:
        - test2.recal: vqsr recalibration table, based on test2_haplotc.ann.vcf.gz
        - test2.recal.idx: vqsr recalibration index, based on test2_haplotc.ann.vcf.gz
        - test2.tranches: vqsr recalibration tranches file, based on test2_haplotc.ann.vcf.gz
        - test2_allele_specific.recal: vqsr allele specific recalibration table, based on test2_haplotc.ann.vcf.gz
        - test2_allele_specific.recal.idx: vqsr allele specific recalibration index, based on test2_haplotc.ann.vcf.gz
        - test2_allele_specific.tranches: vqsr allele specific recalibration tranches file, based on test2_haplotc.ann.vcf.gz
      - test_pon_genomicsdb: Output workspace (directory) from GenomicsdbImport, generated from vcf files in the pon_mutect2_calls subdirectory, used to test CreateSomaticPanelofNormals and GenomicsdbImport, directory has been tar archived to make downloading for tests easier, please remember to untar the directory before using it for testing.
      - test_genomicsdb: Output workspace (directory) from GenomicsdbImport, generated from test.genome.vcf in the gvcf subdirectory, used to test GenotypeGVCFs, directory has been tar archived to make downloading for tests easier.
    - gvcf:
      - test.genome.vcf: Genome vcf corresponding to `test{,.umi}_{1,2}` (normal) reads
      - test2.genome.vcf: Genome vcf corresponding to `test2{,.umi}_{1,2}` (tumor) reads
      - test{,2}.genome.vcf.gz: Bgzipped file based on `test{,2}.genome.vcf` file
      - test{,2}.genome.vcf.gz.tbi: Tbi index based on `test{,2}.genome.vcf.gz` file
      - test{,2}.genome.vcf.idx: Index feature file based on `test{,2}.genome.vcf` file
    - mpileup:
      - test.mpileup.gz: Pileup file correspongind to `test_paired_end_recalibrated_sorted_bam` (normal) generate with `samtools mpileup`
      - test2.mpileup.gz: Pileup file correspongind to `test2_paired_end_recalibrated_sorted_bam` (tumor) generate with `samtools mpileup`
    - broadPeak:
      - test.broadPeak: Genome broadPeak file obtained using MACS2
      - test2.broadPeak: Genome broadPeak file obtained using MACS2, replicate from `test.broadPeak`
    - narrowPeak:
      - test.narrowPeak: Genome narrowPeak file obtained using MACS2
      - test2.narrowPeak: Genome narrowPeak file obtained using MACS2, replicate from `test.narrowPeak`
    - plink
      - test.rnaseq.bed: Plink binaries obtained using test.rnaseq.vcf with plink tool
      - test.rnaseq.bim: Plink binaries obtained using test.rnaseq.vcf with plink tool
      - test.rnaseq.fam: Plink binaries obtained using test.rnaseq.vcf with plink tool
    - vcf:
      - test.rnaseq.vcf: RNAseq vcf corresponding to `test.rnaseq_{1,2}` reads
      - test.genome_21.somatic_sv.vcf: Indels VCF corresponding to `test.paired_end.recalibrated.sorted` and `genome_21.fasta` generated with Manta
      - NA12878*chrM.vcf.gz: mitochondrial variants corresponding to `testdata/NA12878_mito*{1,2}.fq.gz`from the`rarediseases` branch.
      - empty.vcf.gz: The RNAseq VCF with all variants removed
      - empty.vcf.gz.tbi: The index of the empty vcf
    - svd:
      - test.genome.vcf.bed: bed file for markers with format(chr\tpos-1\tpos\trefAllele\taltAllele). Derived from test.genome.vcf and genome.fasta and it was a part of reference stack generated by verifybamid2 --RefVCF test.genome.vcf --Reference genome.fasta
      - test.genome.vcf.mu: matrix file of genotype matrix. Derived from test.genome.vcf and genome.fasta and it was a part of reference stack generated by verifybamid2 --RefVCF test.genome.vcf --Reference genome.fasta
      - test.genome.vcf.UD: matrix file from SVD result of genotype matrix. Derived from test.genome.vcf and genome.fasta and it was a part of reference stack generated by verifybamid2 --RefVCF test.genome.vcf --Reference genome.fasta
      - test.genome.vcf.V: test value file, which was generated by verifybamid2 --RefVCF test.genome.vcf --Reference genome.fasta --BamFile test.paired_end.sorted.bam
    - yak:
      - test.yak: Yak kmer index of 1000 of paternal paired-end reads from the GIAB Ashkenazim trio [RM8392](https://www-s.nist.gov/srmors/view_detail.cfm?srm=8392). These reads were selected from D2_S1_L001_R{1,2}\_001.fastq.gz and D2_S1_L001_R{1,2}\_002.fastq.gz so that they map to `pacbio/fastq/test_hifi.fastq.gz`.
      - test2.yak: Yak kmer index of 1000 of maternal reads from the GIAB Ashkenazim trio [RM8392](https://www-s.nist.gov/srmors/view_detail.cfm?srm=8392). These reads were selected from D3_S1_L001_R{1,2}\_001.fastq.gz and D3_S1_L001_R{1,2}\_001.fastq.gz so that they map to `pacbio/fastq/test_hifi.fastq.gz`.

  - pacbio:

    - bam:
      - alz.bam: raw reads extracted from the [public Alzheimer dataset](https://downloads.pacbcloud.com/public/dataset/IsoSeq_sandbox/2020_Alzheimer8M_subset/alz.1perc.subreads.bam)
      - alz.bam.pbi: pacbio index generated with pbindex
      - alz.ccs.bam: CCS reads generated using pbccs on alz.bam
      - alz.ccs.fl.NEB_5p--NEB_Clontech_3p.bam: set of valid CCS reads generated with LIMA on alz.ccs.bam (keep reads with valid pair of primers, then remove those sequences)
      - alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.bam: set of valid CCS reads generated with isoseq refine (keep reads with a polyA tail, then remove it)
      - alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.bam: set of transcripts generated isoseq cluster
      - alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.bam: set of refined CCS reads not clustered by isoseq cluster
      - alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned.bam: transcripts and singletons aligned on genome2.fa
      - alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned.bam.bai: index file generated with samtools index
    - bed:
      - alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned_tc.bed: first set of gene models generated by TAMA collapse
      - alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned_tc.2.bed: first set of gene models generated by TAMA collapse
    - fasta:
      - alz.ccs.fasta: CCS reads generated using pbccs on alz.bam in fasta format
      - alz.ccs.fasta.gz: CCS reads generated using pbccs on alz.bam in gziped fasta format
      - primers.fasta: NEB Clonetech primers
    - fastq:
      - alz.ccs.fastq: CCS reads generated using pbccs on alz.bam in fastq format
      - alz.ccs.fastq.gz: CCS reads generated using pbccs on alz.bam in gziped fastq format
      - test_hifi.fastq.gz: Reads mapping to a randomly selected contig from the whole genome assembly by [Cheng et al., 2021](https://www.nature.com/articles/s41592-020-01056-5) of the child of the GIAB Ashkenazim trio [RM8392](https://www-s.nist.gov/srmors/view_detail.cfm?srm=8392). The reads were taken from [SRR10382244](https://www.ncbi.nlm.nih.gov/sra/?term=SRR10382244).
    - txt:
      - filelist.txt: A TAMA merge filelist file. It's a 4 columns (bed file, cap status, merging order, id) file listing bed files to merge. The file listed are alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned_tc.bed alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned_tc.2.bed.

  - cooler:

    - cload:
      - hg19:
        - hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz, hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz.px2: hg19 pairix test file and its index file.
        - hg19.GM12878-MboI.pairs.subsample.sorted.possrt.txt.gz, hg19.GM12878-MboI.pairs.subsample.sorted.possrt.txt.gz.tbi: hg19 tabix test file and its index file.
        - hg19.sample1.pairs: hg19 pair test file.
        - hg19.chrom.sizes: hg19 chromosome sizes. Downloaded from [goldenpath](http://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes)
    - merge:
      - toy:
        - toy.symm.upper.2.cool, toy.symm.upper.2.cp2.cool: test file for cooler_merge. Downloaded from [open2c/cooler](https://github.com/open2c/cooler/master/tests/data/toy.symm.upper.2.cool)

  - gene_set_analysis:

    - P53_6samples_collapsed_symbols.gct: a gene cluster text file format (GCT) example
    - P53_6samples.cls: a categorical (e.g tumor vs normal) class file format (CLS) example
    - c1.symbols.reduced.gmx: a GMX (Gene MatriX file format) example

  - array_expression:

    - GSE38751.csv: Sample sheet describing Affy array CEL files
    - GSE38751_RAW.tar: compressed CEL files archive

  - scramble:
    - test.fa: A reference file containing chr3:70000000-70100000 and chr11:418014-438014
    - test.fa.fai: The index of this reference
    - test.bam: A BAM file containing soft-clipped clusters
    - test.bam.bai: The index of the BAM file
    - test.cram: The converted CRAM from the BAM file
    - test.cram.crai: The index of the CRAM file
    - test.bed: A BED file containing only the regions from chr11

- mus_musculus
  - mageck
    - ERR376998.small.fastq.gz and ERR376999.small.fastq.gz downloaded from sourceforge mageck documentation, shortened to only 10k reads
    - design_matrix.txt taken from the mageck documentation tutorial
    - yusa_library.csv taken from the mageck sourceforge, crisprcas9 functional genomics library
    - count_table.csv leukemia mice experiment with crisprcas9 functional genomics
  - gene_set_analysis:
    - mh.all.v2022.1.Mm.symbols.gmt hallmark gene sets, downloaded from [MSigDB](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2022.1.Mm/mh.all.v2022.1.Mm.symbols.gmt) 5/1/2023
    - Mouse_Ensembl_Gene_ID_MSigDB.v2022.1.Mm.chip Ensembl ID to gene symbol mapping in Broad's 'chip' format, suitable for passing to GSEA when using matrices keyed by Ensembl Gene ID. Downloaded from [MSigDB](https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations/mouse/Mouse_Ensembl_Gene_ID_MSigDB.v2022.1.Mm.chip) 5/1/2023
- prokaryotes

  - bacteroides_fragilis
    - genome
      - 'genome.fna.gz': NC_006347 genome downloaded from NCBI Genome
      - 'genome.gbff.gz': NC_006347 genome downloaded from NCBI Genomes in GenBank format
      - 'genome.paf': NC_006347 genome PAF file
      - 'genome.mapping.potential.ARG': DeepARG anti-microbial resistance gene prediction of NC_006347 genome
      - 'genome.gff.gz': NC_006347 genome downloaded from NCBI Genomes in GFF3 format
    - illumina
      - fastq
        - 'test1*{1,2}.fastq.gz': synthetic raw short-read sequencing reads of the genome of the mammalian-gut-residing Bacteroides fragilis* YCH46 bacterium (NC_006347). Originally generated for the [MAG pipeline test dataset](https://github.com/nf-core/test-datasets/tree/mag).
        - 'test2*{1,2}.fastq.gz': synthetic raw short-read sequencing reads of the genome of the mammalian-gut-residing Bacteroides fragilis* YCH46 bacterium (NC_006347). Originally generated for the [MAG pipeline test dataset](https://github.com/nf-core/test-datasets/tree/mag).
      - fasta
        - 'test1.contigs.fa.gz': _de novo_ assembled contigs of the test\minigut_sample_1 FASTQ files by MEGAHIT, generated with nf-core/mag (2.1.0) on default settings
      - bam
        - 'test1.bam': 'test1\_{1,2}.fastq.gz' file aligned with bowtie2 on 'genome.fna.gz'
        - 'test1.sorted.bam': sorted 'test1.bam'
        - 'test1.sorted.bai': index of 'test1.sorted.bam'
        - 'test2.bam': 'test2\_{1,2}.fastq.gz' file aligned with bowtie2 on 'genome.fna.gz'
        - 'test2.sorted.bam': sorted 'test2.bam'
        - 'test2.sorted.bai': index of 'test2.sorted.bam'
    - nanopore
      - fastq
        - 'test.fastq.gz' synthetic raw long-read sequencing reads of the genome of the mammalian-gut-residing _Bacteroides fragilis_ YCH46 bacterium (NC_006347). Originally generated for the [MAG pipeline test dataset](https://github.com/nf-core/test-datasets/tree/mag).

- galaxea_fascicularis
  - hic
    - 'jaGalFasc40_2.pretext': sparse data pretext map of hic contacts

### pangenomics

- homo_sapiens
  - 'pangenome.fa': A FASTA file which contains several related genomes.
  - 'pangenome.fa.gz': A GZIP compressed FASTA file which contains several related genomes.
  - 'pangenome.paf': A PAF file which contains the all versus all pairwise alignments of related genomes.
  - 'pangenome.paf.gz': A GZIP compressed PAF file which contains the all versus all pairwise alignments of related genomes.
  - 'pangenome.seqwish.gfa': A GFA file which contains the pangenome graph induced by `seqwish` encoded in the variation graph model.
  - 'pangenome.smoothxg.gfa': A GFA file which contains the `smoothxg` smoothed pangenome graph.
  - 'pangenome.gfaffix.gfa': A GFA file which was normalized with `gfaffix`.
  - 'pangenome.vcf.gz': A BGZIPPED VCF file whic was created by `vg deconstruct`. Contains the `LV=0` tag in the INFO field. This specifies the snarl level of the variation.
  - odgi
    - 'pangenome.og': A variation graph encoded in the binary ODGI format. It is consumed by `odgi view`.
    - 'pangenome.lay': A binary file which holds the 2D graph layout produced by `odgi layout`. Input for `odgi draw`.

### proteomics

- msspectra

  - 'OVEMB150205_12.raw': Thermo RAW mass spectra file.
  - 'OVEMB150205_14.raw': Thermo RAW mass spectra file.

- database

  - 'yeast_UPS.fasta': FASTA database for Yeast organism.

- maxquant

  - 'mqpar.xml': MaxQuant parameters file.

- Bos_taurus
  - idxml:
    - 'BSA_QC_file.idXML': identifications file for Bovine (BSA) organism coming from an OMSSA search engine.
  - raw:
    - 'BSA_QC_file.raw': spectra file for Bovine (BSA) organism.

### spatialomics

- tiff
  - 'mindagap.mouse_heart.wga.tiff': Exemplary tiff file with black gridlines to fill for MindaGap tool.
- h5
  - 'plant_wga.h5' : Image of rice root stained with Wheat-Germ agglutinin (WGA) from :
  - 'plant_wga_probabilities.h5' : Probability maps from pixel classification workflow (plant_wga.pixel_prob.ilp).
- ilp
  - 'plant_wga.multicut.ilp' : Ilastik project file for multicut. Output format is set to tiff.
  - 'plant_wga.pixel_prob.ilp' : Ilastik project file for pixel classification trained on plant_wga.h5

### generic

- csv
  - 'test.csv': exemplary comma-separated file obtained from [here](https://bioinf.shenwei.me/csvtk/usage/#split)
- notebooks
  - jupyter
    - 'ipython_notebook.ipynb': exemplary jupyter notebook
    - 'ipython_notebook.md': exemplary markdown notebook
  - rmarkdown
    - 'rmarkdown_notebook.Rmd': exemplary R notebook
- tsv
  - 'test.tsv': exemplary tab-separated file obtained from [here](https://bioinf.shenwei.me/csvtk/usage/#split)
- txt
  - 'hello.txt': one-line txt file
- tar
  - 'hello.tar.gz': gzipped tar archive containing a single file without a directory

### Uncategorized

- e_coli_k12_16s.fna: E. coli K-12 16S rRNA
- bac.16S_rRNA.hmm: Bacterial 16S HMM file
