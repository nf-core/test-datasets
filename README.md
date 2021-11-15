# Modules Test Data

This branch of the `nf-core/test-datasets` repository contains all data used for the individual module tests.
There are three main directories: `generic`, `genomics` and `delete_me`. The first contains generic files, the second contains all datasets for genomics tools while the latter contains temporary datasets that will be deleted as better data gets available.

### delete_me

The `delete_me` folder does not adhere to a defined structure as data in this folder should be directly as fast as possible, whenever a more suitable dataset is found that can be added to the `genomics` folder.

### generic

The `generic` folder contains generic files that currently cannot be associated to a genomics category. They are organised by their respective file extension. Also, it contains depreciated files which will be removed in the future and exchanged by files in `genomics`.

### genomics

The genomics folder contains subfolders for all organisms for which test data is available. At the moment, there are three organisms available:

* bacteroides_fragilis
* homo_sapiens
* sarscov2

The three folders are structured in a similar way, with any genome-specific files in `genome` (e.g. fasta, gtf, ...) and technology specific raw-data files in the `illumina`, `nanopore`, `pacbio` and `cooler` subfolders whenever available.
`Genomics` contains all typical data required for genomics modules, such as fasta, fastq and bam files. Every folder in `genomics` corresponds to a single organism. For every data file, a short description about how this file was generated is available either in this description or in the respective subfolder.

If you cannot find suitable test data on this repository, please contact us on the [nf-core Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)). The goal is to create a suitable, small test dataset that should be usable with the available test data and if possible be generated with modules available from `nf-core/modules`. If creating that test data is difficult but you want to add the module first, it is also possible to add a small subset to the `delete_me` folder to get your module tests running, and then add proper test data afterwards. This should be discussed on slack. In order to add test data. For a short description of the workflow for adding new data, have a look at [here](docs/ADD_NEW_DATA.md).

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
    * maltextract
      * 'taxon_list.txt': text file of list NCBI sarscov2 species IDs primarily used for MaltExtract
      * 'ncbi_taxmap.zip': mini-NCBI taxonomy map prmiarily used for MaltExtract
  * genome
    * 'genome.fasta': MT192765.1 genome including (GCA_011545545.1_ASM1154554v1)
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
      * 'test.paired_end.name.sorted.bam': Paired-end bam file sorted by name
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
    * gfa
      * 'assembly.gfa': assembly in Graphical Fragment Assembly (GFA) 1.0 format
      * 'assembly.gfa.bgz': compressed with block-compressed GZIP (BGZF)
      * 'assembly.gfa.gz': compressed with GZIP
      * 'assembly.gfa.zst': compressed with Zstandard (zstd)
    * sra
      * `SRR13255544.tar.gz`: Tar archive containing SRA file obtained from SRR13255544.
      * `SRR11140744.tar.gz`: Tar archive containing SRA file obtained from SRR11140744.
    * vcf
      * 'test.vcf', 'test2.vcf': generated from 'test_paired_end.sorted.bam' using bcftools mpileup, call and filter
      * 'test3.vcf': generated from 'test_single_end.sorted.bam' using bcftools mpileup, call and filter
      * 'test2.targets.tsv.gz' from 'test2.vcf.gz' using bcftools query  and bgzip
      * '*.gz': generated from VCF files using bgzip
      * '.tbi': generated from '.vcf.gz' files using `tabix -p vcf -f <file>`
      * ped
        * 'justhusky_minimal.vcf.gz': minimal combination example of VCF/PED file
        * 'justhusky.ped': minimal combination example of VCF/PED file
        * '.tbi': generated from '.vcf.gz' files using `tabix -p vcf -f <file>`
    * wig
      * 'test.wig.gz'
    * picard
      * 'test.single_end.bam.readlist.txt': text file of a list of two read IDs primarily for picard FilterSamReads  
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
    * BUSCO
      * 'chr22_odb10.tar.gz': BUSCO database 'primates_odb10.2021-02-19.tar.gz' purged of entries not matching 'genome.fasta'.
    * vcf
      * dbsnp: DBSnp file downsampled based on reference position
      * gnomAD: gnomAD file downsampled based on reference position
      * mills_and_1000G: Indels file downsampled based on reference position
    * dict: Sequence dictionary corresponding to `fasta`
    * genome.fasta: Reference fasta based on chr22:16570000-16610000
    * genome2.fasta: Reference fasta based on chr22:16600000-16800000
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
      * 'example_hla_pe.bam': Downsampled BAM file for HLATyping workflow / OptiType module. Using existing data did not work as it misses preparation steps.
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
      * test.pileups.table: Summarises normal allele counts, based on test.paired_end.recalibrated.sorted.bam, used by CalculateContamination
      * test2.pileups.table: Summarises tumor allele counts, based on test2.paired_end.recalibrated.sorted.bam, used by CalculateContamination
      * test_test2_paired_mutect2_calls.artifact-prior.tar.gz: Table of artifact priors, generated from test_test2_paired_mutect2_calls.f1r2.tar.gz, used by FilterMutectCalls
      * test_test2_paired.contamination.table: Table of contamination estimates, generated using test.pileups.table and test2.pileups.table, used by FilterMutectCalls
      * test_test2_paired.segmentation.table: Table of tumor segmentations, generated using test.pileups.table and test2.pileups.table, used by FilterMutectCalls
      * paired_mutect2_calls:
        * test_test2_paired_mutect2_calls.vcf.gz: Output vcf of mutect2 tumor_normal mode based on test.paired_end.recalibrated.sorted.bam (normal) and test2.paired_end.recalibrated.sorted.bam (tumor)
        * test_test2_paired_mutect2_calls.vcf.gz.tbi: Index file for test_test2_paired_mutect2_calls.vcf.gz
        * test_test2_paired_mutect2_calls.vcf.gz.stats: Stats table output along with test_test2_paired_mutect2_calls.vcf.gz
        * test_test2_paired_mutect2_calls.f1r2.tar.gz: Output file generated along with test_test2_paired_mutect2_calls.vcf.gz used by LearnReadOrientationModel to generate artifact_priors
      * test_genomicsdb: Output workspace (directory) from GenomicsdbImport, generated from test.genome.vcf, only one sample used to minimize size, used to test CreateSomaticPanelofNormals and GenomicsdbImport, directory has been tar archived to make downloading for tests easier, please remember to untar the directory before using it for testing.
    * gvcf:
      * test.genome.vcf: Genome vcf corresponding to `test{,.umi}_{1,2}` (normal) reads
      * test2.genome.vcf: Genome vcf corresponding to `test2{,.umi}_{1,2}` (tumor) reads
      * test{,2}.genome.vcf.gz: Bgzipped file based on `test{,2}.genome.vcf` file
      * test{,2}.genome.vcf.gz.tbi: Tbi index based on `test{,2}.genome.vcf.gz` file
      * test{,2}.genome.vcf.idx: Index feature file based on `test{,2}.genome.vcf` file
    * broadPeak:
      * test.broadPeak: Genome broadPeak file obtained using MACS2
      * test2.broadPeak: Genome broadPeak file obtained using MACS2, replicate from `test.broadPeak`
    * narrowPeak:
      * test.narrowPeak: Genome narrowPeak file obtained using MACS2
      * test2.narrowPeak: Genome narrowPeak file obtained using MACS2, replicate from `test.narrowPeak`
    * vcf:
      * test.rnaseq.vcf: RNAseq vcf corresponding to `test.rnaseq_{1,2}` reads
    * yak:
      * test.yak: Yak kmer index of 1000 of paternal paired-end reads from the GIAB Ashkenazim trio [RM8392](https://www-s.nist.gov/srmors/view_detail.cfm?srm=8392). These reads were selected from D2_S1_L001_R{1,2}_001.fastq.gz and D2_S1_L001_R{1,2}_002.fastq.gz so that they map to `pacbio/fastq/test_hifi.fastq.gz`.
      * test2.yak: Yak kmer index of 1000 of maternal reads from the GIAB Ashkenazim trio [RM8392](https://www-s.nist.gov/srmors/view_detail.cfm?srm=8392). These reads were selected from D3_S1_L001_R{1,2}_001.fastq.gz and D3_S1_L001_R{1,2}_001.fastq.gz so that they map to `pacbio/fastq/test_hifi.fastq.gz`.
  * pacbio:
    * bam:
      * alz.bam: raw reads extracted from the [public Alzheimer dataset](https://downloads.pacbcloud.com/public/dataset/IsoSeq_sandbox/2020_Alzheimer8M_subset/alz.1perc.subreads.bam)
      * alz.bam.pbi: pacbio index generated with pbindex
      * alz.ccs.bam: CCS reads generated using pbccs on alz.bam
      * alz.ccs.fl.NEB_5p--NEB_Clontech_3p.bam: set of valid CCS reads generated with LIMA on alz.ccs.bam (keep reads with valid pair of primers, then remove those sequences)
      * alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.bam: set of valid CCS reads generated with isoseq refine (keep reads with a polyA tail, then remove it)
      * alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.bam: set of transcripts generated isoseq cluster
      * alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.bam: set of refined CCS reads not clustered by isoseq cluster
      * alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned.bam: transcripts and singletons aligned on genome2.fa
      * alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned.bam.bai: index file generated with samtools index
    * bed:
      * alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned_tc.bed: first set of gene models generated by TAMA collapse
      * alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned_tc.2.bed: first set of gene models generated by TAMA collapse
    * fasta:
      * alz.ccs.fasta: CCS reads generated using pbccs on alz.bam in fasta format
      * alz.ccs.fasta.gz: CCS reads generated using pbccs on alz.bam in gziped fasta format
      * primers.fasta: NEB Clonetech primers
    * fastq:
      * alz.ccs.fastq: CCS reads generated using pbccs on alz.bam in fastq format
      * alz.ccs.fastq.gz: CCS reads generated using pbccs on alz.bam in gziped fastq format
      * test_hifi.fastq.gz: Reads mapping to a randomly selected contig from the whole genome assembly by [Cheng et al., 2021](https://www.nature.com/articles/s41592-020-01056-5) of the child of the GIAB Ashkenazim trio [RM8392](https://www-s.nist.gov/srmors/view_detail.cfm?srm=8392). The reads were taken from [SRR10382244](https://www.ncbi.nlm.nih.gov/sra/?term=SRR10382244).
    * txt:
      * filelist.txt: A TAMA merge filelist file. It's a 4 columns (bed file, cap status, merging order, id) file listing bed files to merge. The file listed are alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned_tc.bed alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned_tc.2.bed.

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

* bacteroides_fragilis
  * genome
      * 'genome.fna.gz': NC_006347 genome downloaded from NCBI Genome           
      * 'genome.gbff.gz': NC_006347 genome downloaded from NCBI Genomes in GenBank format
  * illumina
    * fastq
      * 'test1_{1,2}.fastq.gz': synthetic raw short-read sequencing reads of the genome of the mammalian-gut-residing Bacteroides fragilis_ YCH46  bacterium (NC_006347). Originally generated for the [MAG pipeline test dataset](https://github.com/nf-core/test-datasets/tree/mag).
      * 'test2_{1,2}.fastq.gz': synthetic raw short-read sequencing reads of the genome of the mammalian-gut-residing Bacteroides fragilis_ YCH46 bacterium (NC_006347). Originally generated for the [MAG pipeline test dataset](https://github.com/nf-core/test-datasets/tree/mag).
    * fasta
      * 'test1.contigs.fa.gz': _de novo_ assembled contigs of the test\minigut\_sample_1 FASTQ files by MEGAHIT, generated with nf-core/mag (2.1.0) on default settings
    * bam
      * 'test1.bam': 'test1_{1,2}.fastq.gz' file aligned with bowtie2 on 'genome.fna.gz'
      * 'test1.sorted.bam': sorted 'test1.bam'
      * 'test1.sorted.bai': index of 'test1.sorted.bam'
      * 'test2.bam': 'test2_{1,2}.fastq.gz' file aligned with bowtie2 on 'genome.fna.gz'
      * 'test2.sorted.bam': sorted 'test2.bam'
      * 'test2.sorted.bai': index of 'test2.sorted.bam'
  * nanopore
    * fastq
      * 'test.fastq.gz' synthetic raw long-read sequencing reads of the genome of the mammalian-gut-residing _Bacteroides fragilis_ YCH46 bacterium (NC_006347). Originally generated for the [MAG pipeline test dataset](https://github.com/nf-core/test-datasets/tree/mag).

## generic

* csv
  * 'test.csv': exemplary comma-separated file obtained from [here](https://bioinf.shenwei.me/csvtk/usage/#split)
* notebooks
  * jupyter
    * 'ipython_notebook.ipynb': exemplary jupyter notebook
    * 'ipython_notebook.md': exemplary markdown notebook
  * rmarkdown
    * 'rmarkdown_notebook.Rmd': exemplary R notebook
* tsv
  * 'test.tsv': exemplary tab-separated file obtained from [here](https://bioinf.shenwei.me/csvtk/usage/#split)
* txt
  * 'hello.txt': one-line txt file

### Uncategorized

* e_coli_k12_16s.fna: E. coli K-12 16S rRNA
* bac.16S_rRNA.hmm: Bacterial 16S HMM file
