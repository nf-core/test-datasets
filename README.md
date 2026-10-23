# test-datasets: `rnasplice`

This branch contains test data to be used for automated testing with the [nf-core/rnasplice](https://github.com/nf-core/rnasplice) pipeline.

## Content of this repository

`samplesheet/samplesheet.csv`: Experiment design file for minimal test dataset

`samplesheet/samplesheet_full.csv`: Experiment design file for full test dataset

`reference/`: Sub-sampled genome reference files (ChrX)

`reference/misopy_index.tar.gz`: Misopy index files derived from `reference/genes_chrX.gff3`

`samplesheet/samplesheet_strandedness_auto.csv`: Experiment design file for the reverse stranded test dataset, with two samples set to `auto`

`testdata/*.fastq.gz`: Subsampled fastq files (ChrX)

`testdata/*_chrX_reverse_{1,2}.fastq.gz`: Reverse stranded fastq files derived from the subsampled fastq files (ChrX)

`scripts/make_reverse_stranded.py`: Script used to derive the reverse stranded fastq files

`testdata/salmon_quant/`: Salmon quantification results for the subsampled fastq files (ChrX)

`testdata/drimseq/`: test files for DRIMSeq

`testdata/dexseq/`: test files for DEXSeq subworkflow

`chrX_subset/`: The minimal test dataset cut down to a few regions of chrX, with every file under 1 MB and the same layout and file names as above (`reference/`, `samplesheet/`, `testdata/`)

`scripts/make_chrx_subset.py`, `scripts/filter_pairs.py`, `scripts/package_chrx_subset.sh`: Scripts used to make `chrX_subset/`

## Minimal test dataset origin

ChrX subsampled *H.sapiens* total RNA paired-end RNA-seq data was obtained from:

> Pertea, M., Kim, D., Pertea, G. et al. Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nat Protoc 11, 1650–1667 (2016). doi: 10.1038/nprot.2016.095. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/27560171/)

Original data source:

> Lappalainen T, Sammeth M, Friedländer MR, 't Hoen PA, Monlong J, Rivas MA, Gonzàlez-Porta M, Kurbatova N, Griebel T, Ferreira PG, Barann M, Wieland T, Greger L, van Iterson M, Almlöf J, Ribeca P, Pulyakhina I, Esser D, Giger T, Tikhonov A, Sultan M, Bertier G, MacArthur DG, Lek M, Lizano E, Buermans HP, Padioleau I, Schwarzmayr T, Karlberg O, Ongen H, Kilpinen H, Beltran S, Gut M, Kahlem K, Amstislavskiy V, Stegle O, Pirinen M, Montgomery SB, Donnelly P, McCarthy MI, Flicek P, Strom TM; Geuvadis Consortium, Lehrach H, Schreiber S, Sudbrak R, Carracedo A, Antonarakis SE, Häsler R, Syvänen AC, van Ommen GJ, Brazma A, Meitinger T, Rosenstiel P, Guigó R, Gut IG, Estivill X, Dermitzakis ET. Transcriptome and genome sequencing uncovers functional variation in humans. Nature. 2013 Sep 26;501(7468):506-11. doi: 10.1038/nature12531. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/24037378/) [SRA](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=study&acc=ERP001942)

Of all samples only the following were chosen for testing: ERR188383, ERR188428, ERR188454 and ERR204916.

### Sampling information

| SRA ID | Sex        | Population | Sample title                                                                  |
|---------------|------------|------------|-------------------------------------------------------------------------------|
| ERR188383     | male       | GBR        | ChrX subsampled total RNA-seq male GBR (British from England)                 |
| ERR188428     | female     | GBR        | ChrX subsampled total RNA-seq female GBR (British from England)               |
| ERR188454     | male       | YRI        | ChrX subsampled total RNA-seq male YRI (Yoruba from Ibadan, Nigeria)          |
| ERR204916     | female     | YRI        | ChrX subsampled total RNA-seq female YRI (Yoruba from Ibadan, Nigeria) |

### Sampling procedure

ChrX fastq files we further sub-sampled to improve speed of testing.

1. The example command below was used to sub-sample the raw paired-end FastQ files to 50,000 reads (see [seqtk](https://github.com/lh3/seqtk)):

  ```console
  seqtk sample -s100 ERR188383_unclass_chrX_1.fastq.gz 50000 | gzip > ERR188383_chrX_1.fastq.gz
  seqtk sample -s100 ERR188383_unclass_chrX_2.fastq.gz 50000 | gzip > ERR188383_chrX_2.fastq.gz
  ```
 2. Ensembl GRCh37 annotation was downloaded from iGenomes (https://ewels.github.io/AWS-iGenomes/) and subsetted for ChrX with the following command:

  ```console
  awk -F"\t" '$1=="X"' genes.gtf > genes_chrX.gtf
  ```
 3. Ensembl GRCh37 ChrX fasta file was downloaded from iGenomes (https://ewels.github.io/AWS-iGenomes/).

### Reverse stranded test dataset

The GEUVADIS libraries are unstranded. To test strandedness inference on stranded data, a reverse stranded version of each sample was derived from the subsampled fastq files. In a pair whose read 1 aligns to the same strand as the transcripts it overlaps, the two mates were swapped, so read 1 is antisense to the transcript, as in a dUTP library. Pairs that overlap transcripts on both strands, or no transcript, were left unchanged. The files keep all 50,000 pairs.

  ```console
  samtools view -f 0x40 -F 0x904 testdata/ERR188383.Aligned.out.bam | cut -f1-6 > ERR188383.read1.tsv
  python3 scripts/make_reverse_stranded.py reference/genes_chrX.gtf ERR188383.read1.tsv \
      testdata/ERR188383_chrX_1.fastq.gz testdata/ERR188383_chrX_2.fastq.gz \
      ERR188383_chrX_reverse_1.fastq.gz ERR188383_chrX_reverse_2.fastq.gz
  ```

About 88% of the pairs end up with read 1 antisense to their transcript. nf-core/rnasplice calls them reverse stranded (Salmon about 96% antisense, RSeQC `infer_experiment.py` about 91%).

### chrX subset

The minimal test dataset builds its STAR and Salmon indices from all of chrX (156 Mb, 2,392 genes), although most reads fall on a few hundred genes. `chrX_subset/` is a smaller version of it, with every file under 1 MB (20 MB in all).

- Reference: the regions around the 30 genes with at least two transcripts and the most reads across the four samples, plus the three genes nf-core/rnasplice plots with MISO by default (ENSG00000004961, ENSG00000005302 and ENSG00000147403). Each gene is padded by 1 kb, regions closer than 10 kb are merged, and a region is widened until no gene crosses its ends. The 30 regions are joined into one 3.4 Mb contig, still named `X`, separated by 1 kb of `N`.
- Annotation: the 82 genes inside the regions, in the new coordinates. Only the `transcript`, `exon` and `CDS` lines and the `gene_id`, `transcript_id`, `gene_name`, `gene_biotype` and `exon_number` attributes are kept, which is what the pipeline uses.
- Reads: 10,000 pairs per sample, drawn with a fixed seed from the pairs whose two mates start in the same region.

  ```console
  # columns 1-8 of the primary alignments of the four samples
  for s in ERR188383 ERR188428 ERR188454 ERR204916; do
      samtools view -F 0x904 testdata/$s.Aligned.out.bam | cut -f1-8
  done > alignments.tsv
  python3 scripts/make_chrx_subset.py X.fa reference/genes_chrX.gtf alignments.tsv 30 \
      ENSG00000004961,ENSG00000005302,ENSG00000147403 10000 subset
  python3 scripts/filter_pairs.py subset/keep_reads.txt \
      testdata/ERR188383_chrX_1.fastq.gz testdata/ERR188383_chrX_2.fastq.gz \
      subset/testdata/ERR188383_chrX_1.fastq.gz subset/testdata/ERR188383_chrX_2.fastq.gz
  ```

The other files come from a run of nf-core/rnasplice (`-profile test` with `--leafcutter --save_dexseq_annotation`) on this reference and these reads, packaged by `scripts/package_chrx_subset.sh`:

- `reference/genes_chrX.gff3` and `reference/misopy_index.tar.gz`: from `GFFREAD_GTF2GFF3` (`gffread -L --keep-genes`) and `MISOPY_INDEX`
- `testdata/*.Aligned.out.bam` and `.bai`: the coordinate sorted STAR alignments; `testdata/*.Aligned.toTranscriptome.out.bam`: the STAR transcriptome alignments. Both keep 4,000 of the pairs of each sample with a spliced mate (`samtools view -N`), to stay under 1 MB and still give the splicing tools junctions to work on
- `testdata/ERR*/`, `testdata/ERR*.tar.gz` and `testdata/salmon_quant/salmon_quant.tar.gz` (with `genes_chrX.tx2gene.tsv` from `GFFREAD_TX2GENE`): the Salmon quantification of the trimmed reads
- `testdata/drimseq/drimseq_testdata.tar.gz`: the tximport DTU object and tx2gene table after STAR and Salmon, with a samplesheet of the four samples
- `testdata/dexseq/genes_chrX.DEXSeq.gff.gz`: the DEXSeq annotation prepared from the GTF
- `testdata/*_chrX_reverse_{1,2}.fastq.gz`: made from these reads and their alignments as described above
- `testdata/*_reverse.Aligned.out.bam` and `.bai`: the coordinate sorted STAR alignments of the reverse stranded reads (run with `strandedness` set to `reverse`), cut down to 4,000 pairs with a spliced mate in the same way, for tests of the stranded LeafCutter path

On this data the STAR index drops from 1.7 GB to 42 MB and the Salmon index from 0.9 GB to 23 MB, STAR needs under 1 GB to align instead of 4.9 GB, and every analysis still reports results.

## Full test dataset origin

*H. sapiens* total RNA paired-end RNA-seq dataset was obtained from:

> Shen S, Park JW, Lu ZX, Lin L, Henry MD, Wu YN, Zhou Q, Xing Y. rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. Proc Natl Acad Sci U S A. 2014 Dec 23;111(51):E5593-601. doi: 10.1073/pnas.1419161111. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/25480548/) [SRA](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP014759)

### Sample information full test dataset

| SRA ID    | Title           | Platform            | Sample title                                                        |
|-----------|-----------------|---------------------|---------------------------------------------------------------------|
| SRR536350 | PC3E sample     | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line E-cadherin positive PC-3 cells |
| SRR536348 | PC3E sample     | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line E-cadherin positive PC-3 cells |
| SRR536352 | PC3E sample     | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line E-cadherin positive PC-3 cells |
| SRR536346 | GS689.Li sample | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line isolated from a secondary metastatic liver tumor after intravenous injection of PC-3 cells into mouse |
| SRR536342 | GS689.Li sample | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line isolated from a secondary metastatic liver tumor after intravenous injection of PC-3 cells into mouse |
| SRR536344 | GS689.Li sample | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line isolated from a secondary metastatic liver tumor after intravenous injection of PC-3 cells into mouse |
