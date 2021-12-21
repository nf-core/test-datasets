# test-datasets: `nanoseq`

This branch contains test data to be used for automated testing with the [nf-core/nanoseq](https://github.com/nf-core/nanoseq) pipeline.

## Basecalling barcoded data

The barcoded data in this repository will be used to test the pipeline from end-to-end. The associated parameters and settings to run the default tests for the pipeline can be found in [`test.config`](https://github.com/nf-core/nanoseq/blob/master/conf/test.config).

### Files

* `samplesheet_bc_dx.csv` - Sample information sheet required for the pipeline
* `fast5/barcoded/` - Subset of fast5 files from direct cDNA Nanopore reads for HepG2 (Liver Cancer) and K562 (Leukemia) cell lines

### Sequencing information

|             |         	 |
|-------------|------------|
| Flow Cell   | FLO-MIN106 |
| Kit         | SQK-DCS109 |
| Barcode Kit | EXP-NBD103 |

## Basecalling non-barcoded data

The non-barcoded data in this repository will be used to test the pipeline without the demultiplexing step. The associated parameters and settings to run the pipeline can be found in [`test_bc_nodx.config`](https://github.com/nf-core/nanoseq/blob/master/conf/test_bc_nodx.config).

### Files

* `samplesheet_bc_nodx.csv` - Sample information sheet required for the pipeline
* `fast5/nonbarcoded/` - Subset of fast5 files from direct cDNA Nanopore reads for the MCF7 (Breast Cancer) cell line

### Sequencing information

|             |         	 |
|-------------|------------|
| Flow Cell   | FLO-MIN106 |
| Kit         | SQK-DCS108 |
| Barcode Kit | None    	 |

## Pre-basecalled and demultiplexed data

The pre-basecalled and demultiplexed data in this repository will be used to test the pipeline without both the basecalling and demultiplexing steps. The associated parameters and settings to run the pipeline can be found in [`test_nobc_nodx.config`](https://github.com/nf-core/nanoseq/blob/master/conf/test_nobc_nodx.config).

### Files

* `samplesheet_nobc_nodx.csv` - Sample information sheet required for the pipeline
* `fastq/demultiplexed/` - FastQ files for barcodes 1 and 2 obtained from [Basecalling barcoded data](#basecalling-barcoded-data)

## Pre-basecalled and nondemultiplexed data

The pre-basecalled and nondemultiplexed data in this repository will be used to test the pipeline without both the basecalling and demultiplexing steps. The associated parameters and settings to run the pipeline can be found in [`test_nobc_dx.config`](https://github.com/nf-core/nanoseq/blob/master/conf/test_nobc_dx.config).

### Files

* `samplesheet_nobc_dx.csv` - Sample information sheet required for the pipeline
* `fastq/nondemultiplexed/` - Non-demuliplexed FastQ files for the nanopore DNA reads from the Hct116 (Colon Cancer) cell line.

### Sequencing information

|             |         	    |
|-------------|---------------|
| Flow Cell   | FLO-MIN106    |
| Kit         | SQK-LSK109    |
| Barcode Kit | NBD103/NBD104 |

## Pre-basecalled and demultiplexed data for variant calling

The pre-basecalled and demultiplexed data in this repository will be used to test the dna variant calling without either of the basecalling or demultiplexing steps. The associated parameters and dettings to run the pipeline can be found in [`test_nobc_nodx_vc.config`](https://github.com/nf-core/nanoseq/blob/master/conf/test_nobc_nodx_vc.config)

### Files

* `samplesheet_nobc_nodx_vc.csv` - Sample information sheet required for the pipeline
* `fastq/demultiplexed/` - Demuliplexed FastQ files for the nanopore DNA reads overlapping the EDIL3 gene for NA12878.

## Aligned data

The aligned data in this repository will be used to test the pipeline without the basecalling, demultiplexing and alignment step. The associated parameters and settings to run the pipeline can be found in [`test_nobc_nodx_noaln.config`](https://github.com/nf-core/nanoseq/blob/master/conf/test_nobc_nodx_noaln.config).

### Files

* `samplesheet_nobc_nodx_noaln.csv` - Sample information sheet required for the pipeline
* `bam/` - Bam files obtained from [Pre-basecalled and nondemultiplexed data](#pre-basecalled-and-nondemultiplexed-data)


## Full-sized test data

The full sized test data in this repository will be used to test the pipeline without the basecalling and demultiplexing step. The associated parameters and settings to run the pipeline can be found in [`test_full.config`](https://github.com/nf-core/nanoseq/blob/master/conf/test_full.config).

### Origin

*H. sapiens* Nanopore cDNA and direct cDNA datasets were obtained from:
 > The Singapore Nanopore Expression Consortium. (2020). The Singapore Nanopore Expression Project (SG-NEx) data pre-release v0.1 (Version v0.1-pre-release) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.4159715

The sample information are listed here: 

| Cell line     | Replicate | Sequencing protocol |
|---------------|-----------|---------------------|
| A549          |	1	        | cDNA                |
| A549          | 2         | direct cDNA         |
| A549          |	3	        | direct cDNA         |
| K562          | 1         | cDNA                |
| K562          | 2         | direct cDNA         |
| K562          |	3	        | direct cDNA         |

### Files

* `samplesheet_full.csv` - Sample information sheet required for the pipeline


## Reference genome

The test-datasets in this repository were derived from human samples. The size of the entire human genome is too large and possibly too excessive to test the functionality of the pipeline from end-to-end. To overcome this, the data was initially mapped to the human genome and visually inspected. Two genes, KCMF1 and EDIL3, were selected to represent the reference genome for different tests.

### Files

* `reference/hg19_KCMF1.fa` - DNA for KCMF1 gene +- 1kb obtained from the `hg19` UCSC human genome assembly
* `reference/GRCh38_EDIL3.fa` - DNA from the EDIL3 gene obtained from the GRCh38 human genome assembly.

### Obtaining DNA sequences

The [UCSC Genome Browser](https://genome.ucsc.edu) and other methods can also be used to obtain the gene interval or the DNA sequence directly. The approach outlined below is more flexible for instances where the reference genome isnt hosted on UCSC or for custom interval sets.

#### Creating a BED file of gene intervals

The intervals for KCMF1 were obtained by:
* loading the "hg19" genome in [IGV](http://software.broadinstitute.org/software/igv/)
* searching for "KCMF1" in the search box
* right-clicking on the gene interval in the "RefSeq Genes" track and selecting "Copy Details To Clipboard"

The intervals for EDIL3 were obtained by:
* loading the "GRCh38" genome in [IGV](http://software.broadinstitute.org/software/igv/)
* searching for "EDIL3" in the search box
* right-clicking on the gene interval in the "RefSeq Genes" track and selecting "Copy Details To Clipboard"

The outputs should look like:

```bash
chr2:85198231-85286595
KCMF1
chr2:85198231-85286595
id = NM_020122
```

```bash
chr5:83940554-84384880
EDIL3
chr5:83940554-84384880
id = NM_005711.5
```

This information was reformatted in order to create [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) files called `hg19_KCMF1.bed` and `GRCh38_EDIL3.bed`, respectively.

```bash
chr2    85198230    85286595    KCMF1   0   +
```

```bash
chr5	83940554	84384880	EDIL3	0	-
```

> NB: The BED format uses a 0-based coordinate system so the start position has been adjusted accordingly.

#### Generate chromosome sizes file for genome

We need to use [BEDTools](https://github.com/arq5x/bedtools2/) to extract the DNA sequence for KCMF1 from the hg19 reference. First, you will need to create a file that represents the sizes of all the chromosomes in the genome using [SAMtools](https://sourceforge.net/projects/samtools/files/samtools/):

```bash
samtools faidx hg19.fa
cut -f 1,2 hg19.fa.fai > hg19.sizes
```

```bash
samtools faidx Homo_sapiens_assembly38.fasta
cut -f 1,2 Homo_sapiens_assembly38.fasta.fai > GRCh38.sizes
```

#### Extend the upstream/downstream regions around KCMF1

Now we have the interval for KCMF1 in BED format we can simply use BEDTools to extend this by 1kb both upstream and downstream.

```bash
bedtools slop -i hg19_KCMF1.bed -g hg19.sizes -b 1000 > hg19_KCMF1.slop_1kb.bed
```

This should create a file with the following contents:

```bash
chr2    85197230    85287595    KCMF1   0   +
```

> NB: The intervals for EDIL3 were not extended so this step was omitted. 

#### Extract DNA sequence from reference

Finally, we can use BEDTools again to extract the DNA sequence for the KCMF1 and EDIL3 genes from the reference genome:

```bash
bedtools getfasta -name -fi hg19.fa -bed hg19_KCMF1.slop_1kb.bed > hg19_KCMF1.fa
```

```bash
bedtools getfasta -name -fi Homo_sapiens_assembly38.fasta -bed GRCh38_EDIL3.bed > GRCh38_EDIL3.fa
```

## Dataset origin

Fast5 files were obtained from the SG-NEx Project public dataset.

The SG-NEx project was initiated by the Genome Institute of Singapore with the aim to generate reference transcriptomes for 5 of the most commonly used cancer cell lines using Nanopore RNA-Seq data.

[Read more about SG-NEx](https://github.com/GoekeLab/sg-nex-data)

As outlined in the section above, there are 2 test-datasets - one each for barcoded and non-barcoded data, respectively.

The full datasets consisted of multiple subdirectories each containing hundreds of fast5 files.

A subset of 150 fast5 files were obtained from a randomly selected sub-directory using the following:

```bash
tar -xvf `ls | shuf -n1` # untar a random subdirectory
cd <untarred_dir>
cp `ls | shuf -n 150` <output_dir> # obtain random files
```
