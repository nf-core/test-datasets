# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
This branch contains data to be used for automated testing with the [nf-core/cageseq](https://github.com/nf-core/cageseq) pipeline.

## Content of this repository

`reference`: ChrIII of the sacCer3 (*Saccharomyces cerevisiae*) genome assembly.

`testdata/{S1,S2}.fastq.gz`: Single-end CAGE reads uniquely mapping to chrIII of sacCer3 (two test samples).

## Test dataset generation

### Dependencies

To generate reference and FASTQ files, you will need the following tools:

```
awk
faToTwoBit
gzip
samtools
seqtk
STAR
tar
wget
```

You will also need R and the `BSgenome` R package.

### Reference files

1. To generate `sacCer3chrIII.fa`, obtain [`chromFa.tar.gz`](https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz) (*Saccharomyces cerevisiae* genome assembly, split by chromosome), unpack it and retain only `chrIII.fa`, renaming it into `sacCer3chrIII.fa`:

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
tar xvzf chromFa.tar.gz
mv chrIII.fa sacCer3chrIII.fa
rm chr*
```

2. To generate `sacCer3chrIII.gtf`, obtain [`sacCer3.ensGene.gtf.gz`](https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz), unpack it and select only chrIII:

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz
gzip -d sacCer3.ensGene.gtf.gz
awk -F'\t' '$1 == "chrIII"' sacCer3.ensGene.gtf > sacCer3chrIII.gtf
rm sacCer3.ensGene.gtf
```

3. To generate `BSgenome.sacCer3ChrIII.Test.sacCer3ChrIII_1.0.tar.gz`, convert `sacCer3chrIII.fa` into the `2bit` format (it will retain the FASTA file) and move the obtained `2bit` file to a dedicated directory:

```
faToTwoBit sacCer3chrIII.fa sacCer3chrIII.2bit
mkdir seq_dir
mv sacCer3chrIII.2bit seq_dir/
```

Then, save the following text

```
Package: BSgenome.sacCer3ChrIII.Test.sacCer3ChrIII
Title: Full chrIII from sacCer3
Description: Full chrIII from sacCer3
Version: 1.0
organism: Saccharomyces cerevisiae
common_name: Yeast
provider: custom
genome: sacCer3ChrIII
release_date: April 2011
source_url: https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/
organism_biocview: Saccharomyces_cerevisiae
BSgenomeObjname: sacCer3ChrIII
SrcDataFiles: sacCer3chrIII.2bit
seqs_srcdir: /path/to/seq_dir
seqfile_name: sacCer3chrIII.2bit
circ_seqs: character(0)
```

into `BSgenome.sacCer3ChrIII.Test.sacCer3ChrIII-seed`, replace `/path/to/seq_dir` with the actual path to `seq_dir` on you machine, open an R session and forge the BSgenome package:

```
library(BSgenome)
forgeBSgenomeDataPkg("BSgenome.sacCer3ChrIII.Test.sacCer3ChrIII-seed")
```

This command will produce the package directory `BSgenome.sacCer3ChrIII.Test.sacCer3ChrIII`. Next, exit R and build the package from the command line, ignoring warnings and notes:

```
R CMD build BSgenome.sacCer3ChrIII.Test.sacCer3ChrIII
R CMD check BSgenome.sacCer3ChrIII.Test.sacCer3ChrIII_1.0.tar.gz
```

In this way, you will obtain and validate `BSgenome.sacCer3ChrIII.Test.sacCer3ChrIII_1.0.tar.gz`.

4. Additionally, generate a STAR index for sacCer3ChrIII to map CAGE reads for further selection of a test subsample:

```
mkdir sacCer3chrIII_star_index
STAR --runMode genomeGenerate \
     --genomeFastaFiles sacCer3chrIII.fa \
     --sjdbGTFfile sacCer3chrIII.gtf \
     --genomeDir sacCer3chrIII_star_index \
     --genomeSAindexNbases 8 \
     --runThreadN 4 
```

You will need the `sacCer3chrIII_star_index` to generate test FASTQ files (see below).

### FASTQ files

*S. cerevisiae* single-end CAGE datasets were obtained from the following publication:

> Börlin, C. S., Cvetesic, N., Holland, P. et al. *Saccharomyces cerevisiae* displays a stable transcription start site landscape in multiple conditions. FEMS Yeast Research, 19, foy128 (2019). https://doi.org/10.1093/femsyr/foy128. [PubMed:30590648](https://pubmed.ncbi.nlm.nih.gov/30590648/). [ArrayExpress:E-MTAB-6650](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-6650?query=E-MTAB-6650).

Sampling information:

| Project Accession | Sample Accession | Experiment Accession | Run Accession | Original Read Count | Subsampled Read Count | Sample Condition                                       |
|-------------------|------------------|----------------------|---------------|---------------------|-----------------------|--------------------------------------------------------|
| PRJEB25845 | SAMEA1054240     | ERX2514501           | ERR2495148    | 9,637,094           | 10,000 | anaerobic conditions, replicate 1 |
| PRJEB25845 | SAMEA1054242     | ERX2514503           | ERR2495150    | 5,295,517           | 10,000 | ethanol limitation, replicate 1 |

To obtain the test FASTQ files, follow the procedure below.

Suppose we generate the FASTQ files in the same directory as the reference files.

1. Download and unpack the original FASTQ files:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR249/ERR2495148/20171106_CAGE_SequencingData_Ana_Rep1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR249/ERR2495150/20171106_CAGE_SequencingData_Eth_Rep1.fastq.gz
gzip -d 20171106_*.gz
```

2. Map the `Ana_Rep1` sample to the yeast chrIII (which may take several minutes) and retain only the BAM file, renaming it into `S1_initial.bam`:

```
runThreadN=32 # set to the number of available cores
STAR \
    --readFilesIn 20171106_CAGE_SequencingData_Ana_Rep1.fastq \
    --genomeDir sacCer3chrIII_star_index \
    --alignEndsType Local \
    --outSAMtype BAM Unsorted \
    --runThreadN $runThreadN
rm Log.* SJ.out.tab
mv Aligned.out.bam S1_initial.bam
```

3. Map the `Eth_Rep1` sample to the yeast chrIII (which may take several minutes) and retain only the BAM file, renaming it into `S2_initial.bam`:

```
runThreadN=32 # set to the number of available cores
STAR \
    --readFilesIn 20171106_CAGE_SequencingData_Eth_Rep1.fastq \
    --genomeDir sacCer3chrIII_star_index \
    --alignEndsType Local \
    --outSAMtype BAM Unsorted \
    --runThreadN $runThreadN
rm Log.* SJ.out.tab
mv Aligned.out.bam S2_initial.bam
```

4. Obtain only uniquely mapped alignments, dump the corresponding reads into FASTQ files, choose only reads starting with a `G` (to maximise the relevant signal) and randomly choose 10,000 test reads per sample:

```
for sample in S1 S2; do
    samtools view -b -q 255 "${sample}_initial.bam" | \
    samtools fastq - > "${sample}_unique.fastq"
    <${sample}_unique.fastq \
    awk 'NR % 4 == 1 { h = $0 } \
            NR % 4 == 2 {sq = $0} \
            NR % 4 == 3 { s = $0 } \
            NR % 4 == 0 && sq ~ /^G/ {print h "\n" sq "\n" s "\n" $0}' > \
    ${sample}_unique_G.fastq
    seqtk sample -s 42 "${sample}_unique_G.fastq" 10000 > "${sample}.fastq"
done
```

Then gzip the obtained FASTQ files and remove the intermediate ones:

```
gzip S1.fastq S2.fastq
rm S1_initial.bam S1_unique.fastq S1_unique_G.fastq 
rm S2_initial.bam S2_unique.fastq S2_unique_G.fastq
```

The obtained `S1.fastq.gz` and `S2.fastq.gz` are the final test datasets.
