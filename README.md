# nfcore/test-datasets: `rifseq`

This branch contains test data to be used for automated testing with the
[nf-core/rifseq](https://github.com/nf-core/rifseq) pipeline.

## Contents of this repository

The test dataset for the *nf-core/rifseq* pipeline contains two minimal RIF-Seq
plates with four barcodes each, a sample sheet, metadata with per-well sample
information, a cDNA FASTA reference of the human chromosome 22 as well as
a minimally subset rRNA FASTA database. The test data comes from a real RIF-Seq
dataset subset to 8000 reads per plate (2000 reads per sample per plate), with
an additional 10 non-RIF-Seq reads per sample and plate. The repository is
divided into the following directories and files:

`testdata/` <br>
    - `test_plate_{A,B}.fastq.gz`: Two plates of RIF-Seq test data <br>
    - `samplesheet.tsv`: The sample sheet <br>
    - `metadata.tsv`: The metadata <br>
`reference/` <br>
        - `Homo_sapiens.GRCh38.cdna.chr22.fa.gz`: Human chromosome 22 cDNA FASTA <br>
        - `silva-euk-18s-id95-minimal.fasta`: Minimally subset rRNA FASTA for SortMeRNA <br>
        - `rrna-db-manifest.fasta`: File manifest for the rRNA FASTA <br>
`script/` <br>
    - `simulate-rifseq-data.py`: Script for generating RIF-Seq test data

### Sample sheet

| plate        | fastq                                  |
| -------------|----------------------------------------|
| test_plate_A | [github address]/test_plate_A.fastq.gz |
| test_plate_B | [github address]/test_plate_B.fastq.gz |

### Metadata

| plate_id     | barcode_id | treatment  | dose  |
|--------------|------------|------------|-------|
| test_plate_A | 1          | DMSO       | 0.001 |
| test_plate_A | 97         | DMSO       | 0.001 |
| test_plate_A | 193        | Vorinostat | 8.3   |
| test_plate_A | 289        | Vorinostat | 8.3   |
| test_plate_B | 1          | DMSO       | 0.001 |
| test_plate_B | 97         | DMSO       | 0.001 |
| test_plate_B | 193        | Vorinostat | 8.3   |
| test_plate_B | 289        | Vorinostat | 8.3   |

## Dataset origin

The test data was generated using the `simulate-rifseq-data.py` Python script,
along with the `Homo_sapiens.GRCh38.cdna.chr22.fa.gz` FASTA reference, like so:

```bash
python script/simulate-rifseq-data.py \
    reference/Homo_sapiens.GRCh38.cdna.chr22.fa \
    testdata/test_plate_A.fastq \
    --seed 42
python script/simulate-rifseq-data.py \
    reference/Homo_sapiens.GRCh38.cdna.chr22.fa \
    testdata/test_plate_B.fastq \
    --seed 24
```

The FASTA reference was generated from the `Homo_sapiens.GRCh38.cdna.all.fa`
file available at [Ensembl's FTP site](ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/cdna/)
and the following code:

```bash
grep -A 1 "chromosome:GRCh38:22" Homo_sapiens.GRCh38.cdna.all.fa \
    | grep -v -- -- \
    > reference/Homo_sapiens.GRCh38.cdna.chr22.fa
```

The rRNA FASTA reference was generated from the `silva-euk-18s-id95.fasta`
available at the [SortMeRNA GitHub](https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases)
and the following code:

```bash
head -100 silva-euk-18s-id95.fasta > reference/silva-euk-18s-id95-minimal.fasta
```
