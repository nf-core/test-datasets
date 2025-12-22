# [nf-core/epitopeprediction](https://github.com/nf-core/epitopeprediction) test data

This branch contains test data and non-free software for the epitopeprediction pipeline.

## Input data

  The pipeline accepts three sources of inputs specified in the `filename` column of the samplesheet to predict (neo-)epitopes:
* `variants/`

  `SnpEff`/`VEP`-annotated `vcf` file containing information about SNPs and indels. Mandatory to predict neo-epitopes.

* `proteins/`

  Fasta file containing protein sequences. Each peptide of the protein is considered for MHC-binding prediction.

* `peptides/`

  Tab-separated file with mandatory header `id` and `sequence` containing peptides considered for MHC-binding prediction.

## Additional data

* `alleles/`

  Allele list or `.txt` file to be provided in the `alleles` column of the samplesheet.

* `biomart_dumps/`

  CSV/TSV dumps of the Ensembl BioMart data. The dump file can be created by querying the [Ensembl Biomart](https://www.ensembl.org/biomart/martview/) for the relevant database and dataset (e.g. `grch37` or `grch38`) and selecting the attributes Protein stable ID (`ensembl_peptide_id`), RefSeq peptide ID (`refseq_peptide`), UniProtKB/Swiss-Prot ID (`uniprotswissprot`), Transcript stable ID (`ensembl_transcript_id`). You can select other genome versions as described above. A list of currently available archives can be found [here](https://www.ensembl.org/info/website/archives/index.html?redirect=no).

## Software

* `software/non-free-software.tar.gpg`

  An encrypted tarball providing non-free software supported by the pipeline to be used in CI tests. The tarball contains the original compressed tarballs for each tool as provided upstream.
  ```
  $ gpg -d software/non-free-software.tar.gpg | tar tv
  gpg: AES256.CFB encrypted data
  gpg: encrypted with 1 passphrase
  -rw-r--r--  0 lk     staff 3624305 Jan 13 13:10 netMHC-4.0a.Linux.tar.gz
  -rw-r--r--  0 lk     staff 8034930 Jan 14 17:49 netMHCII-2.2.Linux.tar.Z
  -rw-r--r--  0 lk     staff 3251741 Jan 15 20:05 netMHCIIpan-3.1a.Linux.tar.gz
  -rw-r--r--  0 lk     staff 1332999 Jan 15 16:18 netMHCpan-4.0a.Linux.tar.gz
  lrwxr-xr-x  0 lk     staff       0 Jan 22 13:06 netmhc.tar.gz -> netMHC-4.0a.Linux.tar.gz
  lrwxr-xr-x  0 lk     staff       0 Jan 22 13:06 netmhcii.tar.Z -> netMHCII-2.2.Linux.tar.Z
  lrwxr-xr-x  0 lk     staff       0 Jan 22 13:06 netmhciipan.tar.gz -> netMHCIIpan-3.1a.Linux.tar.gz
  lrwxr-xr-x  0 lk     staff       0 Jan 22 13:06 netmhcpan.tar.gz -> netMHCpan-4.0a.Linux.tar.gz
  ```
  The passphrase is included as a repository level secret in the pipeline repository.
