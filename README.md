# [nf-core/epitopeprediction](https://github.com/nf-core/epitopeprediction) test data

This branch contains test data and non-free software for the epitopeprediction pipeline.

## Input data

* `variants/`

  Input files to be provided to the pipeline via `--input`
* `proteins/`

  Alternative input files to be provided to the pipeline via `--proteins`
* `peptides/`

  Alternative input files to be provided to the pipeline via `--peptides`

## Additional data

* `alleles/`

  Allele lists to be provided to the pipeline via `--alleles`

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
