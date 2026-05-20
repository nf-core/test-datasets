# PacBio Test Data

This directory contains PacBio test datasets used for local pipeline testing.

## `revio-with-kinetics.bam`

`revio-with-kinetics.bam` is a PacBio long-read Fiber-seq test dataseti, downsized (10 reads)
from the[`fiberseq/fibertools-rs`](https://github.com/fiberseq/fibertools-rs)
repository. It contains kinetic signatures and no m6A tags (A+a). 


Upstream source:
[`tests/data/rebio.bam`](https://github.com/fiberseq/fibertools-rs/blob/main/tests/data/rebio.bam)

## `revio-with-m6a-tags.bam`

`revio-with-m6a-tags.bam` is a PacBio long-read Fiber-seq 100-read subset of the
source dataset in GSE330647, obtaining m6A modification tags (A+a).

