# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
# test-datasets: `fastqrepair`

This branch contains test data to be used for automated testing with the [nf-core/fastqrepair](https://github.com/nf-core/fastqrepair) pipeline.

## Content of this repository
Test data to be used for automated testing with the nf-core pipelines

### Custom test-dataset

>`testdata/test_30reads_R[1/2].fastq.gz`: 30 paired-end reads with `R1` corrupted and `R2` not corrupted and where:
##### R1
- The quality string of the first read (`@NS500299:185:HK57NBGXG:1:11101:16144:1046 2:N:0:CGAGGCTG`) has two nucleaotides less than the sequence string
- Three are blank lines among the first three reads (not an error tough)
- The third read's name begins with "test " instead of "@" (`test @NS500299:185:HK57NBGXG:1:11101:18734:1046 2:N:0:CGAGGCTG`)
##### R2
- It does not contains the first four reads (i.e., `@NS500299:185:HK57NBGXG:1:11101:16144:1046 2:N:0:CGAGGCTG`, `@NS500299:185:HK57NBGXG:1:11101:19673:1046 2:N:0:CGAGGCTG`, `@NS500299:185:HK57NBGXG:1:11101:18734:1046 2:N:0:CGAGGCTG`, `@NS500299:185:HK57NBGXG:1:11101:16513:1046 2:N:0:CGAGGCTG`)
- Fifth and sixth reads are exchanged (if compared with `R1`'s reads order)

### [Bio Data Zoo](https://github.com/omgenomics/bio-data-zoo) test-dataset ([License](https://github.com/omgenomics/bio-data-zoo/blob/main/LICENSE))
>`testdata/quality_mismatch.fastq`: fastq where 2nd read has len(sequence) != len(quality)

>`testdata/truncated_clean.fastq`: fastq where 3rd read is truncated right after the sequence

>`testdata/truncated_halfway.fastq`: fastq where 2nd read is truncatd half-way through the sequence
