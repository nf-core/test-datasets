# test-datasets: evexplorer


# Content of this repository
## ev_small_rna
Contains small-RNA of extracellular vesicles obtained prostate cancer cell lines. Permission for using this dataset has been obtained by Dr. Rodney Ouellette, who's team initially generated the data. The data was generated from cell line media with the vn96 peptide pulldown of Extracellular Vesicles. Total RNA was isolated using the miRvana kit and sequenced using a Trilink library preparation for smallRNA and ion proton sequencers.

## Sub-sampling
The fastq files were subsampled using a ratio of 0.10 using `seqkit`

	parallel 'seqkit sample -p 0.10 -s 2020 {} | pigz > {.}.sampled.fastq.gz' ::: HPREC-Batch*.fastq.gz LNCAP-Batch*.fastq.gz
