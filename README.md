# test-datasets: `nanoseq`

This branch contains test data to be used for automated testing with the [nf-core/nanoseq](https://github.com/nf-core/nanoseq) pipeline.

## Content of this repository

`samplesheet.csv`: Sample information file
`fast5/` : subset of fast5 files from direct cDNA Nanopore reads for MCF7 (Breast Cancer) cell line
`reference/`: Genome reference files (iGenomes `GRCh37` Ensembl release; region `` of chromosome `1` only)

## Dataset origin

Fast5 files were obtained from the SG-NEx Project public dataset.

The SG-NEx project was initiated by the Genome Institute of Singapore with the aim to generate reference transcriptomes for 5 of the most commonly used cancer cell lines using Nanopore RNA-Seq data.

[Read more about SG-NEx](https://github.com/GoekeLab/sg-nex-data)

### Sample information

The data is not barcoded.

| Type	                |   Sample	  |
|-----------------------|-------------|
| Flow Cell       	    | FLO-MIN106	|
| Kit	                  | SQK-DCS108	|
| Barcode Kit	          | None    	  |
