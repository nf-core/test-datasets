# test-datasets: `nanoseq`

This branch contains test data to be used for automated testing with the [nf-core/nanoseq](https://github.com/nf-core/nanoseq) pipeline.

## Content of this repository

| File	                | Description	                                                                                              |
|-----------------------|-----------------------------------------------------------------------------------------------------------|
| `samplesheet.csv`     | Sample information file                                                                                   |
| `fast5/`              | Subset of fast5 files from direct cDNA Nanopore reads for HepG2 (Liver Cancer) and K562 (Leukemia) cell lines                  |
| `reference/`          | Genome reference files (`hg19` UCSC release; cDNA for KCMF1 gene +- 1kb obtained via UCSC Table Browser)  |

## Dataset origin

Fast5 files were obtained from the SG-NEx Project public dataset.

The SG-NEx project was initiated by the Genome Institute of Singapore with the aim to generate reference transcriptomes for 5 of the most commonly used cancer cell lines using Nanopore RNA-Seq data.

[Read more about SG-NEx](https://github.com/GoekeLab/sg-nex-data)

A subset of 150 fast5 files was obtained using the following:
```
tar -xvf `ls | shuf -n1` # to obtain a random directory
cd <untarred_dir>
cp `ls | shuf -n 150` <output_dir> # to obtain random files
```

### Sample information

|     	                |         	  |
|-----------------------|-------------|
| Flow Cell       	    | FLO-MIN106	|
| Kit	                  | SQK-DCS109	|
| Barcode Kit	          | EXP-NBD103	  |
