# test-datasets: `pixelator`

This branch contains test data to be used for automated testing with the [nf-core/pixelator](https://github.com/nf-core/pixelator) pipeline.

## Content of this repository

This repository contains data for the Pixelgen Technologies Molecular Pixelation (PNA) and
Proximity Network Technology(PNA) assays. The data for the respective assay lives in the
subdirectories `mpx` and `pna`.

- `mpx/testdata/` :
    - `mpx/modules` : Testdata for pipeline local module tests
    - `mpx/fastq` : Input fastq files for pipeline level tests
- `mpx/panels/`: panel files to test passing custom panels

- `mpx/samplesheet/samplesheet.csv`: Experiment design file for minimal test dataset.
- `mpx/samplesheet/samplesheet_full.csv`: Experiment design file for full test dataset.
- `mpx/samplesheet/samplesheet_v2.csv`: Experiment design file for minimal test dataset using the v2 panel.
- `mpx/samplesheet/samplesheet_mpx_scsp_v1_immunology1.csv` : Experiment design file for minimal test dataset using the same data as the local module tests
