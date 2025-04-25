# test-datasets: mcmicro

This branch contains test data to be used for automated testing with the [nf-core/mcmicro](https://github.com/nf-core/mcmicro) pipeline.

## Content of this repository
`samplesheets/markers-test.csv`: Markersheet file for minimal test<br>
`samplesheets/samplesheet-test.csv`: Samplesheet file for minimal test<br>
`samplesheets/markers-test_full.csv`: Markersheet file for full test<br>
`samplesheets/samplesheet-test_full.csv`: Samplesheet file for full test<br>

### Test dataset origin

The image files referenced in the samplesheets are human tonsil section 5 μm thick, imaged in four-channel immunofluorescence over three rounds of bleaching and re-staining with different antibodies. Published in https://doi.org/10.1038/s41592-021-01308-y (slide WD-75684-02). These image files each contain a 3 x 3 grid of overlapping four-channel image tiles with dimensions 220 x 180 pixels. The tiles were cropped from a single raw tile taken from one of the original data files. A Python script (https://github.com/labsyspharm/nf-core-test-datasets-build/blob/main/tonsil-cycif/tonsil-cycif.py) was written to extract the synthetic tile fields from an "interesting" region of the single raw tile, generate stage position metadata with random perturbations to approximate the behavior of a real-world microscope stage, and save the images and metadata as an OME-TIFF file.
