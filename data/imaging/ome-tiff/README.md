# Test Data

## Cyclic immunofluorescence images of human tonsil

This data is intended for use in image processing tasks involving immunofluorescence imaging, whole-slide imaging, and cyclic imaging. It is a small and manageable semi-synthetic multi-field-of-view dataset derived from real images -- the image fields are smaller cropped windows from the original image data, but the pixel values themselves have not been changed. Single image fields can be used alone, or stitching and registration could be applied to assemble a full-size image.

### cycif-tonsil-cycle*.ome.tif

Human tonsil section 5 μm thick, imaged in four-channel immunofluorescence over three rounds of bleaching and re-staining with different antibodies. Published in https://doi.org/10.1038/s41592-021-01308-y (slide WD-75684-02). These image files each contain a 3 x 3 grid of overlapping four-channel image tiles with dimensions 220 x 180 pixels. The tiles were cropped from a single raw tile taken from one of the original data files. A Python script (https://github.com/labsyspharm/nf-core-test-datasets-build/blob/main/tonsil-cycif/tonsil-cycif.py) was written to extract the synthetic tile fields from an "interesting" region of the single raw tile, generate stage position metadata with random perturbations to approximate the behavior of a real-world microscope stage, and save the images and metadata as an OME-TIFF file.

### cycif-tonsil-dfp.ome.tif, cycif-tonsil-ffp.ome.tif

Dark-field (dfp) and flat-field (ffp) illumination correction profiles for cycif-tonsil-cycle*.ome.tif. The dfp profile is all zeros and the ffp profile is all ones and thus will not modify the image when applied, but the images are the correct shape for this data.

### cycif-tonsil-channels.csv

A three-column table containing channel number, imaging cycle number, and marker (target of the fluorescently labeled antibody or other stain) for each channel in each imaging cycle, useful for downstream tools expecting channel annotation.
