# Complete Test Dataset

A complete plate of data can be downloaded from the Cell Painting Gallery using the AWS Command Line Interface.
See [AWS CLI documentation](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) for installation instructions.

This is a whole plate of data (plate `BR00117035`) from project `cpg0016-jump`, batch `2021_04_26_Batch1`.

The raw images: \
`aws s3 cp --recursive s3://cellpainting-gallery/cpg0016-jump/source_4/images/2021_04_26_Batch1/images/BR00117035__2021-05-02T16_02_51-Measurement1/ cpg0016-jump/source_4/images/2021_04_26_Batch1/images/BR00117035__2021-05-02T16_02_51-Measurement1/ --no-sign-request`

The pipelines: \
`curl https://raw.githubusercontent.com/broadinstitute/imaging-platform-pipelines/refs/heads/master/JUMP_production/JUMP_analysis_v3.cppipe -o cpg0016-jump/source_4/workspace/pipelines/2021_04_26_Batch1/analysis.cppipe --create-dirs`
`curl https://raw.githubusercontent.com/broadinstitute/imaging-platform-pipelines/refs/heads/master/JUMP_production/JUMP_illum_LoadData_v1.cppipe -o cpg0016-jump/source_4/workspace/pipelines/2021_04_26_Batch1/illum.cppipe --create-dirs`
`curl https://raw.githubusercontent.com/broadinstitute/imaging-platform-pipelines/refs/heads/master/JUMP_production/JUMP_segment_LoadData_v1.cppipe -o cpg0016-jump/source_4/workspace/pipelines/2021_04_26_Batch1/assaydev.cppipe --create-dirs`

The load_data.csv (note that the file paths will need to be edited to match the local download location): \
`aws s3 cp --recursive s3://cellpainting-gallery/cpg0016-jump/source_4/workspace/load_data_csv_orig/2021_04_26_Batch1/BR00117035/ cpg0016-jump/source_4/workspace/load_data_csv_orig/2021_04_26_Batch1/BR00117035/ --no-sign-request`

The illumination correction images: \
`aws s3 cp --recursive s3://cellpainting-gallery/cpg0016-jump/source_4/images/2021_04_26_Batch1/illum/BR00117035/ cpg0016-jump/source_4/images/2021_04_26_Batch1/illum/BR00117035/ --no-sign-request`

The assaydev outputs: \
`aws s3 cp --recursive s3://cellpainting-gallery/cpg0016-jump/source_4/workspace/assaydev/2021_04_26_Batch1/ cpg0016-jump/source_4/workspace/assaydev/2021_04_26_Batch1/ --exclude "*" --include "*BR00117035*" --no-sign-request`

The analysis outputs: \
`aws s3 cp --recursive s3://cellpainting-gallery/cpg0016-jump/source_4/workspace/analysis/2021_04_26_Batch1/BR00117035/ cpg0016-jump/source_4/workspace/analysis/2021_04_26_Batch1/BR00117035/ --no-sign-request`

The backends: \
`aws s3 cp --recursive s3://cellpainting-gallery/cpg0016-jump/source_4/workspace/backend/2021_04_26_Batch1/BR00117035/ cpg0016-jump/source_4/workspace/backend/2021_04_26_Batch1/BR00117035/ --no-sign-request`

The metadata: \
`aws s3 cp --recursive s3://cellpainting-gallery/cpg0016-jump/source_4/workspace/metadata/external_metadata/ cpg0016-jump/source_4/workspace/metadata/external_metadata/ --no-sign-request`
`aws s3 cp --recursive s3://cellpainting-gallery/cpg0016-jump/source_4/workspace/metadata/platemaps/2021_04_26_Batch1/ cpg0016-jump/source_4/workspace/metadata/platemaps/2021_04_26_Batch1/ --no-sign-request`

The profiles: \
`aws s3 cp --recursive s3://cellpainting-gallery/cpg0016-jump/source_4/workspace/profiles/2021_04_26_Batch1/BR00117035/ cpg0016-jump/source_4/workspace/profiles/2021_04_26_Batch1/BR00117035/ --no-sign-request`


If using this data, please cite the original creators of the data ([Chandrasekaran et al., 2023](https://doi.org/10.1101/2023.03.23.534023)) and the Cell Painting Gallery where the data is hosted ([Weisbart et al., 2024](https://doi.org/10.1038/s41592-024-02399-z)).

Chandrasekaran, S. N., Ackerman, J., Alix, E., Michael Ando, D., Arevalo, J., Bennion, M., Boisseau, N., Borowa, A., Boyd, J. D., Brino, L., Byrne, P. J., Ceulemans, H., Ch’ng, C., Cimini, B. A., Clevert, D.-A., Deflaux, N., Doench, J. G., Dorval, T., Doyonnas, R., … Carpenter, A. E. (2023). JUMP Cell Painting dataset: morphological impact of 136,000 chemical and genetic perturbations. In bioRxiv (p. 2023.03.23.534023). https://doi.org/10.1101/2023.03.23.534023

Weisbart, E., Kumar, A., Arevalo, J. et al. Cell Painting Gallery: an open resource for image-based profiling. Nat Methods 21, 1775–1777 (2024). https://doi.org/10.1038/s41592-024-02399-z
