## Data Description

- prestitched
  - 'lunaphore_example.ome.tiff'
    - exemplary stitched OME-TIFF file
    - subset of the `sample_control.r2.SeqIF.tiff` file accessed at https://www.synapse.org/Synapse:syn51471704 from the Wuennemann et al. 2024 study.
    - the file includes the first 3 channels (DAPI, Cy5, TRITC), and the last two channels (TNNT2, WGA) and was subsetted with input_image[[0,1,2,18,19], 15500:16500, 7200:8200]
    - pixel size: 0.23um
  - 'markers-lunaphore_example.csv'
    - marker sheet with correct values for the `lunaphore_example.ome.tiff` file
