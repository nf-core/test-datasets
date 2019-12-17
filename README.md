# test-datasets: `imcyto`

This branch contains test data to be used for automated testing with the [nf-core/imcyto](https://github.com/nf-core/imcyto) pipeline.

## Content of this repository

`inputs/`: Directory containing example input files required for pipeline in `.mcd`, `.txt` and `.tiff` format.  
`plugins/` : Cellprofiler and Ilastik pipeline files.  
`plugins/cp_plugins/` : Plugin files or scripts required to process CellProfiler pipeline modules.  

See [`usage.md`](https://github.com/nf-core/imcyto/blob/master/docs/usage.md) for a more detailed description of the parameters for the pipeline and [`output.md`](https://github.com/nf-core/imcyto/blob/master/docs/output.md) for further information regarding how to create your own CellProfiler/Ilastik pipeline files.
