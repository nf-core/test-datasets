This folder holds test datasets for modules (DemuxEM, GMM-Demux, Seurat/HTO-Demux, Seurat/MULTI-seq, Hashed Drops, CellHashR/BFF and HashSolo) that deconvolve single-cell multiplexed datasets using the hashing demultiplexing method. 
The dataset provided by CellHashR/BFF contains two datasets for HTO and RNA data in different formats (10x mtx, h5 and csv) as required for the different tools.
The dataset is part of the Barnyard Dataset, which can be obtained from GEO with SRA ID: SRR8890636. More information can be found at https://doi.org/10.1093/bioinformatics/btac213, section: Data sources.
For compatibility and testing within nf-core pipelines, the HTO 10x mtx was transformed into CSV as required by DemuxEM and the 10 mtx was filtered using CellHashR (“ProcessCountMatrix”) since only two barcodes were whitelisted for the experiment. 
More information can be found in the CellhashR/BFF repository: https://github.com/BimberLab/cellhashR/tree/master


