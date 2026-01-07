## Where does the data come from?

TOYv1 is subsetted from the dose-response screen of CTRPv2 (37 drugs, 90 cell lines). 
TOYv2 is subsetted from the dose-response screen of GDSC2 (37 drugs, 90 cell lines, 80 cell lines and 32 drugs overlap with TOYv1).

More information can be found on the [preprocess GitHub connected to drevalpy](https://github.com/daisybio/preprocess_drp_data/blob/main/Toy_Data/create_toy.py).

The omics data for TOYv1 is taken from CCLE and heavily subsetted. The preprocessing information for it is at [the GitHub](https://github.com/daisybio/preprocess_drp_data/blob/main/CCLE/README.md).
The omics data for TOYv2 is taken from GDSC and heavily subsetted: [preprocessing information](https://github.com/daisybio/preprocess_drp_data/blob/main/GDSC/README.md).

### Original sources (before preprocessing):

#### TOYv1

|                                                                    | Source                                                                                                                                                                                                                     |
|--------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Gene expression                                                    | [Ghandi et al. in 2019](doi.org/10.1038/s41586-019-1186-3), BioProject PRJNA523380                                                                                                                                         |
| Methylation                                                        | From [DepMap](https://depmap.org/portal/data_page/?tab=allData): Methylation (RRBS), file: CCLE_RRBS_TSS_CpG_clusters_20180614.txt.                                                                                        |
| Mutations                                                          | From [Sanger Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/downloads): Mutation Data -> Mutations All.                                                                                                     | 
| Copy number variation                                              | From [Sanger Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/downloads): Copy Number Data -> Copy Number (SNP6) -> PICNIC absolute copy numbers and GISTIC scores derived from Affymetrix SNP6.0 array data. |
| Proteomics                                                         | From SangerCellModelPassports ([Gonçalves et al. (2021)](https://www.sciencedirect.com/science/article/pii/S1535610822002744))                                                                                             |
| Response                                                           | From [the CTD^2 portal](https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/)                                                                                                                    |
| Drug SMILES (-> used for fingerprints, graphs, MolGNet, ChemBERTa) | From PubChem. Fingerprints, graphs, MolGNet embeddings can be created with the [featurizers](https://github.com/daisybio/drevalpy/tree/development/drevalpy/datasets/featurizer)                                           |
| DIPK_features                                                      | MolGNet see above. gene_list_sel.txt, human_ppi_features.tsv, key_genes.txt from the [DIPK Google Drive](https://drive.google.com/drive/folders/16hP48-noHi3-c_LP9TcZxkwAzqxgR0VB)                                         |

#### TOYv2

|                                                                    | Source                                                                                                                                                                                       |
|--------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Gene expression                                                    | From [GDSC1000 resources](https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html): RMA normalised expression data for cell-lines (Cell_line_RMA_proc_basalExp.txt)            |
| Methylation                                                        | From [GDSC1000 resources](https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html): Pre-Processed beta values for all CpG islands across all the cell-line (F2_METH_CELL_DATA) |
| Mutations                                                          | Identical to TOYv1.                                                                                                                                                                          | 
| Copy number variation                                              | Identical to TOYv1.                                                                                                                                                                          |                                                                                                                                                                                   |
| Proteomics                                                         | Identical to TOYv1.                                                                                                                                                                          |
| Response                                                           | From [the GDSC website](https://www.cancerrxgene.org/downloads/bulk_download): GDSC2-raw-data.                                                                                               |
| Drug SMILES (-> used for fingerprints, graphs, MolGNet, ChemBERTa) | From PubChem, see above.                                                                                                                                                                     |
| DIPK_features                                                      | MolGNet see above.                                                                                                                                                                           |

#### meta/gene_lists

* landmark_genes_reduced.csv: Originally, the 978 landmark genes are from the L1000 assay. Here, just 12 genes occurring in all omics screens.
* drug_target_genes_all_drugs(_proteomics).csv: Originally, the drug target genes are the genes targeted by the drugs used in GDSC, extractable from the GDSC Data Portal (compounds annotation). Here, just the same 12 genes occurring in all omics screens.
* *_intersection.csv: Originally, the omics features of all screens intersecting. Here, just the same 12 genes as above, except for the methylation file. The methylation file just contains all methylation features from TOYv1.intersection(TOYv2).  
