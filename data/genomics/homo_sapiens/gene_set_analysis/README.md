# Test files for gene set analysis

## GSEA test inputs

These were downloaded from [the Broad's test suite](https://data.broadinstitute.org/gsea-msigdb/gsea/test_suite/gpunit/GSEA/v20/input) for the purpose of testing a module wrapping the Broad's GSEA tool:

 - P53_6samples_collapsed_symbols.gct: [gene cluster text file format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GCT:_Gene_Cluster_Text_file_format_.28.2A.gct.29) (essentially expression matrices) 
 - P53_6samples.cls: [Categorical (e.g tumor vs normal) class file format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CLS:_Categorical_.28e.g_tumor_vs_normal.29_class_file_format_.28.2A.cls.29) (describing the classes in .gct files)
 - c1.symbols.reduced.gmx: [GMX: Gene MatriX file format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29) (lists of gene sets) 
