# 10X_Multiome_test_data
Test data for 10X_Multiome (cellranger-arc) analysis

## cellranger-arc mkgtf
1. mus_musculus chr19.filtered.gtf file is from Genode Release M31 (https://www.gencodegenes.org/mouse/). Then filtered to only genes on chromosome 19. It is at `data/genomics/mus_musculus/genome/`
2. Run test with command: 
```
cellranger-arc mkgtf chr19.filtered.gtf chr19.filtered2.gtf \
                   --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lincRNA \
                   --attribute=gene_biotype:antisense \
                   --attribute=gene_biotype:IG_LV_gene \
                   --attribute=gene_biotype:IG_V_gene \
                   --attribute=gene_biotype:IG_V_pseudogene \
                   --attribute=gene_biotype:IG_D_gene \
                   --attribute=gene_biotype:IG_J_gene \
                   --attribute=gene_biotype:IG_J_pseudogene \
                   --attribute=gene_biotype:IG_C_gene \
                   --attribute=gene_biotype:IG_C_pseudogene \
                   --attribute=gene_biotype:TR_V_gene \
                   --attribute=gene_biotype:TR_V_pseudogene \
                   --attribute=gene_biotype:TR_D_gene \
                   --attribute=gene_biotype:TR_J_gene \
                   --attribute=gene_biotype:TR_J_pseudogene \
                   --attribute=gene_biotype:TR_C_gene
```
3. Test runs successfully when output file is the same since only basic gene annotations were downloaded.

## cellranger-arc mkref
1. Mus musculus chr19.fa.gz file downloaded from UCSC Genome Browser website. Decompress the gzipped file with this command:
```
gzip -d data/genomics/mus_musculus/genome/chr19.fa.gz
```
2. Run command:
```
cellranger-arc mkref --config=data/generics/config/cellranger_arc_mkref_test_mm39_chr19_config.json
```
5. Test successful when program completes without error and a GRCm39 folder is created.

## cellranger-arc count
1. Test fastq data were downloaded from SRA database (file ids: SRR18907480 and SRR18907481). Starting 2 M reads from SRR18907480 and 10 M reads from SRR18907481 were initially mapped to mm39 chr19 and only aligned reads were filtered and converted back to fastq files. Files are at `data/genomics/mus_musculus/illumina/10xgenomics/multiome
2. Generate the libraries.csv input file 
```
echo fastqs,sample,library_type > libraries.csv
echo `pwd`/data/genomics/mus_musculus/illumina/10xgenomics/multiome,SRR18907480_chr19_sub,Gene Expression >> libraries.csv
echo `pwd`/data/genomics/mus_musculus/illumina/10xgenomics/multiome,SRR18907481_chr19_sub,Chromatin Accessibility >> libraries.csv
```

4. Run the command:
```
cellranger-arc count --id test_count --reference ./GRCm39 --libraries libraries.csv
```
3. Verify that the output test_count/outs/summary.csv is same as the data/generics/csv/cellranger_arc_count_output_summary.csv file.
