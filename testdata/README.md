
## Make aligned + unaligned bam / tgz files

### Subset bam to just *Ptprc* expression

```
(samtools)
 Wed  8 Apr - 06:12  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  samtools view -bh MACA_18m_M_BAT_52__possorted_genome_bam.bam chr1:138,063,407-138,175,306 > MACA_18m_M_BAT_52__possorted_genome_bam__ptprc.bam
```

### Get unique cell barcodes in PTPRC data

```
(samtools)
 Wed  8 Apr - 06:21  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  samtools view MACA_18m_M_BAT_52__possorted_genome_bam__ptprc.bam | cut -f 17- | grep -Eo '(CB|CB):Z:(([ACGT]+)(\\-1)?)' | sort | uniq > ptprc_cell_barcodes.txt
(samtools)
 ✘  Wed  8 Apr - 06:22  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  head ptprc_cell_barcodes.txt
CB:Z:AAAGTAGGTGTGAATA
CB:Z:AAAGTAGGTTGGTTTG
CB:Z:AAATGCCCAGCGAACA
CB:Z:AAATGCCGTTACCGAT
CB:Z:AACACGTAGCTAAGAT
CB:Z:AACTCCCTCAAACGGG
CB:Z:AACTCTTGTTCTGAAC
CB:Z:AACTTTCTCCCATTTA
CB:Z:AACTTTCTCGAATCCA
CB:Z:AAGACCTAGCGGATCA
```

### Get just the barcode sequence, without the "CB" and such

```
(samtools)
 Wed  8 Apr - 06:22  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  cut -f 3 -d ':' ptprc_cell_barcodes.txt > ptprc_cell_barcodes__barcode_sequence.txt
(samtools)
 Wed  8 Apr - 06:22  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  head ptprc_cell_barcodes__barcode_sequence.txt
AAAGTAGGTGTGAATA
AAAGTAGGTTGGTTTG
AAATGCCCAGCGAACA
AAATGCCGTTACCGAT
AACACGTAGCTAAGAT
AACTCCCTCAAACGGG
AACTCTTGTTCTGAAC
AACTTTCTCCCATTTA
AACTTTCTCGAATCCA
AAGACCTAGCGGATCA
```


### Grep for ptprc barcodes in sorted umi per cell barcode counts


```
 Wed  8 Apr - 06:26  /mnt/data_sm/olga/tabula-microcebus/analyses/kmermaid/tenx-tgz--tabula-muris-senis/10x-fastqs/umis-per-cell 
 olga@lrrr  fgrep --file /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b/ptprc_cell_barcodes__barcode_sequence.txt MACA_18m_M_SCAT_52__aligned__n_umi_per_cell.csv | sort -k2n -t',' -
```

### Use `tail` to just get the top 10 cells in the ptprc data

```
 Wed  8 Apr - 06:26  /mnt/data_sm/olga/tabula-microcebus/analyses/kmermaid/tenx-tgz--tabula-muris-senis/10x-fastqs/umis-per-cell 
 olga@lrrr  fgrep --file /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b/ptprc_cell_barcodes__barcode_sequence.txt MACA_18m_M_SCAT_52__aligned__n_umi_per_cell.csv | sort -k2n -t',' | tail
ACCTTTACAGACGTAG,65151
CGTCAGGAGCTAGTCT,66979
TCAATCTGTTGGACCC,67262
TTTCCTCCAGCGATCC,67799
ATTCTACCATTTCACT,71779
GTAGTCAAGACGCAAC,73761
AAAGTAGGTTGGTTTG,81051
TTGCCGTGTACCGGCT,86128
GGGTTGCAGCTGCCCA,89491
GTCACAATCGAGGTAG,101820
```

Output to file

```
 Wed  8 Apr - 06:29  /mnt/data_sm/olga/tabula-microcebus/analyses/kmermaid/tenx-tgz--tabula-muris-senis/10x-fastqs/umis-per-cell 
 olga@lrrr  fgrep --file /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b/ptprc_cell_barcodes__barcode_sequence.txt MACA_18m_M_SCAT_52__aligned__n_umi_per_cell.csv | sort -k2n -t',' | tail > /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b/top_barcodes_umi_counts.csv
```

### Convert top barcodes to 10x style barcodes

Take only the first column

```
(samtools)
 Wed  8 Apr - 06:38  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  cut -f 1 -d, top_barcodes_umi_counts.csv > top_barcodes__sequence.csv
```

Add CB: and -1 to each barcode seqauence

```
 Wed  8 Apr - 08:10  /mnt/data_sm/olga/tabula-microcebus/analyses/kmermaid/tenx-tgz--tabula-muris-senis/10x-fastqs/umis-per-cell 
 olga@lrrr  fgrep --file /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b/ptprc_cell_barcodes__barcode_sequence.txt MACA_18m_M_SCAT_52__aligned__n_umi_per_cell.csv | sort -k2n -t',' > /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b/ptprc_barcodes_umi_counts.csv
 Wed  8 Apr - 08:10  /mnt/data_sm/olga/tabula-microcebus/analyses/kmermaid/tenx-tgz--tabula-muris-senis/10x-fastqs/umis-per-cell 
 olga@lrrr  head /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b/ptprc_barcodes_umi_counts.csv
CCATGTCCAAGCGATG,1
CATCGGGCACATCCGG,2
GGTGTTACATTTGCTT,2
GCTCTGTTCCTTAATC,3
GTCCTCACAAGGTTTC,3
CAGCAGCGTCAAGCGA,4
ACCGTAACATATGAGA,5
TTAGTTCCATGTCCTC,5
TGCTACCTCCAAACAC,7
TAAGAGAGTCTCCACT,11
```

### Get only unaligned reads for these cells

Save unaligned reads

```
(samtools)
 Wed  8 Apr - 07:29  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  time samtools view -f4 MACA_18m_M_BAT_52__possorted_genome_bam.bam > unaligned.sam              samtools view -f4 MACA_18m_M_BAT_52__possorted_genome_bam.bam > unaligned.sam  160.77s user 4.46s system 91% cpu 2:59.71 total
```

fgrep the unaligned file for the ptprc cell barcodes

```
(samtools)
 Wed  8 Apr - 08:12  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  fgrep --file ptprc_cell_barcodes.txt unaligned.sam > unaligned__ptprc_cell_barcodes.sam
```


### Add header to unaligned reads from only cell barcodes

```
(samtools)
 Wed  8 Apr - 11:11  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  samtools view -H MACA_18m_M_BAT_52__possorted_genome_bam.bam > header.sam
(samtools)
 Wed  8 Apr - 11:11  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  samtools reheader header.sam unaligned__ptprc_cell_barcodes__no_header.bam > unaligned__ptprc_cell_barcodes.bam
```


### Merge aligned and unaligned bams

```
(samtools)
 Wed  8 Apr - 11:12  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  samtools merge mouse_brown_fat_ptprc_plus_unaligned.bam MACA_18m_M_BAT_52__possorted_genome_bam__ptprc.bam unaligned__ptprc_cell_barcodes.bam
```


### Samtools index

```
(samtools)
 Wed  8 Apr - 11:16  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  samtools index mouse_brown_fat_ptprc_plus_unaligned.bam
```

### Filter 10x barcodes.tsv file for Ptprc-aligning ones

```
(samtools)
 Wed  8 Apr - 11:19  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  fgrep --file ptprc_cell_barcodes__barcode_sequence.txt MACA_18m_M_BAT_52__barcodes.tsv > mouse_brown_fat_ptprc_plus_unaligned__barcodes.tsv
(samtools)
 Wed  8 Apr - 11:19  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  head mouse_brown_fat_ptprc_plus_unaligned__barcodes.tsv                             AAAGTAGGTGTGAATA-1
AAAGTAGGTTGGTTTG-1
AAATGCCCAGCGAACA-1
AAATGCCGTTACCGAT-1
AACTCCCTCAAACGGG-1
AACTTTCTCCCATTTA-1
AAGGAGCCACATAACC-1
ACACCAAGTCTCTTAT-1
ACACCCTTCGCCAGCA-1
ACACCGGGTTGCGCAC-1
```

#### Use first 1000 lines of unaligned reads

```
(samtools)
 Wed  8 Apr - 18:15  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  head -n 1000 unaligned__ptprc_cell_barcodes.sam > unaligned__ptprc_cell_barcodes_first1000.sam
(samtools)
 Wed  8 Apr - 18:17  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  samtools -Sb unaligned__ptprc_cell_barcodes_first1000.bam
[main] unrecognized command '-Sb'
(samtools)
 ✘  Wed  8 Apr - 18:17  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  samtools view -Sb unaligned__ptprc_cell_barcodes_first1000.bam
[E::hts_open_format] Failed to open file unaligned__ptprc_cell_barcodes_first1000.bam
samtools view: failed to open "unaligned__ptprc_cell_barcodes_first1000.bam" for reading: No such file or directory
(samtools)
 ✘  Wed  8 Apr - 18:17  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  samtools view -Sb unaligned__ptprc_cell_barcodes_first1000.sam > unaligned__ptprc_cell_barcodes_first1000.bam
(samtools)
 Wed  8 Apr - 18:17  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  samtools view -Sb unaligned__ptprc_cell_barcodes_first1000.sam > unaligned__ptprc_cell_barcodes_first1000__no_header.bam
(samtools)
 Wed  8 Apr - 18:17  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  samtools reheader header.sam unaligned__ptprc_cell_barcodes_first1000__no_header.bam > unaligned__ptprc_cell_barcodes_first1000.bam
```



### Copy files over and move into correct directories

#### Copy files

```
(samtools)
 Wed  8 Apr - 11:20  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz 1☀ 1● 
 olga@lrrr  cp /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b/mouse* .
```


#### Look at directory structure of existing mouse_lung folder

```
(samtools)
 Wed  8 Apr - 11:23  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz 1☀ 1● 
 olga@lrrr  tree mouse_lung
mouse_lung
└── outs
    ├── filtered_gene_bc_matrices
    │   └── placeholder
    │       └── barcodes.tsv
    ├── possorted_genome_bam.bam
    └── possorted_genome_bam.bam.bai
```

#### Copy files into their place for the tgz file


```
(samtools)
 ✘  Wed  8 Apr - 11:24  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz 1☀ 1● 
 olga@lrrr  mkdir -p mouse_brown_fat_ptprc_plus_unaligned/outs/filtered_gene_bc_matrices/placeholder/       (samtools)
 Wed  8 Apr - 11:24  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz 1☀ 1● 
 olga@lrrr  cp mouse_brown_fat_ptprc_plus_unaligned.bam mouse_brown_fat_ptprc_plus_unaligned/outs/possorted_genome_bam.bam
(samtools)
 Wed  8 Apr - 11:24  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz 2☀ 1● 
 olga@lrrr  cp mouse_brown_fat_ptprc_plus_unaligned.bam.bai mouse_brown_fat_ptprc_plus_unaligned/outs/possorted_genome_bam.bam.bai
(samtools)
 Wed  8 Apr - 11:24  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz 2☀ 1● 
 olga@lrrr  cp mouse_brown_fat_ptprc_plus_unaligned__barcodes.tsv mouse_brown_fat_ptprc_plus_unaligned/outs/filtered_gene_bc_matrices/placeholder/
```

#### Compress into tar.gz file

```
(samtools)
 Wed  8 Apr - 11:24  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz 2☀ 1● 
 olga@lrrr  tar -zcvf mouse_brown_fat_ptprc_plus_unaligned.tgz mouse_brown_fat_ptprc_plus_unaligned
mouse_brown_fat_ptprc_plus_unaligned/
mouse_brown_fat_ptprc_plus_unaligned/outs/
mouse_brown_fat_ptprc_plus_unaligned/outs/possorted_genome_bam.bam.bai
mouse_brown_fat_ptprc_plus_unaligned/outs/filtered_gene_bc_matrices/
mouse_brown_fat_ptprc_plus_unaligned/outs/filtered_gene_bc_matrices/placeholder/
mouse_brown_fat_ptprc_plus_unaligned/outs/filtered_gene_bc_matrices/placeholder/mouse_brown_fat_ptprc_plus_unaligned__barcodes.tsv
mouse_brown_fat_ptprc_plus_unaligned/outs/possorted_genome_bam.bam
```

#### Remove files from flat directory

```
(samtools)
 Wed  8 Apr - 11:31  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz ✔ 
 olga@lrrr  rm -rf mouse_brown_fat_ptprc_plus_unaligned.bam*                                                (samtools)
 Wed  8 Apr - 11:31  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz 2‒ 
 olga@lrrr  git rm -rf mouse_brown_fat_ptprc_plus_unaligned.bam*
zsh: no matches found: mouse_brown_fat_ptprc_plus_unaligned.bam*
(samtools)
 ✘  Wed  8 Apr - 11:31  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz 2‒ 
 olga@lrrr  git rm -rf mouse_brown_fat_ptprc_plus_unaligned.bam
rm 'testdata/mouse_brown_fat_ptprc_plus_unaligned.bam'
(samtools)
 Wed  8 Apr - 11:31  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz 1‒1± ⚑ 
 olga@lrrr  git rm -rf mouse_brown_fat_ptprc_plus_unaligned.bam.bai
rm 'testdata/mouse_brown_fat_ptprc_plus_unaligned.bam.bai'
(samtools)
 Wed  8 Apr - 11:31  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz ‒2± ⚑ 
 olga@lrrr  rm mouse_brown_fat_ptprc_plus_unaligned__barcodes.tsvrm: remove regular file 'mouse_brown_fat_ptprc_plus_unaligned__barcodes.tsv'? y
````

#### Move barcodes file to *correct* place

```
(samtools)
 Wed  8 Apr - 11:36  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz 1● 1‒ 
 olga@lrrr  mv mouse_brown_fat_ptprc_plus_unaligned/outs/filtered_gene_bc_matrices/placeholder/mouse_brown_fat_ptprc_plus_unaligned__barcodes.tsv mouse_brown_fat_ptprc_plus_unaligned/outs/filtered_gene_bc_matrices/placeholder/barcodes.tsv
```


#### Re-add smaller bam and bai

```
(samtools)
 Wed  8 Apr - 18:19  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  cd ~/code/nf-core/test-datasets--kmermaid/testdata/mouse_brown_fat_ptprc_plus_unaligned/outs/(samtools)
 Wed  8 Apr - 18:31  ~/code/nf-core/test-datasets--kmermaid/testdata/mouse_brown_fat_ptprc_plus_unaligned/outs   olgabot/kmermaid-unaligned-tgz-v2 1● 
 olga@lrrr  ll
Permissions Size User Group Date Modified Git Name
drwxr-xr-x     - olga czb    8 Apr 11:24   -- filtered_gene_bc_matrices
.rw-r--r--   38M olga czb    8 Apr 11:24   -- possorted_genome_bam.bam
.rw-r--r--   68k olga czb    8 Apr 11:24   -- possorted_genome_bam.bam.bai
(samtools)
 Wed  8 Apr - 18:31  ~/code/nf-core/test-datasets--kmermaid/testdata/mouse_brown_fat_ptprc_plus_unaligned/outs   olgabot/kmermaid-unaligned-tgz-v2 1● 
 olga@lrrr  cp  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b/mouse_brown_fat_ptprc_plus_unaligned_first1000.bam possorted_genome_bam.bam
cp: overwrite 'possorted_genome_bam.bam'? y
(samtools)
 Wed  8 Apr - 18:31  ~/code/nf-core/test-datasets--kmermaid/testdata/mouse_brown_fat_ptprc_plus_unaligned/outs   olgabot/kmermaid-unaligned-tgz-v2 2● 
 olga@lrrr  cp  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b/mouse_brown_fat_ptprc_plus_unaligned_first1000.bam.bai possorted_genome_bam.bam.bai
cp: overwrite 'possorted_genome_bam.bam.bai'? y
```

#### Re-make targz

```
(samtools)
 Wed  8 Apr - 18:32  ~/code/nf-core/test-datasets--kmermaid/testdata   olgabot/kmermaid-unaligned-tgz-v2 3● 
 olga@lrrr  tar -zcvf mouse_brown_fat_ptprc_plus_unaligned.tgz mouse_brown_fat_ptprc_plus_unaligned
mouse_brown_fat_ptprc_plus_unaligned/
mouse_brown_fat_ptprc_plus_unaligned/outs/
mouse_brown_fat_ptprc_plus_unaligned/outs/possorted_genome_bam.bam.bai
mouse_brown_fat_ptprc_plus_unaligned/outs/filtered_gene_bc_matrices/
mouse_brown_fat_ptprc_plus_unaligned/outs/filtered_gene_bc_matrices/placeholder/
mouse_brown_fat_ptprc_plus_unaligned/outs/filtered_gene_bc_matrices/placeholder/barcodes.tsv
mouse_brown_fat_ptprc_plus_unaligned/outs/possorted_genome_bam.bam
```