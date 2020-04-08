
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



```
 ✘  Wed  8 Apr - 06:45  /mnt/data_sm/olga/nextflow-intermediates-lrrr/7a/ed752d5e82cdd39e5fde9c21c3423b 
 olga@lrrr  time samtools view -f4 MACA_18m_M_BAT_52__possorted_genome_bam.bam | fgrep --file top_barcodes__as_tenx_barcode.txt | samtools view -Sbh > MACA_18m_M_BAT_52__possorted_genome_bam__ptprc_top_barcodes_unaligned.bam
samtools view -f4 MACA_18m_M_BAT_52__possorted_genome_bam.bam  161.56s user 5.34s system 99% cpu 2:47.73 total
fgrep --file top_barcodes__as_tenx_barcode.txt  14.12s user 2.02s system 9% cpu 2:47.73 total
samtools view -Sbh >   0.00s user 0.01s system 0% cpu 2:47.73 total
```

