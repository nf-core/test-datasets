# nf-core/mitotags test files

This is the minimal test files to test the nf-core/mitotags package.

## input 

ChrM_testData_R1_001.fastq.gz

```
bedtools bamtofastq -i ChrM_subset.bam -fq ChrM_testData_R1_001.fastq -fq2 ChrM_testData_R2_001.fastq
gzip *.fastq
```

## bam 

ChrM_subset.bam


```
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_atac_possorted_bam.bam
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_atac_possorted_bam.bam.bai
samtools view 10k_PBMC_Multiome_nextgem_Chromium_X_atac_possorted_bam.bam chrM:15000-22000 -b > Chrm_subset.bam
## find the 300 highest chrM read count barcodes and split on them.
gunzip barcodes.tsv
samtools view -h Chrm_subset.bam | head -n 3000 | grep "^@" > ChrM_subset.sam
samtools view Chrm_subset.bam | grep -f barcodes.tsv  >> ChrM_subset.sam
samtools view ChrM_subset.sam -b > ChrM_subset.bam
gzip barcodes.tsv
rm Chrm_subset.bam 
rm ChrM_subset.sam

```

## bai

ChrM_subset.bam.bai

```
samtools index ChrM_subset.bam
```

## barcodes  

barcodes.tsv.gz

```
samtools view Chrm_subset.bam | grep -o "CB:Z:[AGCT]*-1" > barcodes_total.txt
```

R script to get the top 300 barcodes:

```
dat = scan( 'barcodes_total.txt', what=character())
write ( names( sort(table(dat), decreasing=TRUE)[1:300]), file="barcodes.tsv" )
```

Or on the command line (slow!):

```
sort barcodes_total.txt | uniq -c > barcodes_unique.txt
sort -k 2n  barcodes_unique.txt| head -n 300 | cut -f1 > barcodes.tsv
```

Clean up

```
gzip barcodes.tsv
rm barcodes_total.txt
```

## genome 
ChrM.fa.gz
