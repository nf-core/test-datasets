# test-datasets `lncpipe`

nf-core is a collection of high quality Nextflow pipelines.

## Content of this repository

NOTE: A full packaged test data set can be downloand from http://cancerbio.info/pub/lncpipe/testdata.tar.gz

This branch contains test data instruction for the nf-core/lncpipe pipeline.   

## Run test data of lncpipe   

As lncpipe involved integrated analysis of multi samples, we provided a four-sample packaged test data that have two different  experiment conditions. To run the test of lncpipe, plz run the following command step by step :  

```shell
      #prepare test data 
      wget http://cancerbio.info/pub/lncpipe/testdata.tar.gz
      tar -xvzf testdata.tar.gz 
      cd testdata
      #run test command 
      nextflow run nf-core/lncpipe -profile test,docker
```

## Test data content 


```
      ├── design.file
      ├── Fastq
      │   ├── dPDLSCs1_1.fastq.gz
      │   ├── dPDLSCs1_2.fastq.gz
      │   ├── dPDLSCs2_1.fastq.gz
      │   ├── dPDLSCs2_2.fastq.gz
      │   ├── uPDLSCs1_1.fastq.gz
      │   ├── uPDLSCs1_2.fastq.gz
      │   ├── uPDLSCs2_1.fastq.gz
      │   └── uPDLSCs2_2.fastq.gz
      └── Genome
          ├── chr22.fa
          ├── gencode.chr22.gtf
          ├── hisat_index
          │   ├── chr22.1.ht2
          │   ├── chr22.2.ht2
          │   ├── chr22.3.ht2
          │   ├── chr22.4.ht2
          │   ├── chr22.5.ht2
          │   ├── chr22.6.ht2
          │   ├── chr22.7.ht2
          │   └── chr22.8.ht2
          └── lncipedia.chr22.gtf
```  
      
* `design.file` stored the test exprimental design for performing comparision   

```
      D:dPDLSCs1,dPDLSCs2
      u:uPDLSCs1,uPDLSCs2
```  
* `Fastq` folder contains paired, compressed fastq files, also known as raw reads.   

* `Genome` folder contains several files explained below:  

     * `chr22.fa` genome reference of chromosome 22 with fasta format.  
      
     * `hisat_index` folder contains the hisat2 index file build from chr22.fa
      
     * `lncipedia.chr22.gtf` GTF files grepped from [lncpedia_4.0.gtf](https://lncipedia.org/downloads/lncipedia_5_0_hc_hg38.gtf)  
     
     * `gencode.chr22.gtf` GTP files grepped from [GENCODE](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz)
     
