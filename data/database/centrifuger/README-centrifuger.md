#Centrifuger Test Database

##Introduction 
This file shows how the DB files used for centrifuger module were generated in order to incorporate them in the taxprofiler pipeline. 
It follows the steps from https://github.com/nf-core/test-datasets/tree/taxprofiler#databases.

##Download Database
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/533/775/GCF_015533775.1_ASM1553377v1/GCF_015533775.1_ASM1553377v1_genomic.fna.gz # P. roqueforti
```

Include human mitochondrial genome
```bash 
curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=251831106&extrafeat=null&conwithfeat=on&hide-cdd=on'| gzip > NC_012920.1.fa.gz # H. sapiens mito
```

Unzip the files:
```bash
gunzip *.gz
```
##Combine fastas together 
```bash
cat *.{fa,fna} > input-sequences.fna
```

##Download taxonomy files
```bash
centrifuger-download -o taxonomy taxonomy
```

## Create custom seqid2taxid.map
Map P.requefort scafolds to taxID 5082
```bash
grep "^>" GCF_015533775.1*.fna | awk '{print $1"\t5082"}' | sed 's/>//' > seqid2taxid.map
```

Add the human mito to taxid 9606
```bash
echo -e "NC_012920.1\t9606" >> seqid2taxid.map
```

## Build index files  for centrifuger with centrifuger-build - using Docker to avoid platform-specifc seqfaults

```bash
docker run --rm -v $(pwd):/data -w /data quay.io/biocontainers/centrifuger:1.1.0--hf426362_0 centrifuger-build \
 -r input-sequences.fna \
--conversion-table seqid2taxid.map \
--taxonomy-tree taxonomy/nodes.dmp \
--name-table taxonomy/names.dmp \
-t 4 \
-o test-db-centrifuger/test-db-centrifuger
```


##Create a tarball
```bash
tar czvf test-db-centrifuger.tar.gz test-db-centrifuger
```


