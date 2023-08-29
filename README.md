# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

- For nf-core/taxprofiler CI test information see [here](#taxprofiler-ci-test-specific-information)
- For nf-core/taxprofiler AWS full test information see [here](#taxprofiler-aws-full-test-specific-information)

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

1.  [Add a new test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
2.  [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Downloading test data

Due the large number of large files in this repository for each pipeline, we highly recommend cloning only the branches you would use.

```bash
git clone <url> --single-branch --branch <pipeline/modules/branch_name>
```

To subsequently clone other branches[^1]

```bash
git remote set-branches --add origin [remote-branch]
git fetch
```

## Taxprofiler CI test specific information

### FASTQ

The main CI test data used for nf-core/taxprofiler is from [Maixner et al. (2021) _Curr. Bio._](https://doi.org/10.1016/j.cub.2021.09.031), with ENA project accession ID: PRJEB44507. The following selected libraries were all sequenced on an Illumina MiSeq, and were selected due to their small size (~1million reads, <100MB) and known mixture of (gut) bacteria, (ancient human) eukaryotes, and (yeast) fungi (according to the results of the paper).

- ERX5474937
- ERX5474932
- ERX5474930
- ERX5474936

Data was downloaded with nf-core/fetchNGS 1.5 (with Nextflow 21.10.06):

```bash
nextflow run nf-core/fetchngs --input maixner2021_acc_codes.txt --input_type sra
```

FASTQ files are stored under `data/fastq/`

Test data for long reads with ENA project accession ID: PRJEB29152. They were subsampled with seqtk 1.3-r106:

```bash
seqtk sample ERR3201952.fastq.gz 10000 > ERR3201952.fastq.gz
```

### FASTA

One of the files was converted to FASTA file with seqtk 1.3-r106

```bash
seqtk seq -a  ERX5474930_ERR5766174_1.fastq.gz > ERX5474930_ERR5766174_1.fa.gz
```

FASTA files are stored under `data/fasta/`

### Databases

An abundant species found in the above dataset is _Penicillium roqueforti_.
The genome and translations of P. roqueforti were downloaded with:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/533/775/GCF_015533775.1_ASM1553377v1/GCF_015533775.1_ASM1553377v1_genomic.fna.gz # P. roqueforti
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/533/775/GCF_015533775.1_ASM1553377v1/GCF_015533775.1_ASM1553377v1_translated_cds.faa.gz # P. roqueforti

```

In addition we include the human mitochondrial genome.

```bash
curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=251831106&extrafeat=null&conwithfeat=on&hide-cdd=on'| gzip > NC_012920.1.fa.gz # H. sapiens mito
curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta_cds_aa&id=251831106&extrafeat=null&conwithfeat=on&hide-cdd=on'| gzip > NC_012920.1.faa.gz # H. sapiens mito
```

All four were unzipped

```bash
gunzip *.gz
```

Note all test database files should be placed in an archive containing the following structure

```
test-db-<tool>.tar.gz
└── test-db-<tool>/
     ├── <file1>
     └── <file2>
```

#### MALT

MALT Version 0.6.1

```bash
wget https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Feb2022.db.zip
unzip megan-nucl-Feb2022.db
malt-build -i *.{fna,fa} -s DNA -d taxprofiler-testdb -t 8 -st 4 -a2t megan-nucl-Feb2022.db
```

#### Kraken2

Kraken Version 2.1.2

```bash
kraken2-build --download-taxonomy --db taxprofiler-testdb
kraken2-build --add-to-library ../raw/GCF_000146045.2_R64_genomic.fna --db taxprofiler-testdb/
kraken2-build --add-to-library ../raw/NC_012920.1.fa --db taxprofiler-testdb/
kraken2-build --build --db taxprofiler-testdb/
kraken2-build --clean --db taxprofiler-testdb/
```

#### Bracken

Bracken Version 2.7

```bash
nextflow run https://github.com/Midnighter/kraken2-bracken-test-db -profile docker
```

#### Centrifuge

Centrifuge version 1.0.4

Downloaded taxonomy files

```
centrifuge-download -o taxonomy taxonomy
```

Made custom seqid2taxid.map

```
NC_001133.9    4392
NC_012920.1    9606
NC_001134.8    4392
NC_001135.5    4392
NC_001136.10    4392
NC_001137.3    4392
NC_001138.5    4392
NC_001139.9    4392
NC_001140.6    4392
NC_001141.2    4392
NC_001142.9    4392
NC_001143.9    4392
NC_001144.5    4392
NC_001145.3    4392
NC_001146.8    4392
NC_001147.6    4392
NC_001148.4    4392
NC_001224.1    4392
```

Combined the two FASTAs together

```
cat *.{fa,fna} > input-sequences.fna
```

Then build the CF database files

```bash
centrifuge-build -p 4 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp input-sequences.fna taxprofiler_cf
```

#### DIAMOND

Diamond Version 2.0.15

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip

## warning: large file!
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz

## warning: takes a long time!
cat ../raw/*.faa | diamond makedb -d testdb-diamond --taxonmap prot.accession2taxid.FULL.gz --taxonnodes nodes.dmp --taxonnames names.dmp

rm *dmp *txt *gz *prt *zip
```

#### KrakenUniq

KrakenUniq version 1.0.0

> ⚠️  This database _does not_ use the specified files used in the other databases, as this built into a database that was tool large.

This database includes the SARS-CoV2 genome used on the nf-core/modules test-datasets repository (NCBI Accession: MT192765.1).

It was generated using the nf-core/module KRAKENUNIQ_BUILD module.

#### ganon

ganon version 1.5.1

```bash
ganon build-custom --threads 4 --input *.fa --db-prefix test-db-ganon --verbose -x ncbi --write-info-file --ncbi-sequence-info --ncbi-file-info -e fa --input-target sequence
```

#### kmcp

kmcp version 0.9.1

```bash
mkdir gtdb-genomes

## Copy the downloaded fasta files for Penicillium roqueforti and Human genome mitochondral to folder gtdb-genomes.

## Rename the file for Penicillium roqueforti to match the seqid2taxid.map file
mv GCF_015533775.1_ASM1553377v1_genomic.fna.gz NW_024067565.1.fna.gz 

kmcp compute -k 21 -n 10 -l 150 -O tmp-k21-n10-l150 -I gtdb-genomes
kmcp index -f 0.3 -n 1 -j 32 -I tmp-k21-n10-l150 -O gtdb.kmcp 
```


## Taxprofiler AWS Full Test specific-information

### FASTQ

The main AWS full test data used for nf-core/taxprofiler is from [Meslier et al. (2022) _Sci. Data_](https://doi.org/10.1038/s41597-022-01762-z), with ENA project accession ID: PRJEB52977. The following selected libraries were all sequenced on an Illumina HiSeq 3000 and ONT Minion R9.

They were selected as a benchmarking dataset containing a semi-complex microbial community with strains that have known reference genomes, and with multiple sequencing runs one of sample of the Illumina dataset. ENA Experiment IDs are as follows

- ONT MiniION R9
    - ERX9314125
    - ERX9314126
    - ERR9765782
- Illumina HiSeq 3000
    - ERX9314116
    - ERX9314117
    - ERX9314118 (x2 runs)

FASTQ files for the `samplesheet_full.tsv` are stored on the [EBI ENA servers](https://www.ebi.ac.uk/ena/browser/view/PRJEB52977)

### FASTA

FASTA files for use in database construction were identified based on Supp. Table 1 from [Meslier (2022)](https://doi.org/10.1038/s41597-022-01762-z), which were copy-pasted into an empty text file called `meslier2022_supptab1.tsv`.

We then used the AWK combined with the NCBI Datasets package (v14.7.0) to download the reference genomes and protein translations of each strain from the file.

```bash
awk -F'\t' '{print $2}' meslier2022_supptab1.tsv | xargs -I '{}' datasets download genome accession {} --include genome,protein --filename meslier2022_fasta/{}.zip
```

In this case I had one failure for assembly `GCA_000009225.1`, which has since been replaced with `GCA_931907645.1`, which was downloaded manually using the command above. A further five accessions have been suppressed with no replacement, and thus these were not included in the databases.

Once downloaded, we need to unpack the `datasets` and rename the protein translation for each file to make them unique.

```bash
for i in meslier2022_fasta/*zip; do
    f_basename=$(basename $i)
    unzip $i ncbi_dataset/data/*/*.f* -d meslier2022_fasta/
    mv meslier2022_fasta/ncbi_dataset/data/${f_basename%%.zip}/protein.faa meslier2022_fasta/ncbi_dataset/data/${f_basename%%.zip}/${f_basename%%.zip}.faa
done
```

### Databases

The following sections describe how each database within `databases_full.csv` were constructed, following on from the full test FASTA download in the section above.

> ⚠️ Be aware all commands include thread/CPU and/or memory parameters which may not apply to your machine.

#### Kraken2

The following steps were performed using Kraken2 v2.1.2

First create a working directory

```bash
mkdir -p meslier2022/kraken2
```

Download NCBI taxonomy files

```bash
kraken2-build --download-taxonomy --db meslier2022/kraken2
```

Index all the FASTA files to the Kraken2 database

```bash
find meslier2022_fasta/ -name '*.fna' -print0 | xargs -0 -I{} -P 8 -n1 kraken2-build --add-to-library {} --db meslier2022/kraken2
```

Build the database

```bash
kraken2-build --build --threads 32 --db /raven/ptmp/jfellowsy/databases/taxprofiler_full_test/meslier2022/kraken2
```

Copy the `seqid2taxid.map` for reuse in other profiler database construction

```bash
cp meslier2022/kraken2/seqid2taxid.map .
```

#### Bracken

The following steps were performed using Bracken (v2.8).

Bracken requires an existing Kraken2 database to build upon, for which we use the database built in the previous section.

We also need to specify the average read length, which according to the ENA was sequenced paired-end 150, so we will select 150 bp.

Make a working directory

```bash
mkdir -p meslier2022/bracken
```

Make a symbolic link from the Kraken2 database

```bash
cd meslier2022/bracken
ln -s ../kraken2/* .
cd ../../
```

Build the database

```bash
bracken-build -t 72 -l 150 -d meslier2022/bracken
```

#### KrakenUniq

The following steps were performed using KrakenUniq (v1.0.3).

Make a working directory

```bash
mkdir -p meslier2022/krakenuniq/library
```

Download NCBI taxonomy files

```bash
krakenuniq-download --db meslier2022/krakenuniq/ taxonomy
```

Make the FASTA file directory following structure required for the build command, and symlink in the FASTAs

```
mkdir -p meslier2022/krakenuniq/library
cd meslier2022/krakenuniq/library
ln -s ../../../meslier2022_fasta/ncbi_dataset/data/*/*.fna .
cd ../../../
```

Copy the saved Kraken2 `seqid2taxid.map` file into the corresponding KrakenUniq location

```bash
cp seqid2taxid.map meslier2022/krakenuniq/
```

Run the build command. Note due to a suspected misconfiguration in the `conda` package used, we had to manually call the `jellyfish` `bin/` path

```bash
krakenuniq-build --build --threads 72 --db meslier2022/krakenuniq --work-on-disk --jellyfish-bin $(which jellyfish)
```

#### MALT

The following steps were performed using MALT (v0.6.1).

Make a working directory

```bash
mkdir -p meslier2022/malt
```

Download and unpack the required taxonomy files

```bash
wget https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Feb2022.db.zip
unzip megan-nucl-Feb2022.db
```

Run the build command

```bash
malt-build -J-Xmx490G -i meslier2022_fasta/ncbi_dataset/data/*/*.fna -s DNA -d meslier2022/malt -t 72 -st 16 -a2t megan-nucl-Feb2022.db
```

#### DIAMOND

The following steps were performed used DIAMOND (v2.0.15).

Make a working directory

```bash
mkdir -p meslier2022/diamond
```

Download and unpack the required taxonomy files.

> ⚠️ The accession2taxid file is very large!

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
```

Run the build command

```bash
cat fasta/*.faa | diamond makedb --threads 72 -v --log -d meslier2022/diamond/diamond --taxonmap prot.accession2taxid.FULL.gz --taxonnodes nodes.dmp --taxonnames names.dmp
```

#### Centrifuge

The following steps were performed used Centrifuge (v1.0.4).

Make a working directory

```bash
mkdir -p meslier2022/centrifuge
```

Download and unpack the required taxonomy files

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip new_taxdump.zip
```

Combine all reference sequences into a single FASTA

```bash
cat meslier2022_fasta/ncbi_dataset/data/*/*.fna > centrifuge_sequences.fna
```

Run the build command using the saved Kraken2 `seqid2taxid.map` file from the corresponding section above

```bash
centrifuge-build -p 32 --conversion-table seqid2taxid.map --taxonomy-tree nodes.dmp --name-table names.dmp centrifuge_sequences.fna meslier2022/centrifuge/centrifuge
```

#### Kaiju

The following steps were performed used Kaiju (v1.9.2).

Make a working directory

```bash
mkdir -p meslier2022/kaiju
```

Download and unpack the required taxonomy files

> ℹ️ You can skip this step if you also download the same files already with Centrifuge

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip new_taxdump.zip
```

Combine all protein translation sequences into a single amino acid FASTA

```bash
cat meslier2022_fasta/ncbi_dataset/data/*/*.faa > kaiju_sequences.faa
```

We now must make modifications to the combined amino acid FASTA file to ensure the headers are formatted expected by Kaiju.

The headers are expected to be the numeric NCBI taxon identifiers of the protein sequences, which can optionally be prefixed by another identifier (e.g. a counter) followed by an underscore.

First we extract all the Kaiju headers

```bash
cut -d ' ' -f1 kaiju_sequences.faa > kaiju_sequences_accession.faa
grep ">" kaiju_sequences_accession.faa > kaiju_sequences_accession.txt
sed -i 's/>//g' kaiju_sequences_accession.txt
```

We then use the following R (v4.2.2) commands with the taxonomizr package (0.7.1) to pull the NCBI Taxonomic ID from the NCBI Protein accession IDs in the FASTA file.

We firstly install the taxonmizer package

```r
install.packages("taxonomizr")
library(taxonomizr)
```

Download the necessary nodes and names files from NCBI

```r
getNamesAndNodes()
```

Download and load the prot file from NCBI. This may take some time as it is a large file

```r
getAccession2taxid(types='prot')
read.names.sql('names.dmp','accessionTaxa.sql')
read.nodes.sql('nodes.dmp','accessionTaxa.sql')
read.accession2taxid(list.files('.','prot.accession2taxid.gz'),'accessionTaxa.sql')
```

Load the FASTA headers, bind with taxonomy ID information and save the file.

```r
kaiju_sequences_taxid <- read.csv("kaiju_sequences_accession.txt",header=FALSE)
taxaId < -accessionToTaxa(kaiju_sequences_taxid$V1,"accessionTaxa.sql")
write.csv(taxaId,file="kaiju_sequences_taxId.txt",row.names=F)
```

Exit the R session

```r
quit()
```

We can then create a new FASTA with the correct numeric IDs using AWK

```bash
paste kaiju_sequences_accession.txt kaiju_sequences_taxId.txt >  kaiju_accession_taxid.txt
awk 'FNR==NR {f2[$1]=$2;next} /^>/ { for (i in f2) { if (index(substr($1,2), i)) { print ">"f2[i]; next } } }1' kaiju_accession_taxid.txt kaiju_sequences_accession.faa > kaiju_sequences.faa
```

We can copy the previously downloaded NCBI taxonomy files into the Kaiju working directory


```bash
cp nodes.dmp names.dmp meslier2022/kaiju/
```

And build the database

```bash
kaiju-mkbwt -n 32 -a ACDEFGHIKLMNPQRSTVWY -o meslier2022/kaiju/kaiju kaiju_sequences.faa
kaiju-mkfmi meslier2022/kaiju/kaiju
```

#### MetaPhlAn3

The following steps were performed used MetaPhlAn3 (v3.1.0).

```bash
metaphlan --install --bowtie2db meslier2022/metaphlan3/
```

#### mOTUs

The following steps were performed used MetaPhlAn3 (v3.0.3).

```bash
motus downloadDB
```

You then will need to find the location of the downloaded file in the mOTUs installation directory, e.g.,

```console
<conda_installion>/envs/motus/lib/python3.9/site-packages/motus/db_mOTU/
```

#### ganon

Create directory, then run command. Note taxonomy files will be automatically downloaded for you.

```
mkdir -p meslier2022/ganon

ganon build-custom --threads 4 --input meslier2022_fasta/ncbi_dataset/data/GCA_*/*.fna --db-prefix meslier2022/ganon/ganon --verbose -x ncbi --write-info-file --ncbi-sequence-info --ncbi-file-info -e fa --input-target sequence
```

#### kmcp

Create directory, then copy the input files to the directory and rename them.

```
mkdir -p meslier2022/kmcp
cd kmcp
cp meslier2022_fasta/ncbi_dataset/data/GCA_*/*.fna .
```

Open an empty file `rename_files.sh` and add the following piece of code:

```bash
#!/bin/bash
rename_fna_files_in_directory() {
    local dir="$1"

    # Iterate through .fna files in the directory
    for file in "$dir"/*/*_genomic.fna; do
        if [[ -f "$file" ]]; then
            new_name="${file%_ASM*}.fna"
            mv "$file" "$new_name"
        fi
    done
}

# Iterate through subdirectories
for folder in */; do
    if [[ -d "$folder" ]]; then
        rename_fna_files_in_directory "$folder"
    fi
done
```

```
chmod +x rename_files.sh
./rename_files.sh
```

After the files have been renamed, gather the filenames in a text file. You will need those in order to run the `kmcp compute`

```bash
ls *.fna > list.txt

## The names.dmp, nodes.dmp and seqid2taxid.map are needed for kmcp profile step
cp meslier2022/kraken2/seqid2taxid.map .

mkdir taxdump

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip new_taxdump.zip

cp names.dmp taxdump/
cp nodes.dmp taxdump/

kmcp compute -k 21 -n 10 -l 150 -O tmp-k21-210-l150 -i list.txt
kmcp index -I tmp-k21-210-l150/ --threads 8 --num-hash 1 --false-positive-rate 0.3 --out-dir refs.kmcp
```

## Database Archive Creation

To make the compressed TAR, we must make sure all symlinks are followed as necessary. It is recommended to run the cleanup commands below _prior_ to archiving, however it is critical that Bracken archiving is performed BEFORE running the Kraken2 cleanup.

```bash
tar -hzcvf <toolname>.tar.gz meslier2022_<toolname>/
```

These archives are then used in the `database_full.csv` file.

### Cleanup

For KrakenUniq

```bash
rm meslier2022/krakenuniq/*.{log,counts,tsv,map,jdb,txt} meslier2022/krakenuniq/database0.kbd meslier2022/krakenuniq/library/ meslier2022/krakenuniq/taxonomy
```

For Kraken2

> ⚠️ Only do once Bracken/KrakenUniq databases are built, if they are required

```bash
kraken2-build --db meslier2022/kraken2/ --clean
```

For MALT

```bash
rm *.db *.db.zip
```

For DIAMOND

```bash
rm dmp readme.txt *prt taxdmp.zip  prot.accession2taxid.FULL.gz
```

Centrifuge

```bash
new_taxdump.zip
```

KAIJU

```bash
rm meslier2022/kaiju/*.{bwt,sa}
```

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
