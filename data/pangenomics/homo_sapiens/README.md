# Test Data

## Table of contents

- [Data Fetching](#data-fetching)
- [Data Generation](#data-generation)

## Data Fetching

The initial FASTA file can be optained from the [PGGB](https://github.com/pangenome/pggb) repository:

```sh
git clone https://github.com/pangenome/pggb
cd pggb
```

The FASTA is at `data/HLA/V-352962.fa.gz`.

## Data Generation

We use `PGGB` to generate all the output files from the initial FASTA:

```sh
./pggb -i ~data/HLA/V-352962.fa.gz -p 70 -s 3000 -G 20000 -n 10 -t 1 -v -U -V 'gi|568815592:#' -o nf-core -M -C cons,100,1000,10000 -L
```

Let's take a look at the output:

```sh
tree nf-core/
nf-core/
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.12-10-2021_13:30:39.log
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.12-10-2021_13:30:39.params.yml
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.cons@0__y_0_1000000.gfa
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.cons@10000__y_0_1000000.gfa
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.cons@1000__y_0_1000000.gfa
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.cons@100__y_0_1000000.gfa
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.fix.affixes.tsv.gz
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.fix.gfa
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.gfa
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.gi_568815592.vcf
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.og
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.og.lay
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.og.lay.draw_mqc.png
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.og.lay.draw.png
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.og.viz_depth.png
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.og.viz_inv.png
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.og.viz_mqc.png
├── V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.og.viz_slide.png
├── V-352962.fa.gz.7f683ec.2ff309f.seqwish.gfa
├── V-352962.fa.gz.7f683ec.2ff309f.seqwish.gfa.prep.gfa
└── V-352962.fa.gz.7f683ec.wfmash.paf

0 directories, 21 files
```

We create a new folder for the test data and copy our required files there:

```sh
mkdir pangenome-test-data
cp data/HLA/V-352962.fa.gz pangenome-test-data/pangenome.fa.gz
cp nf-core/V-352962.fa.gz.7f683ec.wfmash.paf pangenome-test-data/pangenome.paf
cp nf-core/V-352962.fa.gz.7f683ec.2ff309f.seqwish.gfa pangenome-test-data/pangenome.seqwish.gfa
cp nf-core/V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.gfa pangenome-test-data/pangenome.smoothxg.gfa
cp nf-core/V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.fix.gfa pangenome-test-data/pangenome.gfaffix.gfa
cp nf-core/V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.og pangenome-test-data/pangenome.og
cp nf-core/V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.og.lay pangenome-test-data/pangenome.og.lay
gunzip -c pangenome-test-data/pangenome.fa.gz > pangenome-test-data/pangenome.fa
gzip -k pangenome-test-data/pangenome.paf
cp nf-core/V-352962.fa.gz.7f683ec.2ff309f.9799452.smooth.gi_568815592.vcf pangenome-test-data/pangenome.vcf
bgzip pangenome-test-data/pangenome.vcf
```

That's it :)