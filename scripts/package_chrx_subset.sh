#!/bin/bash
# Package chrX_subset/ from a pipeline run on the subset. Usage: package.sh subset_dir run_dir bam_dir out_dir
set -euo pipefail
V=$1; RUN=$2; B=$3; P=$4
R=$RUN/results; W=$RUN/work
task() { (cd $W; grep -l "### name: '.*$1" */*/.command.run | head -1 | xargs dirname); }
rm -rf $P; mkdir -p $P/reference $P/testdata/salmon_quant $P/testdata/drimseq $P/testdata/dexseq $P/samplesheet
cp $V/reference/X.fa.gz $V/reference/genes_chrX.gtf $P/reference/
mi=$(task MISOPY_INDEX); cp -L $W/$mi/genes_chrX_genes.gff3 $P/reference/genes_chrX.gff3
tmp=$(mktemp -d); cp -rL $W/$mi/index $tmp/misopy_index; (cd $tmp && tar czf $P/reference/misopy_index.tar.gz misopy_index); rm -rf $tmp
gzip -9 -c $W/$(task DEXSEQ_PREPAREANNOTATION)/DEXSeq.gff > $P/testdata/dexseq/genes_chrX.DEXSeq.gff.gz
for s in ERR188383 ERR188428 ERR188454 ERR204916; do
    cp $V/testdata/${s}_chrX_{1,2}.fastq.gz $P/testdata/
    # $B holds the BAMs cut down to 4,000 pairs with a spliced mate (samtools view -N)
    cp $B/$s.Aligned.out.bam $B/$s.Aligned.out.bam.bai $B/$s.Aligned.toTranscriptome.out.bam $P/testdata/
    cp -rL $R/salmon/$s $P/testdata/$s; (cd $P/testdata/$s && tar czf ../$s.tar.gz .)
done
tmp=$(mktemp -d); for s in ERR188383 ERR188428 ERR188454 ERR204916; do cp -rL $R/salmon/$s $tmp/; done
cp -L $W/$(task GFFREAD_TX2GENE)/genes_chrX.tx2gene.tsv $tmp/
(cd $tmp && tar czf $P/testdata/salmon_quant/salmon_quant.tar.gz ERR188383 ERR188428 ERR188454 ERR204916 genes_chrX.tx2gene.tsv); rm -rf $tmp
tmp=$(mktemp -d); cp $R/star_salmon/tximport/salmon.merged.txi.dtu.rds $R/star_salmon/tximport/tximport.tx2gene.tsv $tmp/
printf "sample,condition\nERR188383,GBR\nERR188428,GBR\nERR188454,YRI\nERR204916,YRI\n" > $tmp/samplesheet.csv
(cd $tmp && tar czf $P/testdata/drimseq/drimseq_testdata.tar.gz salmon.merged.txi.dtu.rds samplesheet.csv tximport.tx2gene.tsv); rm -rf $tmp
