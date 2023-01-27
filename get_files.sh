# cat conf/test.config | awk -F',' '{ for( i=1; i<=NF; i++ ) print $i }' | grep https | sort -u 
# https://github.com/nf-core/test-datasets/raw/rnaseq/reference/gfp.fa.gz
git checkout origin/rnaseq -- reference/gfp.fa.gz
# https://github.com/nf-core/test-datasets/raw/rnaseq/reference/bbsplit_fasta_list.txt
git checkout origin/rnaseq -- reference/bbsplit_fasta_list.txt
# https://github.com/nf-core/test-datasets/raw/rnaseq/reference/genome.fasta
git checkout origin/rnaseq -- reference/genome.fasta
# https://github.com/nf-core/test-datasets/raw/rnaseq/reference/genes.gff.gz
git checkout origin/rnaseq -- reference/genes.gff.gz
# https://github.com/nf-core/test-datasets/raw/rnaseq/reference/genes.gtf.gz
git checkout origin/rnaseq -- reference/genes.gtf.gz
# https://github.com/nf-core/test-datasets/raw/rnaseq/reference/hisat2.tar.gz
git checkout origin/rnaseq -- reference/hisat2.tar.gz
# https://github.com/nf-core/test-datasets/raw/rnaseq/reference/rsem.tar.gz
git checkout origin/rnaseq -- reference/rsem.tar.gz
# https://github.com/nf-core/test-datasets/raw/rnaseq/reference/salmon.tar.gz
git checkout origin/rnaseq -- reference/salmon.tar.gz
# https://github.com/nf-core/test-datasets/raw/rnaseq/reference/transcriptome.fasta
git checkout origin/rnaseq -- reference/transcriptome.fasta

# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/samplesheet/v3.10/samplesheet_test.csv
git checkout origin/rnaseq -- samplesheet/v3.10/samplesheet_test.csv
# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz
git checkout origin/rnaseq -- testdata/GSE110004/SRR6357070_1.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz
git checkout origin/rnaseq -- testdata/GSE110004/SRR6357070_2.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357071_1.fastq.gz
git checkout origin/rnaseq -- testdata/GSE110004/SRR6357071_1.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357071_2.fastq.gz
git checkout origin/rnaseq -- testdata/GSE110004/SRR6357071_2.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_1.fastq.gz
git checkout origin/rnaseq -- testdata/GSE110004/SRR6357072_1.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_2.fastq.gz
git checkout origin/rnaseq -- testdata/GSE110004/SRR6357072_2.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357073_1.fastq.gz
git checkout origin/rnaseq -- testdata/GSE110004/SRR6357073_1.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357074_1.fastq.gz
git checkout origin/rnaseq -- testdata/GSE110004/SRR6357074_1.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357075_1.fastq.gz
git checkout origin/rnaseq -- testdata/GSE110004/SRR6357075_1.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_1.fastq.gz
git checkout origin/rnaseq -- testdata/GSE110004/SRR6357076_1.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_2.fastq.gz
git checkout origin/rnaseq -- testdata/GSE110004/SRR6357076_2.fastq.gz
