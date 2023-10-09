# REFERENCES

# cat conf/test/test* | grep params.test_data | cut -d "=" -f 2 | sort -u

# params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/chr21/germlineresources/dbsnp_138.hg38.vcf.gz
#  params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/vcf/dbsnp_146.hg38.vcf.gz
#  params.test_data['homo_sapiens']['genome']['genome_21_chromosomes_dir']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/chr21/sequence/chromosomes.tar.gz
#  params.test_data['homo_sapiens']['genome']['genome_21_fasta']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/chr21/sequence/genome.fasta
#  params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/chr21/sequence/multi_intervals.bed
#  params.test_data['homo_sapiens']['genome']['genome_fasta']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/genome.fasta
#  params.test_data['homo_sapiens']['genome']['genome_interval_list']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/genome.interval_list
#  params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/genome.multi_intervals.bed
#  params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/chr21/germlineresources/gnomAD.r2.1.1.vcf.gz
#  params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/vcf/gnomAD.r2.1.1.vcf.gz
#  params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/chr21/germlineresources/mills_and_1000G.indels.hg38.vcf.gz
#  params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/vcf/mills_and_1000G.indels.vcf.gz
#  params.test_data['homo_sapiens']['genome']['ngscheckmate_bed']
git checkout origin/modules -- data/genomics/homo_sapiens/genome/chr21/germlineresources/SNP_GRCh38_hg38_wChr.bed

# TEST DATA

# cat tests/csv/3.0/* | awk -F',' '{ for( i=1; i<=NF; i++ ) print $i }' | grep https | sort -u

# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam.bai
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam.bai
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam.bai
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam.bai
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/cram/test2.paired_end.recalibrated.sorted.cram
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/cram/test2.paired_end.recalibrated.sorted.cram
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/cram/test2.paired_end.recalibrated.sorted.cram.crai
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/cram/test2.paired_end.recalibrated.sorted.cram.crai
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/cram/test.paired_end.recalibrated.sorted.cram
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/cram/test.paired_end.recalibrated.sorted.cram
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/cram/test.paired_end.recalibrated.sorted.cram.crai
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/cram/test.paired_end.recalibrated.sorted.cram.crai
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/cram/test.paired_end.sorted.cram
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/cram/test.paired_end.sorted.cram
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/cram/test.paired_end.sorted.cram.crai
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/cram/test.paired_end.sorted.cram.crai
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_1.fastq.gz
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/fastq/test_1.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test2_1.fastq.gz
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/fastq/test2_1.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test2_2.fastq.gz
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/fastq/test2_2.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_2.fastq.gz
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/fastq/test_2.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test.umi_1.fastq.gz
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/fastq/test.umi_1.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test.umi_2.fastq.gz
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/fastq/test.umi_2.fastq.gz
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test.baserecalibrator.table
git checkout origin/modules -- data/genomics/homo_sapiens/illumina/gatk/test.baserecalibrator.table
# https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/vcf/test.vcf.gz
git checkout origin/modules -- data/genomics/sarscov2/illumina/vcf/test.vcf.gz
