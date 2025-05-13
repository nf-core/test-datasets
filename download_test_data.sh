# Download somalier sites
wget https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz -O data/sites.hg38.vcf.gz
tabix data/sites.hg38.vcf.gz

# Download CoLoRSdb 
wget https://zenodo.org/records/11511513/files/CoLoRSdb.GRCh38.v1.0.0.pbsv.jasmine.vcf.gz -O data/CoLoRSdb.GRCh38.v1.0.0.pbsv.jasmine.vcf.gz
wget https://zenodo.org/records/11511513/files/CoLoRSdb.GRCh38.v1.0.0.pbsv.jasmine.vcf.gz.tbi -O data/CoLoRSdb.GRCh38.v1.0.0.pbsv.jasmine.vcf.gz.tbi

# Download reference genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -O data/GRCh38_no_alt_analysis_set.fasta.gz
gunzip -c data/GRCh38_no_alt_analysis_set.fasta.gz > data/GRCh38_no_alt_analysis_set.fasta

# Download ONT HG002 haplotagged bam
aws s3 cp s3://ont-open-data/giab_2023.05/analysis/variant_calling/hg002_sup_60x/hg002.haplotagged.bam data/hg002.haplotagged.bam --no-sign-request 
aws s3 cp s3://ont-open-data/giab_2023.05/analysis/variant_calling/hg002_sup_60x/hg002.haplotagged.bam.bai data/hg002.haplotagged.bam.bai --no-sign-request 

# Add aligned PacBio HG002, HG003, HG004 manually, from Nallo run, as
# data/HG002_aligned_haplotagged.bam (m84011_220902_175841_s1.hifi_reads.bam)
# data/HG003_aligned_haplotagged.bam (m84010_220919_235306_s2.hifi_reads.bam)
# data/HG004_aligned_haplotagged.bam (m84010_220919_232145_s1.hifi_reads.bam)
# data/HG002_aligned_haplotagged.bam.bai
# data/HG003_aligned_haplotagged.bam.bai
# data/HG004_aligned_haplotagged.bam.bai

# Download VEP cache

wget https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz -O data/homo_sapiens_vep_110_GRCh38.tar.gz