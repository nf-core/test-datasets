# Using ncbi datasets to download genomes
./datasets download genome accession GCF_024498275.1 GCF_000027325.1 GCF_000733995.1 GCF_017654545.1  --include genome
unzip ncbi_dataset.zip 
mv ncbi_dataset/data/GCF_0*/*.fna .
rm -rf md5sum.txt ncbi_dataset* README.md 
# Simulate 100 5000-bp single end reads
find *.fna | xargs -I{} echo "./wgsim  -N 100 -1 5000 {} {}.fq /dev/null && gzip {}.fq" 
find *.fna | xargs -I{} echo "./wgsim  -N 100 -1 5000 {} {}.fq /dev/null && gzip {}.fq" |sh
