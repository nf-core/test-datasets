wget -nd ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR044/ERR044595/ERR044595_1.fastq.gz
wget -nd ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR044/ERR044595/ERR044595_2.fastq.gz
wget -nd ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR064/ERR064912/ERR064912_1.fastq.gz
wget -nd ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR064/ERR064912/ERR064912_2.fastq.gz

mkdir downsampled;
for f in ERR044595_1.fastq.gz ERR044595_2.fastq.gz ERR064912_1.fastq.gz ERR064912_2.fastq.gz; do
    seqtk sample $f 1000000 | gzip > downsampled/$(basename $f | sed -e 's,\(_[12]\).fastq.gz,_1M\1.fastq.gz,'); 
done

mv downsampled/* .
rm ERR044595_1.fastq.gz ERR044595_2.fastq.gz ERR064912_1.fastq.gz ERR064912_2.fastq.gz
