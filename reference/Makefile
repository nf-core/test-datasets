uncompressed=$(wildcard *.gtf *.gff *.fa *.fasta)

gz:
	for F in ${uncompressed}; do echo $$F ; gzip -c $$F > $$F.gz ; done
