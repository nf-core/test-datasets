# Centrifuger test database (cfr_files_short.tar.gz)

Files from: https://github.com/mourisl/centrifuger/tree/master/example


1. Truncated reference to 500bp per sequence (to create small db):

awk '/^>/{if(seq) print substr(seq,1,500); print; seq=""} !/^>/{seq=seq$0} END{print substr(seq,1,500)}' ref.fa > small_ref.fa

2. Built index inside Docker on native x86 (Codespaces):

docker run -v $(pwd):/data -w /data quay.io/biocontainers/centrifuger:1.1.0--hf426362_0 centrifuger-build --conversion-table ref_seqid.map --taxonomy-tree  nodes.dmp --name-table names.dmp -r small_ref.fa -o refseq_abv

3. Packaged:

tar -czf cfr_files_short.tar.gz refseq_abv.*.cfr

WARNING: Do NOT build on macOS under Rosetta - produces corrupt .cfr files.

