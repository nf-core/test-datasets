### `Homo_sapiens.GRCh38.102_chr4_subset.gff3`
Full file (`Homo_sapiens.GRCh38.102.gff3`) obtained from running the pipeline with `--build_references`. Then subset with:
``` bash
awk '$1 == "4" && $4 >= 1700000 && $5 <= 54900000' Homo_sapiens.GRCh38.102.gff3 > Homo_sapiens.GRCh38.102_chr4_subset.gff3
```

### `Homo_sapiens.GRCh38.102.chr4_transcripts.fasta`

1. Get IDs from desired genomic regions (`transcript_ids.txt`)
2. Subset full fasta (`Homo_sapiens.GRCh38.102.fasta`) obtained from running the pipeline with flag `--build_references`. Then run:
``` bash
seqkit grep -f transcript_ids.txt Homo_sapiens.GRCh38.102.fasta >  Homo_sapiens.GRCh38.102.chr4_transcripts.fasta
```