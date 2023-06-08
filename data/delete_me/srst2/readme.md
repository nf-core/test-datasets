### To get fastq files:

Using sra tools `sratoolkit/2.10.9` do the following

```
prefetch SRR9067271
cd SRR9067271/
fastq-dump SRR9067271.sra --split-files
gzip SRR9067271_1.fastq 
gzip SRR9067271_2.fastq
```

More details on the sequencing is found [here](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR9067271)

### Getting MLST_DB on 4/29/2022:

Download each gene fasta file at https://pubmlst.org/bigsdb?db=pubmlst_bfragilis_seqdef&page=downloadAlleles  
Then cat them together with `cat *.fas > MLST_DB.fasta`
Lastly, change the underscore to a dash in the sequence names with `sed -i 's/_/-/' MLST_DB.fasta`

### MLST Profiles 4/29/2022:

Copy the MLST Scheme https://pubmlst.org/bigsdb?db=pubmlst_bfragilis_seqdef&page=downloadProfiles&scheme_id=1 to get the profiles.csv file. 

