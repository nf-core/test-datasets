# test-datasets: `dualrnaseq`
This branch contains test data to be used for automated testing with the [nf-core/dualrnaseq](https://github.com/nf-core/dualrnaseq) pipeline.

## Test dataset

The test dataset for the dualrnaseq pipeline consists of subsets of genomes and annotation files of Human (GRCh38.p13) and *Salmonella* Typhimurium (FQ312003), which can be found in the `references/` directory and simulated paired-end reads located in the `PE_reads/` folder. 

## Read simulation procedure

### Reference files

Dual RNA-seq test reads were simulated using a minimal chimeric transcriptome of Human and *Salmonella* Typhimurium.

The transcriptome fasta file for each organism was created using the genome fasta and gff annotation files stored in the `references/` directory. These files were further merged into the chimeric `host_pathogen_transcriptome.fasta` file.  

### Read simulations

Paired end data were simulated using the R package [Polyester (v1.24.0)](https://bioconductor.org/packages/release/bioc/html/polyester.html) following the [Polyester documentation](https://bioconductor.org/packages/release/bioc/vignettes/polyester/inst/doc/polyester.html).

```{r}
library(polyester)

# read minimal chimeric transcriptome file
fasta_file = 'host_pathogen_transcriptome.fasta'
fasta = readDNAStringSet(fasta_file)

# define fold changes
fold_changes = c(4, 4, 1/4,rep(1, 12), 3/5)

# define coverage
readspertx = round(20 * width(fasta) / 100)

# simulate reads:
simulate_experiment('host_pathogen_transcriptome.fasta', reads_per_transcript=readspertx,
    num_reps=c(1,2), fold_changes=fold_changes, outdir='simulated_reads_test')

```

Next, the reformat.sh script from BBTools (v38.79) was used to convert read fast files into fastq files. 

```bash
for filename in simulated_reads_test/*.fasta;
do
	base_name=`basename $filename .fasta`
	out_name="fastq/${base_name}.fq"
	reformat.sh in=$filename out=${out_name} qfake=35
	gzip ${out_name}
done
```
