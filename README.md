# test-datasets: `rarevariantburden`

This branch contains test data to be used for automated testing with the [nf-core/rarevariantburden](https://github.com/nf-core/rarevariantburden) pipeline.

## Content of this repository
`case/`: Input files needed for the pipeline from the case/patient cohort

`case/samples.1KG.chr21-22.vcf.gz`: Joined called and VQSR applied vcf file from the case cohort, here we have a very small subset of joined called vcf file containg part of chr 21 and chr 22 exome region. Here we used 25 test samples from 1000 Genomes Project (build GRCh37) as case data. This is the input of pipeline parameter `caseJointVCF`

`case/samples.1KG.chr21-22.vcf.gz.tbi`: The tabix index file for the joined called vcf file

`case/samples.txt`: One column text file containing list of samples, one sample ID per line. Here we used 25 test samples from 1000 Genomes Project (build GRCh37) as case data. This is the input of pipeline parameter `caseSample`

`case/annotationFiles.csv`: List of annotated files for test profile. This is the input of pipeline parameter `caseAnnotatedVCFFileList`

`case/annotationGDSFiles.csv`: List of annotated GDS files for test profile. This is the input of pipeline parameter `caseAnnotationGDSFileList`

`case/genotypeGDSFiles.csv`: List of genotype GDS files for test profile. This is the input of pipeline parameter `caseGenotypeGDSFileList`

`case/21.annotated.vcf.gz`: Pre-annotated case VCF file for chr 21, we used the pre-annotated file to skip the annotation steps in the pipelie for testing purpose.

`case/21.annotated.vcf.gz.tbi`: The tabix index file for the pre-annotated case VCF file for chr 21

`case/21.annotated.vcf.gz.gds`: The GDS format file for the pre-annotated case VCF file for chr 21, we used the GDS file to skip the VCF to GDS conversion steps in the pipelie for testing purpose.

`case/22.annotated.vcf.gz`: Pre-annotated case VCF file for chr 22, we used the pre-annotated file to skip the annotation steps in the pipelie for testing purpose.

`case/22.annotated.vcf.gz.tbi`: The tabix index file for the pre-annotated case VCF file for chr 22

`case/22.annotated.vcf.gz.gds`: The GDS format file for the pre-annotated case VCF file for chr 22, we used the GDS file to skip the VCF to GDS conversion steps in the pipelie for testing purpose.

`case/21.biallelic.leftnorm.ABCheck.vcf.gz.gds`: The GDS format for the left normalized case VCF file for chr 21, we used this to skip the normalization and convert nomalized VCF file to GDS format steps in the pipeline.

`case/22.biallelic.leftnorm.ABCheck.vcf.gz.gds`: The GDS format for the left normalized case VCF file for chr 22, we used this to skip the normalization and convert nomalized VCF file to GDS format steps in the pipeline.

`case/casePopulation.txt`: The predicted ancestry for each sample, if not specified pipeline will estimate the ancestry/ethnicity of each sample using gnomAD classifier. This is the input of pipeline parameter `casePopulation`

`control/`: Input files needed for the pipeline from the control dataset (here we used gnomAD v2 exome data as control dataset)

`control/controlAnnotationGDS.csv`: List of pre-annotated and GDS converted VCF file from gnomAD v2 exome.

`control/controlGenotypeGDS.csv`: List of normalized and GDS converted VCF file from gnomAD v2 exome.

`control/annotation/chr21.annovar.vep.vcf.gz.gds`: The pre-annotated and GDS converted VCF file from gnomAD v2 exome. Here we used a small exome region for chr 21 as the control data.

`control/annotation/chr22.annovar.vep.vcf.gz.gds`: The pre-annotated and GDS converted VCF file from gnomAD v2 exome. Here we used a small exome region for chr 22 as the control data.

`control/genotypeCount/gnomad.21.vcf.bgz.gds`: The normalized and GDS converted VCF file from gnomAD v2 exome containg only the genotype count. Here we used a small exome region for chr 21 as the control data.

`control/genotypeCount/gnomad.22.vcf.bgz.gds`: The normalized and GDS converted VCF file from gnomAD v2 exome containg only the genotype count. Here we used a small exome region for chr 22 as the control data.

`control/reference.fasta.gz`: The reference file, here we used a small subset of chr 21 and chr 22 exome region as the reference.

`control/reference.fasta.gz.fai`: The FAI index file for the reference.

`control/reference.fasta.gz.gzi`: The GZI index file for the reference.

`control/coverage10x.bed.gz`: File containing coverage information for the gnomAD v2 exome data.

`control/full_vs_gnomAD.p0.05.OR1.ignoreEthnicityInLD.rds`: File containing precomputed high LD (Linkage disequilibrium) variant pairs in controls (gnomAD v2 exome data)

`control/gnomAD.exclude.allow.segdup.lcr.v3.txt.gz`: File containing list of variants which will be exclude during the burden test, here we provided a list of variants in gnomAD control data which failed VQSR filter

`control/stratified_config_gnomad.txt`:  File specifying the AC AN info ID which will be used for stratification, here we used ethnicity stratified analysis and used AC AN info field IDs for each population used in the gnomAD VCF files.
