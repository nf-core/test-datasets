#Using same IDs as present in the genotype file we make an artificial phenotyp and covariate file

#path to genotype vcf file
pathToVcf <- "data/data_shrink_chunk_4500/chr22.vcf.bgz"

#Get sample name
header <- system(paste("zcat ", pathToVcf, " | grep '#CHROM'", sep=""), intern=TRUE)
header <- unlist(strsplit(header, split="\t"))
#sampleIDs
ids <- header[10:length(header)]

#make a phenotype
#Here, a not too bad tutorial on different techniques on how to simulate data
# https://aosmith.rbind.io/2018/08/29/getting-started-simulating-data/

#Here I used the blog's proposed way of simulating data for a regression analysis
#We will use the generated data slightly different, but hopefully good enough to
# actually get som results
n <- length(ids)
set.seed(16)
y = rnorm(n = n, mean = 0, sd = 1)
x1 = runif(n = n, min = 1, max = 2)
x2 = runif(n = n, min = 200, max = 300)

# Write this to one phenodata file and one covardata file 
# first column, unique ids, second column family ids, remaining columns are 
# phenotyp or covariate columns (here individual IDs are family IDs)
phe <- data.frame(ids=ids, fam=ids, pheno=y)
cov <- data.frame(ids=ids, fam=ids, cov1=x1, cov2=x2)

#write to file
outDir <- "data/data_phenotypes_and_covariates"
write.table(phe, file = paste(outDir,"/example1.pheno", sep=""), append = FALSE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = FALSE, qmethod = c("escape", "double"),
                 fileEncoding = "")
write.table(cov, file = paste(outDir,"/example1.covar", sep=""), append = FALSE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = FALSE, qmethod = c("escape", "double"),
                 fileEncoding = "")


