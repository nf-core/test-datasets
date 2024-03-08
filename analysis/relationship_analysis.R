## Charge library
library(optparse)

## Parse the arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input file (plink.genome)"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output directory (analysis)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## Read data
print(paste("Reading file:", opt$input))
df <- read.table(opt$input, sep = "", header = T)

# Compute for each individuals the number of individuals with PI_HAT > 0.2
lst_ind <- lapply(unique(df$FID1), function(x) sum(df$FID1 == x & df$PI_HAT > 0.2))
names(lst_ind) <- unique(df$FID1)

# convert to data.frame
df_ind <- data.frame(individual = names(lst_ind), count = unlist(lst_ind))

# Display the 20 individuals with the lowest number of individuals with PI_HAT > 0.25
print("Individuals with the lowest number of individuals with PI_HAT > 0.25")
head(df_ind[order(df_ind$count, decreasing = F),],20)

# Selection of the individuals
ind_selected <- c("NA12878", "NA19401", "NA20359")
print(paste("Selected individuals:", paste(ind_selected, collapse = ", ")))

# Show the individuals with PI_HAT > 0.25 for the selected individuals
print("Individuals with PI_HAT > 0.25 for the selected individuals")
df_ind <- df[df$FID1 %in% ind_selected & df$PI_HAT > 0.25, c("IID1", "IID2", "PI_HAT")]
df_ind

# Write the output
write.csv(df_ind, file = paste0(opt$output, "/rel_ind_sel.csv"), row.names = F)
write.table(
  ind_selected,
  file = paste0(opt$output, "/selected_individuals.txt"),
  row.names = F, col.names = F, quote = F
)
write.table(
  unique(c(ind_selected, df_ind$IID2)),
  file = paste0(opt$output, "/all_rel_individuals.txt"),
  row.names = F, col.names = F, quote = F
)
