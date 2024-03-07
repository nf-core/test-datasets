## Analyse the relationship between individuals
df <- read.table("plink.genome", sep = "", header = T)

# Compute for each individuals the number of individuals with PI_HAT > 0.2
lst_ind <- lapply(unique(df$FID1), function(x) sum(df$FID1 == x & df$PI_HAT > 0.2))
names(lst_ind) <- unique(df$FID1)

# convert to data.frame
df_ind <- data.frame(individual = names(lst_ind), count = unlist(lst_ind))

# Display the 20 individuals with the lowest number of individuals with PI_HAT > 0.2
head(df_ind[order(df_ind$count, decreasing = F),],20)
