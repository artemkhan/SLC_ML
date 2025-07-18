library(data.table)
library("optparse")
library(tidyverse)

Sys.time()
# Parsing the arguments
option_list = list(
  make_option("--index", action="store", default=NA, type='character',
              help="Path to dataframe of gwas and eQTL summary statistics [required]")
)


opt = parse_args(OptionParser(option_list=option_list))
i = as.numeric(unlist(opt$index))

#### Iterates through the tissues and chooses the min gene-metabolite pair for each tissue
acc1 <- c() #### Significant gene-metabolite pairs and their tissues
dir <- "/ru-auth/local/home/akhan01/PheWeb/PheWAS_Processed/Results_Summary/"
list_files <- list.files(dir)
list_files <- list_files[grepl("Canadian", list_files)]
df <- c()
# the first one would be already sliced
# start the loop
print("Start the loop")
print(Sys.time())
k = 0
for (j in 1:7) {
  k = i*7 + j
  curr <- readRDS(paste(dir, list_files[k], sep = ""))
  print(list_files[k])
  curr <- curr[, c("metab_name", "gene", "gene_name", "pvalue")]
  df <- rbind(df, curr)
  print(nrow(df))
}

# slicing
print("Slicing")
print(Sys.time()) 
df$pair <- paste(df$gene_name, df$metab_name, sep = "+")
df <- df[, .SD[which.min(pvalue)], by = "pair"]
print(Sys.time()) 

# Save
saveRDS(df, paste("20240911_MinPVal_CLSA_TWAS/20240911_CLSA_MinTWASPval", i, ".Rds", sep = ""))

print("Complete")
print(Sys.time()) 

