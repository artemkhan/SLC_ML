library(data.table)
library("optparse")
library(tidyverse)


#### Iterates through the tissues and chooses the min gene-metabolite pair for each tissue
acc1 <- c() #### Significant gene-metabolite pairs and their tissues
print("Combined Intermediate")
dir <- "/ru-auth/local/home/akhan01/PheWeb/PheWAS_Processed/Results_DeepGeneMAP/20240911_MinPVal_CLSA_TWAS/"
list_files <- list.files(dir)
df <- c()
# the first one would be already sliced
# start the loop
print("Start the loop")
print(Sys.time())


for (i in 1:length(list_files)) {
  curr <- readRDS(paste(dir, list_files[i], sep = ""))
  print(list_files[i])
  df <- rbind(df, curr)
  print(nrow(df))
}

# slicing
print("Slicing")
print(Sys.time()) 
df <- df[, .SD[which.min(pvalue)], by = "pair"]
print(Sys.time()) 

# Save
saveRDS(df, paste("20240911_MinPVal_CLSA_TWAS//20240911_CLSA_MinTWASPval_Across_All_Combined.Rds", sep = ""))

print("Complete")
print(Sys.time()) 

