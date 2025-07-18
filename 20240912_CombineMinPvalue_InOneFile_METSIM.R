library(doParallel)
registerDoParallel(6)
library(data.table)
library(openxlsx)
library(tidyverse)

### METSIM 
dir <- "/ru-auth/local/home/akhan01/PheWeb/PheWAS_Processed/Results_METSIM_MinPval_March2024/"
files <- list.files(dir)
n = length(list.files(dir))
acc <- foreach (i = 1:n, .combine = "rbind") %dopar% {
  df <- read.table(paste(dir, files[i], sep = ""))
  colnames(df)[grepl( "Gene_name", colnames(df))] <- "Gene_Name"
  colnames(df)[grepl("LOGPVALUE", colnames(df))] <- "LOGPVALUE"
  met <- gsub("20240325_METSIM_MinPVal_Output_", "", files[i])
  met <- gsub(".txt", "", met)
  df <- cbind(metabolite = rep(met, nrow(df)), df)
  df }

acc <- subset(acc, select = -c(NEA, EA, N, EAC, MAF) )

metabolites <- fread("Reference files/tab_metabolites.txt")

acc <- merge(metabolites, acc, by.x = "CID", by.y = "metabolite")

print("Saving METSIM")
saveRDS(as_tibble(acc), "Results_DeepGeneMAP/20240912_MinPValue_CLSA_METSIM/20240912_METSIM_MinPValue_Combined.Rds")

