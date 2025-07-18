library(tidyverse)
library(readxl)

# Read the MinPvalue
aa <- readRDS("0240912_METSIM_MinPValue_Combined.Rds")
# read TWAS file
dir <- "20240429_METSIM_MinTWASPval_Across_All_Combined.Rds"
df <- readRDS(dir)
# add the pair column 
aa$pair <- paste(aa$Gene_Name, aa$BIOCHEMICAL_NAME, sep = "+")

# Merge by pair
merged <- merge(aa, df, by = "pair")

# Filter essential columns
merged <- merged[, c("BIOCHEMICAL_NAME", "Gene_Name", "LOGPVALUE", "gene", "pvalue")]
saveRDS(merged, "Results_DeepGeneMAP/20240913_CombineTWASMin_DF/20240913_METSIM_TWAS_MinPVal_Combined.Rds")