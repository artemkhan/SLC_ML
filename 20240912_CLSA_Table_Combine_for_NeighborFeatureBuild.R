library(tidyverse)

# Read the MinPvalue
aa <- readRDS("20240912_CLSA_MinPValue_Combined.Rds")
# read TWAS file
dir <- "20240911_MinPVal_CLSA_TWAS/20240911_CLSA_MinTWASPval_Across_All_Combined.Rds"
df <- readRDS(dir)
# add the pair column 
aa$pair <- paste(aa$Gene_name, aa$`Reported trait`, sep = "+")

# Merge by pair
merged <- merge(aa, df, by = "pair")

# Filter essential columns
merged <- merged[, c("metab_name", "Gene_name", "p_value", "gene", "pvalue")]
saveRDS(merged, "20240913_CombineTWASMin_DF/20240913_CLSA_TWAS_MinPVal_Combined.Rds")
