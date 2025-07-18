library(openxlsx)

CLSA <- readRDS("20240919_FeaturesCombined/20240919_CLSA_Combined.Rds")
METSIM <- readRDS("20240919_FeaturesCombined/20240919_METSIM_Combined.Rds")


## Harmonize metabolites 
harmon <- read.xlsx("Reference files/20230807_HarmonizationMetabolitesProper.xlsx")
harmon <- na.omit(harmon)

CLSA <- merge(CLSA, harmon, by.x = "metab_name", by.y = "CLSA.cohort.(orignal.names).x")
METSIM <- merge(METSIM, harmon, by.x = "metab_name", by.y = "Yin.et.al.,.2022.(original.names).x")


CLSA$pair <- paste(CLSA$gene_name, CLSA$harmonized.metabolite.names, sep = "+")
METSIM$pair <- paste(METSIM$gene_name, METSIM$harmonized.metabolite.names, sep = "+")

colnames(CLSA) <- paste(colnames(CLSA), "CLSA", sep = "_")
colnames(METSIM) <- paste(colnames(METSIM), "METSIM", sep = "_")
merged <- merge(CLSA, METSIM, by.x = "pair_CLSA", by.y = "pair_METSIM")
merged <- unique(merged)

merged_new <- merged[, !names(merged) %in% c("gene_name_METSIM", "metab_name_METSIM", "CLSA.cohort.(orignal.names).x_METSIM", "Yin.et.al.,.2022.(original.names).x_CLSA", "gene_METSIM",
                     "harmonized.metabolite.names_METSIM",  "Yin.et.al.,.2022.(original.names).x_CLSA")]
### merged -> is the largest one
saveRDS(merged_new, "20240919_FeaturesCombined/20240919_Final_Combined_LessColumns.Rds")

merged_without_unknown <- merged_new[!(grepl("X-", merged_new$metab_name_CLSA)),]

saveRDS(merged_without_unknown, "20240919_FeaturesCombined/20240919_Final_Combined_LessColumns_without_unknown.Rds")
