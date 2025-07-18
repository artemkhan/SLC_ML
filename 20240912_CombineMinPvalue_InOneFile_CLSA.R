library(doParallel)
registerDoParallel(6)
library(data.table)
library(openxlsx)
library(tidyverse)
dir = "/ru-auth/local/home/akhan01/PheWeb/PheWAS_Processed/Results_Canadian_MinPval_20240326/"
files <- list.files(dir)

### CLSA Combine 
acc <- foreach (i = 1:length(files), .combine = "rbind") %dopar% {
   df <- read.table(paste(dir, files[i], sep = ""))
   met <- gsub("20240326_Canadian_MinPVal_Output_", "", files[i])
   met <- gsub(".txt", "", met)
   df <- cbind(metabolite = rep(met, nrow(df)), df)
   df }
 
acc <- subset(acc, select = -c(effect_allele, other_allele, effect_allele_frequency) )
metabolites <- fread("Reference files/Canadian_MetabList.csv")
metabolites <- metabolites[metabolites$`Ancestry Category` == "European",]
metabolites <- subset(metabolites, select = c(`Reported trait`, `Study Accession`))
acc <- merge(metabolites, acc, by.x = "Study Accession", by.y = "metabolite")
 
print("Saving CLSA")
saveRDS(as_tibble(acc), "20240912_MinPValue_CLSA_METSIM/20240912_CLSA_MinPValue_Combined.Rds")
 
