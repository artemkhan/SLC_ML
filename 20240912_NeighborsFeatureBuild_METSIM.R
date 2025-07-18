# list of features
library(doParallel)
library(tidyr)
library(tidyverse)
registerDoParallel(3)
library(readxl)
library(data.table)
library(openxlsx)
setwd("Neighbor_Features/")
library("optparse")

# Parsing the arguments
option_list = list(
  make_option("--start", action="store", default=NA, type='character',
              help="Path to dataframe of gwas and eQTL summary statistics [required]"))

print("Modified version")
opt = parse_args(OptionParser(option_list=option_list))
start <- opt$start
start <- as.numeric(start)
start_n <- start*10000 + 1
end <- start*10000 + 10000
print(start_n)

# START HERE
#dir <- "20240913_CombineTWASMin_DF/20240913_CLSA_TWAS_MinPVal_Combined.Rds"
# CHANGE for METSIM
dir <- "20240913_CombineTWASMin_DF/20240913_METSIM_TWAS_MinPVal_Combined.Rds"

print(dir)
df <- readRDS(dir)
harmon <- read.xlsx("Reference files/20230807_HarmonizationMetabolitesProper.xlsx")
harmon <- na.omit(harmon)
df <- df[df$BIOCHEMICAL_NAME %in% harmon$`Yin.et.al.,.2022.(original.names).x`,]
df$LOGPVALUE <- 10^(-df$LOGPVALUE)
colnames(df)[1:5] <- c("metab_name", "gene_name", "p_value", "gene", "pvalue")

metabolite_list <- unique(df$metab_name)
if (end > nrow(df)) { end <- nrow(df)}
rxn_map <- read.xlsx("Reference files/20240319_RxnMap.xlsx")

Sys.time()
cadd_scores <- fread("/ru-auth/local/home/akhan01/PheWeb/PheWAS_Processed/Results_CADD/20230307_CADD_MaxScore_PerGene.tsv")
cadd_scores <- cadd_scores[cadd_scores$score > 0,]
together <- df[start_n:end,]
together <- merge(together, cadd_scores, by.x = "gene", by.y = "gene", all.x = TRUE)
gene_list <- unique(together$gene)
dict <- read.xlsx("Reference files/20240326_ChemicalSimilarity/20240418_ChemicalDictionary_Final.xlsx")

## So that no need to perform harmonization but would choose appropriate for each dataset gene naming
if (grepl("METSIM", dir)) { 
  colnames(rxn_map)[colnames(rxn_map) == "BIOCHEMICAL_NAME"] <- "column_metab" } else {
    colnames(rxn_map)[colnames(rxn_map) == "CLSA.cohort.(orignal.names).x"] <- "column_metab"  }

# Check by the rxn map similar 
find_neighbor_by_rxnmap <- function(metabolite_of_interest, gene_of_interest) {
  rxns <- filter(rxn_map, column_metab == metabolite_of_interest)
  filter_rxns <- filter(rxn_map, `Reactions.catalyzed.by.enzyme` %in% rxns$`Reactions.catalyzed.by.enzyme`)
  subset <- df[(df$gene_name == gene_of_interest) & (df$metab_name %in% filter_rxns$column_metab) & !(df$metab_name == metabolite_of_interest), ]
  return(list(min(subset$p_value, na.rm = TRUE), min(subset$pvalue, na.rm = TRUE), nrow(subset)) ) }

# Check the chemical map
find_neighbor_by_chemical_map <- function(metabolite_of_interest, gene_of_interest) {
  ## Check if there are compounds that contain the metabolite of interest
  mets_containing_moi <- metabolite_list[grepl(tolower(metabolite_of_interest), tolower(metabolite_list), fixed = TRUE)]
  mets_containing_moi <- mets_containing_moi[  !(mets_containing_moi == metabolite_of_interest)]
  acc <- c()
  acc_t <- c()
  ## Check if contains the metabolite of interest
  for (i in 1:(length(metabolite_list))) { acc_t <- c(acc_t, grepl(tolower(metabolite_list[i]), metabolite_of_interest, ignore.case = TRUE, fixed = TRUE)) }
  ## Check the dictionary
  for (i in 1:(nrow(dict))) { acc <- c(acc, grepl(tolower(dict$manual[i]), tolower(metabolite_of_interest), ignore.case = TRUE, fixed = TRUE)) }
  roots <- dict$manual[acc]
  acc1 <- c()
  for (j in roots) {  acc1 <- c(acc1, metabolite_list[grepl(tolower(j), tolower(metabolite_list), fixed = TRUE)]) }
  mets_list <- unique(c(acc1, mets_containing_moi, metabolite_list[acc_t]))
  subset <-  df[(df$gene_name == gene_of_interest) & (df$metab_name %in% mets_list) & !(df$metab_name == metabolite_of_interest), ]
  return(list(min(subset$p_value, na.rm = TRUE), min(subset$pvalue, na.rm = TRUE), nrow(subset)) ) }

# Start the loop
print("Start the loop")
n = nrow(together)

Sys.time()
acc <- foreach (i = 1:n, .combine = "rbind") %dopar% {
  curr <- together[i,]
  add <- find_neighbor_by_rxnmap(curr$metab_name, curr$gene_name)
  add1 <- find_neighbor_by_chemical_map(curr$metab_name, curr$gene_name)
  c(add[[1]], add[[2]], add[[3]], add1[[1]], add1[[2]], add1[[3]]) }
Sys.time()

## Form the table
acc <- as_tibble(acc)
ml <- as_tibble(together)
ml <- cbind(ml, neighb = as.numeric(acc$V1))
ml <- cbind(ml, neighb2 = as.numeric(acc$V2))
ml <- cbind(ml, neighb3 = as.numeric(acc$V3))
ml <- cbind(ml, neighb4 = as.numeric(acc$V4))
ml <- cbind(ml, neighb5 = as.numeric(acc$V5))
ml <- cbind(ml, neighb6 = as.numeric(acc$V6))
ml <- as_tibble(ml)

colnames(ml)[4:12] <- c("Min pval", "TWAS", "CADD", "Neighb_Rxn_MinP", "Neighb_TWAS", "Neighb_Rxn_N", "Neighb_Chem_MinP", "Neighb_Chem_TWAS", "Neighb_Chem_N")
print("Start saving")
#saveRDS(as_tibble(ml), paste("20240913_Neighbors_CLSA/CLSA_", start, ".Rds", sep = ""))
saveRDS(as_tibble(ml), paste("20240913_Neighbors_METSIM/METSIM_", start, ".Rds", sep = ""))
print("Completed")

