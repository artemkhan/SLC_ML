# select the min pvalue 
library(tidyverse)
library(doParallel)
library(openxlsx)
library(data.table)
registerDoParallel(5)
library("optparse")

Sys.time()
# Parsing the arguments
option_list = list(
  make_option("--index", action="store", default=NA, type='character',
              help="Path to dataframe of gwas and eQTL summary statistics [required]")
)


opt = parse_args(OptionParser(option_list=option_list))
i = as.numeric(unlist(opt$index))

all_files <- list.files(path = "/ru-auth/local/home/akhan01/PheWeb/PheWAS/Canadian_SummaryData/", full.names = TRUE, recursive = TRUE, pattern = "gz")
print(all_files[i])
genome <- read.xlsx("Reference files/20240325_Ensembl96_Coordinates.xlsx")
metname <- gsub("/ru-auth/local/home/akhan01/PheWeb/PheWAS/Canadian_SummaryData//", "", all_files[i])
metname <- gsub("_buildGRCh38.tsv.gz", "", metname)

list_genes <- tibble(genome$V1, genome$V4, genome$V5, genome$`gene_name`)

metabolite <- fread(all_files[i], header = TRUE)
metabolite <- as_tibble(metabolite)

print("Iterating through the genes")
Sys.time()
accum <- foreach (i=1:nrow(list_genes), .combine = "rbind") %dopar% {
    gene_name <- list_genes[i, 4]
    start <- as.numeric(list_genes[i, 2])
    end <- as.numeric(list_genes[i, 3] )
    chrom <- as.numeric(list_genes[i, 1])
    look <- metabolite[(metabolite$base_pair_location > (start - 15000)) & (metabolite$base_pair_location < (end + 15000) ) 
                       & (metabolite$chromosome == chrom),]
    look <- look[!(is.na(look$chromosome)),] 
    if (!(nrow(look) == 0)) { cbind(gene_name, look[look$p_value == min(look[["p_value"]]),]) 
    }  }

Sys.time()

colnames(accum) <- c("Gene_name", colnames(metabolite))                     
first <- paste("Results_Canadian_MinPval_20240326/20240326_Canadian_MinPVal_Output", metname, sep = "_")
full <- paste(first, ".txt", sep = "")

print("Writing an Output")
print(full)
# write an output
write.table(accum, full)
print("Completed")

Sys.time()
