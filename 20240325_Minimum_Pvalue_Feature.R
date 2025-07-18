# select the min pvalue 
library(tidyverse)
library(doParallel)
library(data.table)
library(openxlsx)
registerDoParallel(4)
library("optparse")

Sys.time()
# Parsing the arguments
option_list = list(
  make_option("--index", action="store", default=NA, type='character',
              help="Path to dataframe of gwas and eQTL summary statistics [required]")
)


opt = parse_args(OptionParser(option_list=option_list))
i = as.numeric(unlist(opt$index))

s <- list.files("/ru-auth/local/home/akhan01/PheWeb/PheWAS/share.sph.umich.edu/metsim_metabqtl/")
all_files <- list.files(path = "/ru-auth/local/home/akhan01/PheWeb/PheWAS/share.sph.umich.edu/metsim_metabqtl", full.names = TRUE, recursive = TRUE, pattern = "gz")
print(all_files[i])
genome <- read.xlsx("Reference files/20240325_Ensembl96_Coordinates.xlsx")

list_genes <- tibble(genome$V1, genome$V4, genome$V5, genome$`gene_name`)

metabolite <- fread(all_files[i], header = TRUE)

print("Iterating through the genes")
accum <- foreach (i=1:nrow(list_genes), .combine = "rbind") %dopar% {
    gene_name <- list_genes[i, 4]
    start <- as.numeric(list_genes[i, 2])
    end <- as.numeric(list_genes[i, 3] )
    chrom <- as.numeric(list_genes[i, 1])
    look <- metabolite[(metabolite$BEG > (start - 15000)) & (metabolite$END < (end + 15000) ) & (metabolite$CHROM == chrom),]
    if (!(nrow(look) == 0)) { cbind(gene_name, look[look$LOGPVALUE == max(look[["LOGPVALUE"]]),]) 
    }  
}

metname <- gsub("/ru-auth/local/home/akhan01/PheWeb/PheWAS/share.sph.umich.edu/metsim_metabqtl/", "", all_files[i])
metname <- gsub("\\_.*", "", metname)
colnames(accum) <- c(paste("Gene_name", metname), "CHROM", "BEG", "END", "MARKER_ID", "NEA","EA", "N", "EAC", "MAF", "BETA", "SEBETA", paste("LOGPVALUE", metname))                     
first <- paste("Results_METSIM_MinPval_March2024/20240325_METSIM_MinPVal_Output", metname, sep = "_")
full <- paste(first, ".txt", sep = "")

print("Writing an Output")
print(full)
# write an output
write.table(accum, full)
print("Completed")

Sys.time()
