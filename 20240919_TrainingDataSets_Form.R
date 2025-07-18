set.seed(42)
library(tidyverse)
library(arrow)
## Chose positive associatiions
df <- readRDS("20240919_FeaturesCombined/20240919_Final_Combined_LessColumns_without_unknown.Rds")
train_associations <- read.xlsx("Reference files/20240311_TrainingAssociations.xlsx")
df$pair <- paste(df$metab_name_CLSA, df$gene_name_CLSA, sep = "+")

# remove the ones that are not unique 
df_unique <- df %>%
  group_by(pair) %>%
  filter(n() == 1) %>%
  ungroup()

df <- df_unique
train_pos <- df[df$pair %in% train_associations$value,]
train_pos$label <- rep(1, nrow(train_pos))
# Negative random associations 
train_neg <- df[!(df$pair %in% train_associations$value),]
train_neg <- df[sample(1:nrow(df), nrow(train_pos)),]
train_neg$label <- rep(-1, nrow(train_neg))
train_final <- rbind(train_pos, train_neg)


ind = !(colnames(train_final) %in% 
          c("pair_CLSA", "metab_name_CLSA", "gene_CLSA", "gene_name_CLSA", "harmonized.metabolite.names_CLSA", "pair"))
train_final <- train_final[, ..ind]

colnames(train_final)[colnames(train_final) == "label"] <- "target" 
write.csv(train_final, "20240919_FeaturesCombined/20240923_TrainingDataSet_Clean.csv")


## Form on the unseen data
## Chose positive associatiions

# read train
train <- readRDS("20240919_FeaturesCombined/20240923_TrainingDataSet_Clean.Rds")
test <- df[!(df$pair %in% train$pair), ]

#test_final <- test[, !(names(test) %in% 
#                                 c("pair_CLSA", "metab_name_CLSA", "gene_CLSA", "gene_name_CLSA", "harmonized.metabolite.names_CLSA", "pair"))]

write.csv(test_final, "20240919_FeaturesCombined/20240923_TestPredictDataSet_Clean.csv")
ind <- c(!(colnames(test) %in% c("pair_CLSA", "metab_name_CLSA", "gene_CLSA", "gene_name_CLSA", "harmonized.metabolite.names_CLSA")))
test <- test[, ..ind]

write_parquet(test, "20240919_FeaturesCombined/test_with_pairs_Clean.parquet")

