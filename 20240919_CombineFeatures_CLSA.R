files <- list.files("/20240913_Neighbors_CLSA/")
acc <- c()
n <- length(files)
for (i in 1:n) {
  r <- readRDS(paste("20240913_Neighbors_CLSA/", files[i], sep = ""))
  acc <- rbind(acc, r)
}

saveRDS(acc, "20240919_FeaturesCombined/20240919_CLSA_Combined.Rds")

print("Completed")
