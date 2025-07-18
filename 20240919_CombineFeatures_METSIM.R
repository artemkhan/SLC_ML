files <- list.files("20240913_Neighbors_METSIM/")
acc <- c()
n <- length(files)
for (i in 1:n) {
  r <- readRDS(paste("20240913_Neighbors_METSIM/", files[i], sep = ""))
  acc <- rbind(acc, r)
}

saveRDS(acc, "20240919_FeaturesCombined/20240919_METSIM_Combined.Rds")

print("Completed")
