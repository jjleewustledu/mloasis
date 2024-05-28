# install.packages("anticlust")
library(anticlust)
# install.packages("dplyr")
library(dplyr)
# install.packages("rapportools")
library(rapportools)
setwd("/Volumes/PrecunealSSD/Singularity/OASIS3/NMF_FDG")

subjects = read.csv("anticlust_cn_repeat.csv")
continuous.vars <- subjects[, c("age")]
categorical.vars <- subjects[, c("sex", "apoe4")]

N <- 50
for (i in 1:N) {
  v <- anticlustering(continuous.vars, 2, categories=categorical.vars)
  subjects$Split <- v
  selectA = subjects$Split == 1 & !is.empty(c("Filename"))
  repA = subjects[selectA, c("Filename")]
  selectB = subjects$Split == 2 & !is.empty(c("Filename"))
  repB = subjects[selectB, c("Filename")]
  write.table(repA, paste0("table_cn_repA", i, ".csv"), row.names = F, col.names = F, sep = ',', quote = F)
  write.table(repB, paste0("table_cn_repB", i, ".csv"), row.names = F, col.names = F, sep = ',', quote = F)
}


