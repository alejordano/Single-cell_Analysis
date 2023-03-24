#! /home/aconforte/.conda/envs/aobase/bin/Rscript

library(ggplot2)
library(dplyr)
library(stringr)
library(igraph)
library(data.table)

df <- c()

setwd("/data/aconforte/edgeAnalysis/Results")

# Get list of all networks build for subpopulations at both timepoints:
lista <- list.files(pattern = "_correlations")
patient <- sapply(strsplit(lista, split = "_", fixed = TRUE), function(x) (x[1]))
cluster <- sapply(strsplit(lista, split = "_", fixed = TRUE), function(x) (x[2]))
timepoint <- sapply(strsplit(lista, split = "_", fixed = TRUE), function(x) (x[3]))
sample <- paste0(patient, "_", cluster)

# Select samples with networks available for a subpopulation at both diagnosis and relapse
table <- as.data.frame(cbind(patient, cluster, timepoint, sample))
table2 <- table %>% group_by(sample) %>% mutate(ntp = n_distinct(timepoint))
fil <- as.data.frame(unique(table2[, 4:5]))
fil2 <- fil %>% filter(ntp == 2)
table4 <- table %>% filter(sample %in% fil2$sample)
list <- paste0(table4$sample, "_", table4$timepoint, "_correlations")

# Get all interactions with weight > 0.6 in a Df format
for (w in seq_along(list)) {
  print("reading file")
  file <- fread(list[w])
  file <- as.data.frame(file)
  print(paste0("working on ", list[w], " / ", w, " out of ", length(list)))
  colnames(file) <- gsub("-", "\\.", colnames(file))
  file$V1 <- gsub("-", "\\.", file$V1)
  print("Substituted . by - ")
  file2 <- file[, -1]
  rownames(file2) <- file[, 1]
  m <- as.matrix(file2)
  m[lower.tri(m, diag = TRUE)] <- NA
  file3 <- as_data_frame(graph_from_adjacency_matrix(m, weighted = TRUE))
  file4 <- file3 %>% filter(!is.na(weight) & (abs(weight) > 0.6))
  print("building DF")
  file4$pat <- sapply(strsplit(list[w], split = "_", fixed=TRUE), function(x) (x[1]))
  file4$cluster <- sapply(strsplit(list[w], split = "_", fixed=TRUE), function(x) (x[2])) 
  file4$timepoint <- sapply(strsplit(list[w], split = "_", fixed=TRUE), function(x) (x[3])) 
  df <- rbind(df, file4)
}

write.csv(df, "./tableInteractionsAllPaired")