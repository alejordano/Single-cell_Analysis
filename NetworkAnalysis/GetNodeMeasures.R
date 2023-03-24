library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)

setwd("/data/aconforte/edgeAnalysis")

data <- readRDS("/home/aconforte/Documents/data/AllMerged_NormSubsetCluster.Rds")

data$patient <- sapply(strsplit(data$batch, split = "_", fixed=TRUE), function(x) (x[1]))
data$timepoint <- sapply(strsplit(data$batch, split = "_", fixed=TRUE), function(x) (x[2]))
data$group <- paste0(data$patient, "_", data$clusters)

subpop <- unique(data$group)

#Get number of cells to filter subpopulations:
#Relapse
sp <- c()
cellsDg <- c()
cellsR <- c()

for (i in seq_along(subpop)){
  #number of cells
  data2 <- subset(data, group == subpop[i])
  cellsDg <- append(cellsDg, table(data2$timepoint)["Dg"])
  cellsR <- append(cellsR, table(data2$timepoint)["R"])
  sp <- append(sp, subpop[i])
}

table <- as.data.frame(cbind(sp, cellsDg, cellsR))
table <- table %>% mutate_at(c(2:3), as.numeric)

#filter for good quality paired networks
table2 <- table %>% filter(cellsDg > 100 & cellsR > 100)
sps <- unique(table2$sp)

#Get quantile for all genes at diagnosis:
degree <- c()
betwe <- c()
close <- c()
page <- c()
gene <- c()
sample <- c()
timepoint <- c()

for (i in seq_along(sps)) {
  mea <- read.csv(paste0("./Results/", sps[i], "_Dg_centrality"))
  mea <- mea[order(mea$Degree), ]
  mea$indexD <- seq(1, nrow(mea))
  mea <- mea[order(mea$Betweenness), ]
  mea$indexB <- seq(1, nrow(mea))
  mea <- mea[order(mea$Closeness), ]
  mea$indexC <- seq(1, nrow(mea))
  mea <- mea[order(mea$PAGErank), ]
  mea$indexP <- seq(1, nrow(mea))
  for (y in seq_along(mea)) {
    degree <- append(degree, mea$indexD[y] / nrow(mea))
    betwe <- append(betwe, mea$indexB[y] / nrow(mea))
    close <- append(close, mea$indexC[y] / nrow(mea))
    page <- append(page, mea$indexP[y] / nrow(mea))
    gene <- append(gene, mea$X[y])
    sample <- append(sample, sps[i])
    timepoint <- append(timepoint, "Dg")
  }
}
df <- as.data.frame(cbind(degree, betwe, close, page, gene, sample, timepoint))

#Get quantile for all genes at diagnosis:
degree <- c()
betwe <- c()
close <- c()
page <- c()
gene <- c()
sample <- c()
timepoint <- c()

for (i in seq_along(sps)){
  mea <- read.csv(paste0("./Results/", sps[i], "_R_centrality"))
  mea <- mea[order(mea$Degree), ]
  mea$indexD <- seq(1, nrow(mea))
  mea <- mea[order(mea$Betweenness), ]
  mea$indexB <- seq(1, nrow(mea))
  mea <- mea[order(mea$Closeness), ]
  mea$indexC <- seq(1, nrow(mea))
  mea <- mea[order(mea$PAGErank), ]
  mea$indexP <- seq(1, nrow(mea))
  for (y in seq_along(mea)) {
    degree <- append(degree, mea$indexD[y] / nrow(mea))
    betwe <- append(betwe, mea$indexB[y] / nrow(mea))
    close <- append(close, mea$indexC[y] / nrow(mea))
    page <- append(page, mea$indexP[y] / nrow(mea))
    gene <- append(gene, mea$X[y])
    sample <- append(sample, sps[i])
    timepoint <- append(timepoint, "R")
  }
}
df2 <- as.data.frame(cbind(degree, betwe, close, page, gene, sample, timepoint))

df3 <- rbind(df, df2)

write.csv(df3, "/data/network/AllSubpopAllNodeMeasures.csv")