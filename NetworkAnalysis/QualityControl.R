library(ggplot2)
library(dplyr)
library(stringr)
library(igraph)
library(ggpubr)

setwd("/home/aconforte/Documents/Iacono")

data <- readRDS("./data/AllMerged_NormSubsetCluster.Rds")
data$patient <- sapply(strsplit(data$batch, split = "_", fixed = TRUE), function(x) (x[1]))
data$timepoint <- sapply(strsplit(data$batch, split = "_", fixed = TRUE), function(x) (x[2]))
data$subpop <- paste0(data$patient, "_", data$clusters)
subpop <- unique(data$subpop)

#Get number of vertices and edges:
#relapse
sp <- c()
cells <- c()
ed <- c()
ve <- c()
ratio <- c()
for (i in seq_along(subpop)) {
  #number of cells
  data2 <- subset(data, group == subpop[i] & timepoint == "R")
  sp <- append(sp, subpop[i])
  cells <- append(cells, ncol(data2))
  if (ncol(data2) > 50) {
    #number of nodes, edges and ratio
    g <- read.graph(paste0("./Results/", subpop[i], "_R_graph.gml"), format = c("gml"))
    ed <- append(ed, ecount(g))
    ve <- append(ve, vcount(g))
    ratio <- append(ratio, ecount(g) / vcount(g))
  } else {
    print(paste0(subpop[i], " number of cells: ", ncol(data2)))
    ed <- append(ed, "0")
    ve <- append(ve, "0")
    ratio <- append(ratio, "0")
  }
}
table1 <- cbind(sp, cells, ed, ve, ratio)
table1$tp <- "R"

#diagnosis
sp <- c()
cells <- c()
ed <- c()
ve <- c()
ratio <- c()
for (i in seq_along(subpop)) {
  #number of cells
  data2 <- subset(data, group == subpop[i] & timepoint == "Dg")
  sp <- append(sp, subpop[i])
  cells <- append(cells, ncol(data2))
  if (ncol(data2) > 50) {
    #number of nodes, edges and ratio
    g <- read.graph(paste0("./Results/", subpop[i], "_Dg_graph.gml"), format = c("gml"))
    ed <- append(ed, ecount(g))
    ve <- append(ve, vcount(g))
    ratio <- append(ratio, ecount(g) / vcount(g))
  } else {
    print(paste0(subpop[i], " number of cells: ", ncol(data2)))
    ed <- append(ed, "0")
    ve <- append(ve, "0")
    ratio <- append(ratio, "0")
  }
}
table2 <- cbind(sp, cells, ed, ve, ratio)
table2$tp <- "Dg"

table <- rbind(table1, table2)

table <- table %>% mutate_at(c(2:5), as.numeric)

#Get degree distribution. Is it a power law?
#SupplementaryFigure2 : Degree distribution
dfv <- table %>% filter(ve < 500 & ratio != 0)
dfv <- dfv[order(dfv$ve), ]
par(mfrow = c(4, 4))
for (i in seq_along(dfv)){ # had to run maximum 16 plots per time.
  dg <- read.csv(paste0("./Results/", dfv$sp[i], "_", dfv$tp[i], "_centrality"))
  x <- dg$Degree
  hist(x, main = paste0(dfv$sp[i], "_", dfv$tp[i], "\nRatio E/V = ", dfv$ratio[i], "\nN vertices = ", dfv$ve[i]), xlab = "Degree")
}