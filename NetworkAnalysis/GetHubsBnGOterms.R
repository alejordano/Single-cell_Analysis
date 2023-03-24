library(dplyr)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

setwd("~/Documents/Iacono/")

final <- fread("/data/aconforte/edgeAnalysis/Results/tableInteractionsAllPaired")

dgfinal <- final[(timepoint == "Dg"), ]
relfinal <- final[(timepoint == "R"), ]

## previously identified hubs and bn:
nodes <- c("SCO1", "ATPAF2", "SRFBP1", "SLC7A6", "NNT", "COIL", "M6PR", "CDIP1", "EMID1",
          "DALRD3", "MTAP", "CEP152", "PUS1", "MS4A7", "HSPB11", "PRORP", "MDH1", "SVBP",
          "GUCD1", "TTC17", "CHM", "MRPL10", "PSENEN", "TOP1MT", "NEMP1", "INPP4A", "UBE2L6",
          "BRCC3", "DBR1", "KIF22", "MRPL42", "NOTCH2NLC")
#Diagnosis
dgfinal$from2 <- str_replace(dgfinal$from, "^MT[.]", "")
dgfinal$from2 <- str_replace(dgfinal$from2, "[.]AS1", "")
dgfinal$to2 <- str_replace(dgfinal$to, "^MT[.]", "")
dgfinal$to2 <- str_replace(dgfinal$to2, "[.]AS1", "")

dffinal <- c()
for (i in seq_along(nodes)) {
  print(paste0("Working on dg hub ", i, " out of ", length(nodes)))
  df <- dgfinal[from2 == nodes[i]]
  list <- unique(df$to2)
  df <- dgfinal[to2 == nodes[i]]
  list <- append(list, unique(df$from2))
  pos <- enrichGO(gene = list, universe = names(geneList),
        keyType = "SYMBOL", OrgDb = org.Hs.eg.db, ont = "BP",
        pAdjustMethod = "BH", pvalueCutoff  = 0.01,
        qvalueCutoff  = 0.05, readable = FALSE)
  if (!is.null(pos)) {
    pos2 <- simplify(pos, cutoff = 0.7, by = "p.adjust",
            select_fun = min, measure = "Wang", semData = NULL)
    pos3 <- as.data.frame(pos2@result)
    pos3$regulator <- nodes[i]
    dffinal <- rbind(dffinal, pos3)
  } else {
    pos3 <- cbind(rep("NA", n = 9), nodes[i])
    dffinal <- rbind(dffinal, pos3)
  }
}

write.csv(dffinal, "/data/network/GOEA_NodesHubsBNDgAllGOSimplify")

### relapse
relfinal$from2 <- str_replace(relfinal$from, "^MT[.]", "")
relfinal$from2 <- str_replace(relfinal$from2, "[.]AS1", "")
relfinal$to2 <- str_replace(relfinal$to, "^MT[.]", "")
relfinal$to2 <- str_replace(relfinal$to2, "[.]AS1", "")

dffinal <- c()

for (i in seq_along(nodes)) {
  print(paste0("Working on r hub ", i, " out of ", length(nodes)))
  df <- relfinal[from2 == nodes[i]]
  list <- unique(df$to2)
  df <- relfinal[to2 == nodes[i]]
  list <- append(list, unique(df$from2))
  pos <- enrichGO(gene = list, universe = names(geneList),
      keyType = "SYMBOL", OrgDb = org.Hs.eg.db, ont = "BP",
      pAdjustMethod = "BH", pvalueCutoff  = 0.01,
      qvalueCutoff  = 0.05, readable = FALSE)
  if (nrow(pos) != 0) {
    pos2 <- simplify(pos, cutoff = 0.7, by = "p.adjust",
            select_fun = min, measure = "Wang", semData = NULL)
    pos3 <- as.data.frame(pos2@result)
    pos3$regulator <- nodes[i]
    dffinal <- rbind(dffinal, pos3)
  } else {
    pos3 <- append(rep("NA", n= 9), nodes[i])
    dffinal <- rbind(dffinal, pos3)
  }
}

write.csv(dffinal, "/data/network/GOEA_NodesHubsBNRAllGOSimplify")
