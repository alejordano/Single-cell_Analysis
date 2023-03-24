library(dplyr) 
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

setwd("~/Documents/Iacono/") 

## Get number of interactions for each node specific from Dg and Relapse, in and out.
final <- fread("/data/aconforte/edgeAnalysis/Results/tableInteractionsAllPaired")

dgfinal <- final[(timepoint == "Dg"), ]
relfinal <- final[(timepoint == "R"), ]

#Diagnosis
dgfinal3 <- dgfinal[, .(Ninterations = uniqueN(to)), by = c("from", "sample")]
dgfinal4 <- dgfinal[, .(Ninterations = uniqueN(from)), by = c("to", "sample")]

colnames(dgfinal3) <- c("Gene", "Sample", "N_Out")
colnames(dgfinal4) <- c("Gene", "Sample", "N_In")

dgfinal5 <- full_join(dgfinal3, dgfinal4, by = c("Gene", "Sample"))
dgfinal5$N_In <- ifelse(is.na(dgfinal5$N_In), 0, dgfinal5$N_In)
dgfinal5$N_Out <- ifelse(is.na(dgfinal5$N_Out), 0, dgfinal5$N_Out)
#option: dgfinal7 <- setnafill(dgfinal5, nan = NA, fill = 0, cols = c("N_Out", "N_In"))
dgfinal5$connections <- dgfinal5$N_Out + dgfinal5$N_In
dgfinal6 <- dgfinal5[, .(NSamples = uniqueN(Sample)), by = c("Gene")]

## Get Nodes present in all subpopulations
#hist(dgfinal6$NSamples)
hfn <- dgfinal6[dgfinal6$NSamples == 49]
hfn <- unique(hfn$Gene)

##Get the number of interactions for each gene
dgfinal8 <- dgfinal5 %>% filter(Gene %in% hfn)

##plot the median N_out and N_In - select hubs
dgfinal9 <- dgfinal8 %>% group_by(Gene) %>% mutate(MedianReg = median(connections))
dgfinal9 <- unique(dgfinal9[,c(1, 6)])

## filter nodes into connections with 95% confidence.
dghubs <- dgfinal9 %>% filter(MedianReg > quantile(dgfinal9$MedianReg, .95))

##Get activity of target genes
geneList <- unique(c(dgfinal$from, dgfinal$to))

##hubs
dghubslist <- unique(dghubs$Gene)
dghubslist <- str_replace(dghubslist, "^MT[.]", "")
dghubslist <- str_replace(dghubslist, "[.]AS1", "")

dgfinal$from2 <- str_replace(dgfinal$from, "^MT[.]", "")
dgfinal$from2 <- str_replace(dgfinal$from2, "[.]AS1", "")
dgfinal$to2 <- str_replace(dgfinal$to, "^MT[.]", "")
dgfinal$to2 <- str_replace(dgfinal$to2, "[.]AS1", "")


dffinal <- c()
for (i in seq_along(dghubslist)){
  print(paste0("Working on dg hub ", i, " out of ", length(dghubslist)))
  df <- dgfinal[from2 == dghubslist[i]]
  list<- unique(df$to2)
  df <- dgfinal[to2 == dghubslist[i]]
  list <- append(list, unique(df$from2))
  pos <- enrichGO(gene = list, universe = names(geneList), keyType = "SYMBOL", OrgDb = org.Hs.eg.db, ont = "BP", 
                  pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = FALSE) #readable needs to be false to run with "SYMBOL"
  if (!is.null(pos)){
    pos2 <- simplify(pos, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
    pos3 <- as.data.frame(pos2@result)
    pos3$regulator <- dghubslist[i]
    dffinal <- rbind(dffinal, pos3)
  }else{
    pos3 <- cbind(rep("NA", n= 9), dghubslist[i])
    dffinal <- rbind(dffinal, pos3)
  }
  
}

write.csv(dffinal, "/data/aconforte/edgeAnalysis/GOEA_ConnectionsDgAllGOSimplify")

### do everything again for relapse
dgfinal3 <- relfinal[, .(Ninterations = uniqueN(to)), by = c("from", "sample")]
dgfinal4 <- relfinal[, .(Ninterations = uniqueN(from)), by = c("to", "sample")]

colnames(dgfinal3) <- c("Gene", "Sample", "N_Out")
colnames(dgfinal4) <- c("Gene", "Sample", "N_In")

dgfinal5 <- full_join(dgfinal3, dgfinal4, by = c("Gene", "Sample"))
dgfinal5 <- setnafill(dgfinal5, nan = NA, fill = 0, cols = c("N_Out", "N_In"))
dgfinal5$connections <- dgfinal5$N_Out + dgfinal5$N_In
dgfinal6 <- dgfinal5[, .(NSamples = uniqueN(Sample)), by = c("Gene")]

## Get Nodes present in all subpopulations
hfnrel <- dgfinal6[dgfinal6$NSamples == 49]
hfnrel <- unique(hfnrel$Gene)

##Get the number of interactions for each gene
dgfinal8 <- dgfinal5 %>% filter(Gene %in% hfnrel)

##plot the median N_out and N_In - select hubs
dgfinal9 <- dgfinal8 %>% group_by(Gene) %>% mutate(MedianReg = median(connections))
dgfinal9 <- unique(dgfinal9[,c(1,6)])

## filter nodes into In and Out with 95% confidence.
relhubs <- dgfinal9 %>% filter(MedianReg > quantile(dgfinal9$MedianReg, .95))

#### Check target list activity
geneList <- unique(c(relfinal$from, relfinal$to))

##Hubs
dghubslist <- unique(relhubs$Gene)
dghubslist <- str_replace(dghubslist, "^MT[.]", "")
dghubslist <- str_replace(dghubslist, "[.]AS1", "")

relfinal$from2 <- str_replace(relfinal$from, "^MT[.]", "")
relfinal$from2 <- str_replace(relfinal$from2, "[.]AS1", "")
relfinal$to2 <- str_replace(relfinal$to, "^MT[.]", "")
relfinal$to2 <- str_replace(relfinal$to2, "[.]AS1", "")

dffinal <- {}

for (i in seq_along(dghubslist)){
  print(paste0("Working on r hub ", i, " out of ", length(dghubslist)))
  df <- relfinal[from2 == dghubslist[i]]
  list <- unique(df$to2)
  df <- relfinal[to2 == dghubslist[i]]
  list <- append(list, unique(df$from2))
  pos <- enrichGO(gene = list, universe = names(geneList), keyType = "SYMBOL", OrgDb = org.Hs.eg.db, ont = "BP", 
                  pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = FALSE) #readable needs to be false to run with "SYMBOL"
  if (nrow(pos) != 0){
    pos2 <- simplify(pos, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
    pos3 <- as.data.frame(pos2@result)
    pos3$regulator <- dghubslist[i]
    dffinal <- rbind(dffinal, pos3)
  }else{
    pos3 <- append(rep("NA", n= 9), dghubslist[i])
    dffinal <- rbind(dffinal, pos3)
  }
}

write.csv(dffinal, "/data/aconforte/edgeAnalysis/GOEA_ConnectionsRAllGOSimplify")
