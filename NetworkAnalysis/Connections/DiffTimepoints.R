library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

setwd("~/Documents/Iacono/")

final <- fread("/data/aconforte/edgeAnalysis/Results/tableInteractionsAllPaired")

dgfinal <- final[(timepoint == "Dg"), ]
relfinal <- final[(timepoint == "R"), ]

#diagnosis
dgfinal3 <- dgfinal[, .(Ninterations = uniqueN(to)), by = c("from", "sample")]
dgfinal4 <- dgfinal[, .(Ninterations = uniqueN(from)), by = c("to", "sample")]

colnames(dgfinal3) <- c("Gene", "Sample", "N_Out")
colnames(dgfinal4) <- c("Gene", "Sample", "N_In")

dgfinal5 <- full_join(dgfinal3, dgfinal4, by = c("Gene", "Sample"))
dgfinal5 <- setnafill(dgfinal5, nan = NA, fill = 0, cols = c("N_Out", "N_In"))
dgfinal5$connections <- dgfinal5$N_Out + dgfinal5$N_In
dgfinal6 <- dgfinal5[, .(Nnodes = uniqueN(Gene)), by = c("Sample")]
dg <- merge(dgfinal5, dgfinal6, by = "Sample")
dg$timepoint <- "Dg"


#Relapse
rfinal3 <- relfinal[, .(Ninterations = uniqueN(to)), by = c("from", "sample")]
rfinal4 <- relfinal[, .(Ninterations = uniqueN(from)), by = c("to", "sample")]

colnames(rfinal3) <- c("Gene", "Sample", "N_Out")
colnames(rfinal4) <- c("Gene", "Sample", "N_In")

rfinal5 <- full_join(rfinal3, rfinal4, by = c("Gene", "Sample"))
rfinal5 <- setnafill(rfinal5, nan = NA, fill = 0, cols = c("N_Out", "N_In"))
rfinal5$connections <- rfinal5$N_Out + rfinal5$N_In
rfinal6 <- rfinal5[, .(Nnodes = uniqueN(Gene)), by = c("Sample")]
r <- merge(rfinal5, rfinal6, by = "Sample")
r$timepoint <- "R"

dfplot <- rbind(dg, r)
dfplot$ratio <- dfplot$connections / dfplot$Nnodes
#write.csv(dfplot, "tableFigure4.csv")

#plot Figure 4:
plt1 <- ggplot(dfplot, aes(y = Nnodes, x = timepoint, fill = timepoint)) + geom_boxplot() + theme_bw() +
  theme(legend.position = "none") + stat_compare_means(method = "t.test", comparisons = list(c("Dg", "R")), aes(label = ..p.signif..), label.x = 1.5) + 
  xlab("") + ylab("Number of nodes")

plt2 <- ggplot(dfplot, aes(y = connections, x = timepoint, fill = timepoint)) + geom_boxplot() + theme_bw() + 
  theme(legend.position = "none") + stat_compare_means(method = "t.test", comparisons = list(c("Dg", "R")), aes(label = ..p.signif..), label.x = 1.5) + 
  xlab("") + ylab("Number of edges")

plt3 <- ggplot(dfplot, aes(y = ratio, x = timepoint, fill = timepoint)) + geom_boxplot() + theme_bw() + 
  stat_compare_means(method = "t.test", comparisons = list(c("Dg", "R")), aes(label = ..p.signif..), label.x = 1.5) + 
  xlab("") + ylab("Ratio edges/node")


panel1 <- ggarrange(plt1, plt2, plt3, ncol = 3, labels = c("A", "B", "C"))
ggsave("p4.png", width = 6, height = 4)