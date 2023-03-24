library(ggplot2)
library(dplyr)
library(ggpubr)
library(reshape2)
library(patchwork)

setwd("/data/network")

table <- read.csv("AllSubpopAllNodeMeasures.csv", row.names = 1)

#Restructure data so centrality measures of nodes at Dg and R are columns
tableR <- table %>% filter(timepoint == "R")
tableDg <- table %>% filter(timepoint == "Dg")

#Create filtering to get changes for each node at each subpopulation.
tableR$filtering <- paste0(tableR$gene, "_", tableR$sample)
tableDg$filtering <- paste0(tableDg$gene, "_", tableDg$sample)
tableR2 <- tableR %>% filter(filtering %in% tableDg$filtering)
tableDg2 <- tableDg %>% filter(filtering %in% tableR$filtering)

###Differences R - Dg for each measurement
Dif <- merge(tableR2, tableDg2, by = c("sample", "gene"))
Dif$degree <- Dif$degree.x - Dif$degree.y
Dif$betwe <- Dif$betwe.x - Dif$betwe.y
Dif$close <- Dif$close.x - Dif$close.y

#get significant differences: 25% as thresholds to obtain transition between quantiles
dif2 <- Dif %>%
  group_by(gene) %>%
  mutate(updegree = (sum(degree > 0.25) / n_distinct(sample)) * 100,
         downdegree = (sum(degree < -0.25) / n_distinct(sample)) * 100,
         upbetwee = (sum(betwe > 0.25) / n_distinct(sample)) * 100,
         downbetwee = (sum(betwe < -0.25) / n_distinct(sample)) * 100,
         upclose = (sum(close > 0.25) / n_distinct(sample)) * 100,
         downclose = (sum(close < -0.25) / n_distinct(sample)) * 100,
         percSubopop = (n_distinct(sample) / 54) * 100)

dif2 <- as.data.frame(dif2)
#write.csv(dif2, "./MeasuresDifferencesAllNodesAllSubpop.csv")

dif3 <- unique(dif2[, c(2, 18:24)])
dif4 <- dif3[order(-dif3$percSubopop, -dif3$updegree, dif3$downdegree), ]

#filter for those present in more than 25% of the subpopulations
dif4 <- dif3[(dif3$percSubopop > 25), ]

###Betwenness analysis:
dif5 <- dif4 %>% filter(dif4$upbetwee > dif4$downbetwee) # change signal here to determine increase/decrease
dif5$difbetwee <- dif5$upbetwee - dif5$downbetwee
dif5 <- dif5[order(-dif5$difbetwee), ] # order by difference or % of populations that increased/decreased centrality
pdeg <- head(dif5, n = 50)
pdeg <- pdeg[, c(1, 4, 5, 8:9)]
pdeg$gene <- factor(pdeg$gene, levels = pdeg$gene)

#melt data frame into long format
df <- melt(pdeg, id.vars = c("gene", "percSubopop", "difbetwee"), variable.name = "series")
df2 <- df %>% filter(value < 15 & series == "downbetwee")
df3 <- df %>% filter(value > 45 & series == "upbetwee")

list <- intersect(df2$gene, df3$gene)
df4 <- df %>% filter(df$gene %in% list)

#Plot figure 6 A
betp <- ggplot(df4, aes(forcats::fct_rev(gene), value, fill = series)) +
        geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1) +
        ylab("Percentage of subpopulations") + xlab("") + ggtitle("Betweenness") +
        theme_linedraw() +
        theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length = unit(0, "pt")) +
        scale_fill_discrete(name = "", labels = c("Increased", "Decreased")) + coord_flip()

###Degree analysis:
dif5 <- dif4 %>% filter(dif4$updegree > dif4$downdegree) # change signal here to determine increase/decrease
dif5$difdegree <- dif5$updegree - dif5$downdegree
dif5 <- dif5[order(-dif5$difdegree), ] # order by difference or % of populations that increased/decreased centrality
pdeg <- head(dif5, n = 50)
pdeg <- pdeg[-8, c(1:3, 8:9)]
pdeg$gene.x <- factor(pdeg$gene, levels = pdeg$gene)

#melt data frame into long format
df <- melt(pdeg,  id.vars = c("gene", "percSubopop", "difdegree"), variable.name = "series")
df2 <- df %>% filter(value < 15 & series == "downdegree")
df3 <- df %>% filter(value > 50 & series == "updegree")

list <- intersect(df2$gene.x, df3$gene.x)
df4 <- df %>% filter(df$gene.x %in% list)

#Figure 5A
degp <- ggplot(df4, aes(forcats::fct_rev(gene), value, fill = series)) + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1) +
  ylab("Percentage of subpopulations") + xlab("") + ggtitle("Degree") + theme_linedraw() +
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length = unit(0, "pt")) +
  scale_fill_discrete(name = "", labels = c("Increased", "Decreased")) + coord_flip()
