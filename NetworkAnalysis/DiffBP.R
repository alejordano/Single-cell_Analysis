# Load libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(wordcloud)
library(tm)
library(simplifyEnrichment)
library(ggtext)
library(ggpubr)

setwd("data/network/")

cdg <- read.csv("GOEA_NodesAllSubpopHubsBNDgAllGOSimplify",
        row.names = 1) %>% mutate(timepoint = "Dg")
cr <- read.csv("GOEA_NodesAllSubpopHubsBNRAllGOSimplify",
      row.names = 1) %>% mutate(timepoint = "R")

#hubs list
hubs <- c("SCO1", "ATPAF2", "SRFBP1", "SLC7A6", "NNT", "COIL", "M6PR", "CDIP1", "EMID1")

#bn list
bn <- (unique(append(unique(cdg$regulator), unique(cr$regulator))))
bn <- bn[(bn %in% hubs) == FALSE]
bn <- append(bn, "SVBP")

#### HUBS ###
cdg2 <- cdg %>% filter(regulator %in% hubs)
cr2 <- cr %>% filter(regulator %in% hubs)

##Get clusters
#reference
clusters <- read.csv("ClusterManuallyConnections.csv")

#clustering (get those biological processes that are only in one timepoint per regulator)
data <- rbind(cdg2, cr2)
data2 <- data %>%
        group_by(regulator, Description) %>%
        filter(n_distinct(timepoint) == 1)
data3 <- data2[, c(1, 2, 10, 11)]

dataCluster <- left_join(data3, clusters, by = "Description")
dataCluster3 <- dataCluster %>% filter(is.na(cluster)) # 81 GO terms not inserted in any of the 14 clusters created. Do it manually.
#write.csv(dataCluster, "./HubsAllSubpopClusteringCuration.csv")

data4 <- read.csv("./HubsAllSubpopClusteringCuration.csv")
data4 <- data4[, -1]

# Proportions
data4 <- data4 %>%
        group_by(ClusterNumber) %>%
        mutate(ClusterSize = n_distinct(Description))
data4 <- data4 %>%
        group_by(regulator, ClusterNumber, timepoint) %>%
        mutate(Prop = n_distinct(Description) * 100 / ClusterSize)

#determine Nreg and median of the percentage
data4 <- data4 %>%
        group_by(ClusterNumber, timepoint) %>%
        mutate(NReg = n_distinct(regulator), Median = median(Prop), average = mean(Prop))

#Plot Figure 5 B and C
p2 <- ggplot(data4, aes(x = timepoint, y = Prop, fill = timepoint)) +
  geom_boxplot() +
  facet_wrap(~ cluster) + theme_bw(16) +
  theme(legend.position = "none", strip.text = element_text(size = 8)) +
  ylim(0, 100) + xlab("") + ylab("Percentages") +
  stat_compare_means(comparisons = list(c("Dg", "R")), aes(label = ..p.signif..), label.y = 80)

data4$ClusterLevels <- factor(data4$cluster, levels = c(sort(unique(data4$cluster), decreasing = TRUE)))
My_Theme = theme(axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10))

p3 <- ggplot(data4, aes(x = timepoint, y = ClusterLevels, fill = Prop)) +
      geom_tile() +
      facet_wrap(~ regulator) + theme_classic() + My_Theme + xlab("") + ylab("Cluster") +
      guides(fill = guide_legend(title = "Percentages(%))"))

figure <- ggarrange(NULL, p2, p3, labels = c("A", "B", "C"), ncol = 3, nrow = 1)

#### BN ####
cdg2 <- cdg %>% filter(regulator %in% bn)
cr2 <- cr %>% filter(regulator %in% bn)

##Get clusters
#reference
clusters

#clustering (get those biological processes that are only in one timepoint per regulator)
data <- rbind(cdg2, cr2)
data2 <- data %>%
        group_by(regulator, Description) %>%
        filter(n_distinct(timepoint) == 1)
data3 <- data2[, c(1, 2, 10, 11)]

dataCluster <- left_join(data3, clusters, by = "Description")
dataCluster2 <- dataCluster %>% filter(is.na(cluster)) # 223 GO terms not inserted in any of the 14 clusters created. Do it manually.
#write.csv(dataCluster, "./BNAllSubpopClusteringCuration.csv")

data4 <- read.csv("./BNAllSubpopClusteringCuration.csv")
data4 <- data4[, -1]

### Proportions
data4 <- data4 %>%
        group_by(ClusterNumber) %>%
        mutate(ClusterSize = n_distinct(Description))
data4 <- data4 %>%
        group_by(regulator, ClusterNumber, timepoint) %>%
        mutate(Prop = n_distinct(Description) * 100 / ClusterSize)

#determine Nreg and median of the percentage
data4 <- data4 %>%
        group_by(ClusterNumber, timepoint) %>%
        mutate(NReg = n_distinct(regulator), Median = median(Prop), average = mean(Prop))

#Plot Figure 6 B and C
p2 <- ggplot(data4, aes(x = timepoint, y = Prop, fill = timepoint)) +
      geom_boxplot() +
      facet_wrap(~ cluster) + theme_bw(16) +
      theme(legend.position = "none", strip.text = element_text(size = 8)) +
      ylim(0, 100) + xlab("") + ylab("Percentages") +
      stat_compare_means(comparisons = list(c("Dg", "R")), aes(label = ..p.signif..), label.y = 80)

data4$ClusterLevels <- factor(data4$cluster, levels = c(sort(unique(data4$cluster), decreasing = TRUE)))
My_Theme = theme(axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 7))

p3 <- ggplot(data4, aes(x = timepoint, y = ClusterLevels, fill = Prop)) +
      geom_tile() + facet_wrap(~ regulator) +
      theme_classic() + My_Theme + xlab("") + ylab("Cluster") +
      guides(fill = guide_legend(title = "Percentages(%))"))

figure2 <- ggarrange(NULL, p2, p3, labels = c("A", "B", "C"), ncol = 3, nrow = 1)