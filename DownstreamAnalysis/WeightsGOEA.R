library(dplyr)
library(ggplot2)

path <- "/data/GOEA/"

# determine proportion of cells in each cluster and timepoint
table <- read.csv("./data/Source/subclusteringsummary.csv")

pats <- unique(table$pat)

tablefinal <- c()
for (i in seq_along(pats)) {
  table2 <- table %>% filter(pat == pats[i])
  total <- nrow(table2)
  cntclst <- table2 %>% group_by(clusters) %>% summarise(n = n())
  cntclsttp <- table2 %>% group_by(clusters, batch) %>% count()
  table3 <- merge(cntclst, cntclsttp, by = "clusters")
  table3 <- table3 %>% mutate(batchperc = (n.y / total) * 100, patient = pats[i])
  tablefinal <- rbind(tablefinal, table3)
}

#### Relapse analysis ####
tablef <- tablefinal %>% filter(batch == "R") %>% mutate(sample = paste0(patient, "_", clusters))

# upload revigo data.
data <- read.csv(paste0(path, "rrvgoResult_CharAllUpFilterRPsRd0"))

data <- data[, -c(1, 5, 8)]

dat <- data %>% mutate(sample = paste0(patient, "_", cluster))
mydat <- cbind(term = dat$parentTerm, sample = dat$sample, patient = dat$patient)
mydat <- as.data.frame(mydat)

df <- mydat[!duplicated(mydat), ]

## insert weigths on df
tab <- as.data.frame(cbind(batchperc = tablef$batchperc, sample = tablef$sample))
dff <- merge(df, tab)
dff$batchperc <- as.numeric(dff$batch)
scores <- dff %>%
          group_by(term) %>%
          mutate(weights = sum(batchperc), nClusters = n_distinct(sample),
          nPatient = n_distinct(patient))

write.csv(scores, paste0(path, "summaryrrvgoweigthedfreqRd0"))

### Plot Figure 2 B and C
# highest weights
data <- scores
topw <- head(unique(sort(data$weights, decreasing = TRUE)), n = 10)
dfpl <- data %>% filter(weights %in% topw)
a <- ggplot(dfpl, aes(y = reorder(term, weights), x = weights, color = nPatient, size = nClusters)) +
  geom_point() + scale_color_gradient(low = "blue", high = "red") +
  ggtitle("GOEA") + xlab("Weights") +
  ylab("") + theme(text = element_text(size = 15))

# highest number of patients
topw <- head(unique(sort(data$nPatient, decreasing = TRUE)), n = 3)
dfpl <- data %>% filter(nPatient %in% topw)
b <- ggplot(dfpl, aes(y = reorder(term, nPatient), x = nPatient, colour = nPatient, size = nClusters)) +
  geom_point() + scale_color_gradient(low = "blue", high = "red") +
  ggtitle("GOEA") + xlab("Number of Patient") +
  ylab("") + theme(text = element_text(size = 15))

#Plot figure 2:
panel <- plot_grid(pl, a, b, ncol = 3, nrow = 1, labels = c("A", "B", "C"))
ggsave(file = "Figure2.png", panel, bg = "white", height = 15, width = 60, units = "cm")

# different profile for different patients (Supplementary figure 3):
pat <- unique(data$patient)
plot_list <- list()
for (i in seq_along(pat)) {
dfpl <- data %>% filter(patient == pat[i])
topw <- head(unique(sort(dfpl$weights, decreasing = TRUE)), n = 10)
dfpl2 <- dfpl %>% filter(weights %in% topw)
c <- ggplot(dfpl2, aes(y = reorder(term, weights), x = weights, color = nPatient, size = nClusters)) +
  geom_point() + scale_color_gradient(low = "blue", high = "red") +
  ggtitle(paste0("GOEA - Patient ",  i)) + xlab("Weights") +
  ylab("") + theme(text = element_text(size = 10))
plot_list[[i]] = c
}
ggsave(file = "SupplementaryFigure3.png",
      arrangeGrob(grobs = plot_list, ncol = 3),
      height = 50, width = 50, units = "cm")

#### MRD Analysis ####
tablef <- tablefinal %>% filter(batch == "MRD") %>% mutate(sample = paste0(patient, "_", clusters))

### Get revigo data.
data2 <- read.csv(paste0(path, "rrvgoResult_CharAllUpFilterRPsMRDd0"))
data2 <- data2[, -c(1, 5, 8)]

dat2 <- data2 %>% mutate(sample = paste0(patient, "_", cluster))
mydat <- cbind(term = dat2$parentTerm, sample = dat2$sample, patient = dat2$patient)
mydat <- as.data.frame(mydat)

df2 <- mydat[!duplicated(mydat), ]

## insert weigths on df
tab2 <- as.data.frame(cbind(batchperc = tablef$batchperc, sample = tablef$sample))
dff2 <- merge(df, tab2)
dff2$batchperc <- as.numeric(dff2$batch)
scores2 <- dff2 %>%
          group_by(term) %>%
          mutate(weights = sum(batchperc), nClusters = n_distinct(sample),
          nPatient = n_distinct(patient))

write.csv(scores2, paste0(path, "summaryrrvgoweigthedfreqMRDd0"))

## Plot Figure 3:
#Highest weight
#dataMRD <- read.csv("summaryrrvgoweigthedfreqMRDd0", row.names=1)
dataMRD <- scores2

# highest weights
topw <- head(unique(sort(dataMRD$weights, decreasing = TRUE)), n = 10)
dfpl <- dataMRD %>% filter(weights %in% topw)
a <- ggplot(dfpl, aes(y = reorder(term, weights), x = weights, color = nPatient, size = nClusters)) +
  geom_point() + scale_color_gradient(low = "blue", high = "red") +
  ggtitle("GOEA") + xlab("Weights") +
  ylab("") + theme(text = element_text(size = 15))
a

# With data from relapse. Are those changes also present at relapse or are they transient?
dfpl2 <- left_join(dfpl, data, by = c("sample", "term"))
colnames(dfpl2)[12] <- "nPatient_Relapse"
colnames(dfpl2)[11] <- "nClusters_Relapse"
b <- ggplot(dfpl2, aes(y = reorder(term, weights.x), x = weights.y, color = nPatient_Relapse, size = nClusters_Relapse)) +
  geom_point() + scale_color_gradient(low = "blue", high = "red") +
  ggtitle("GOEA") + xlab("Weights") +
  ylab("") + theme(text = element_text(size = 15))
b

#Get DEGs related to cytokine
library(org.Hs.eg.db)
library(GO.db)

go_id <- GOID(GOTERM[Term(GOTERM) == "cytokine-mediated signaling pathway"])
allegs <- get(go_id, org.Hs.egGO2ALLEGS)
genes <- unlist(mget(allegs, org.Hs.egSYMBOL))
genes

go_id <- GOID(GOTERM[Term(GOTERM) == "cytokine production"])
allegs <- get(go_id, org.Hs.egGO2ALLEGS)
genes1 <- unlist(mget(allegs, org.Hs.egSYMBOL))
genes1

ctkgenes <- unique(c(genes, genes1))

#get proportion of cells with high expression of DE cytokine genes.
deg <- read.csv("/data/TableAllUpFDRMRDd0", row.names = 1)

deg1 <- deg %>% filter(X %in% ctkgenes)
deg1$sample <- paste0(deg1$patient, "_", deg1$cluster)
df <- deg1 %>%
      group_by(X) %>%
      mutate(NPat = n_distinct(patient), NCluster = n_distinct(sample))
df <- as.data.frame(df)
df1 <- unique(df[, c(3, 11:12)])
topw <- head(unique(sort(df$NCluster, decreasing = TRUE)), n = 5)
dfpl3 <- df1 %>% filter(NCluster %in% topw)
dfpl3 <- as.data.frame(dfpl3)
c <- ggplot(dfpl3, aes(y = reorder(X, NPat), x = NPat, size = NCluster)) +
  geom_point() +
  ggtitle("DEGs") + xlab("Number of Patients") +
  ylab("") + theme(text = element_text(size = 15))
c

##panel
panel2 <- plot_grid(a, b, c, ncol = 3, nrow = 1, labels = c("A", "B", "C"), label_size = 20)
ggsave(file = "Figure5.png", panel2, bg = "white", height = 15, width = 70, units = "cm")