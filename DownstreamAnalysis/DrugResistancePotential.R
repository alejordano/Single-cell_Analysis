library(dplyr)
library(ggplot2)
library(Seurat)
library(reshape)
library(ggpubr)

### Are cells with high expression of Mitochondrial Respiration genes already present at diagnosis?
#get all genes related to MT respiration
library(org.Hs.eg.db)
library(GO.db)

go_id <- GOID(GOTERM[Term(GOTERM) == "oxidative phosphorylation"])
allegs <- get(go_id, org.Hs.egGO2ALLEGS)
genes <- unlist(mget(allegs, org.Hs.egSYMBOL))

go_id <- GOID(GOTERM[Term(GOTERM) == "proton transmembrane transport"])
allegs <- get(go_id, org.Hs.egGO2ALLEGS)
genes1 <- unlist(mget(allegs, org.Hs.egSYMBOL))

go_id <- GOID(GOTERM[Term(GOTERM) == "ATP synthesis coupled electron transport"])
allegs <- get(go_id, org.Hs.egGO2ALLEGS)
genes2 <- unlist(mget(allegs, org.Hs.egSYMBOL))

oxgenes <- unique(c(genes, genes1, genes2))

#get proportion of cells with high expression of DE MT genes.
deg <- read.csv("/data/Source/TableAllUpFDRRd0", row.names = 1)

d <- readRDS("./AllMerged_NormSubsetCluster.Rds")
d$pat <- sapply(strsplit(d$batch, split = "_", fixed = TRUE), function(x) (x[1]))

Diag <- c()
Rem <- c()
Rel <- c()
Pati <- c()
Gene <- c()

for (w in seq_along(pats)) {
  deg1 <- deg %>% filter(patient == pats[w])
  list <- unique(deg1$X)
  listf <- list[list %in% oxgenes]
  print(paste0("Working on patient ", pats[w]))

  #get oxphos genes' expression
  geneexp <- data1@assays[["RNA"]]@data
  geneexp2 <- geneexp[rownames(geneexp) %in% listf, ]
  geneexp3 <- as.data.frame(t(as.matrix(geneexp2)))

  #get dataframe with umap coordinates
  df <- read.csv(paste0("./JN/Df_", pats[w]))
  rownames(df) <- df$X

  #df final with umap plots and gene expression
  dff <- merge(df, geneexp3, by = "row.names")
  dff <- dff [, -1]
  dff$timepoint <- sapply(strsplit(dff$batch, split = "_", fixed = TRUE), function(x) (x[2]))
  lim <- length(dff)
  print(paste0("Dff fine with ", lim - 1, " genes"))

  for (i in 7:(lim - 1)) {
    dist <- as.numeric(ifelse(dff[, i] > 0, dff[, i], "na"))
    dist <- dist[!is.na(dist)]
    # determine high expression as 5% highest end of the gene expression distribution accross all timepoints
    a <- qnorm(0.95, mean = mean(dist), sd = sd(dist))
    dff[, (lim + (i - 6))] <- ifelse(dff[, i] > a, as.character("high"), as.character("other"))
    colnames(dff)[(lim + (i - 6))] <- paste0(colnames(dff)[i], "_hl")
  }
  print(paste0("Cuts determined"))

  melt_data <- melt(dff, id = colnames(dff)[1:lim])
  melt_data$value <- factor(melt_data$value, levels = c("high", "other"))
  gen <- unique(melt_data$variable)
  gtt <- colnames(melt_data)[7:lim - 1]
  print(paste0("Melt data okay"))

  for (i in seq_along(gen)) {
    print(paste0("Getting proportion of cells. ", i, "/", length(gen)))
    dff1 <- melt_data %>% filter(variable == gen[i])
    #Proportion of cells with high expression of each MT respiration gene at each timepoint.
    dff2 <- dff1 %>% filter(timepoint == "Dg")
    dg <- round((sum(dff2$value == "high") / nrow(dff2)) * 100, digits = 2)
    dff4 <- dff1 %>% filter(timepoint == "MRD")
    mrd <- round((sum(dff4$value == "high") / nrow(dff4)) * 100, digits = 2)
    dff3 <- dff1 %>% filter(timepoint == "R")
    r <- round((sum(dff3$value == "high") / nrow(dff3)) * 100, digits = 2)
    Diag <- append(Diag, dg)
    Rem <- append(Rem, mrd)
    Rel <- append(Rel, r)
    Pati <- append(Pati, pats[w])
    Gene <- append(Gene, as.character(gen[i]))

    #UMAP plot Figure 7C:
    df_layer_1 <- dff1 %>% filter(value == "other")
    df_layer_2 <- dff1 %>% filter(value == "high")
    p <- ggplot() +  geom_point(data = df_layer_1, aes(CorrectedUMAP1, CorrectedUMAP2), colour = "darkseagreen2", size = 1, alpha = 5 / 10) +
      geom_point(data = df_layer_2, aes(CorrectedUMAP1, CorrectedUMAP2), colour = "blue4", size = 1, alpha = 5 / 10) + xlab("UMAP1") + ylab("UMAP2") +
      facet_wrap(. ~ batch) + theme_classic() +
      ggtitle(paste0(gtt[i], " Dg = ", dg, " MRD = ", mrd, " R = ", r)) + theme(plot.title = element_text(size = 10))
    plot_list[[i]] <- p
  }
  jpeg(file = paste0("Results/GenesPanel_", pats[w], ".jpeg"), width = 3800, height = 2500, res = 300)
  multiplot(plotlist = plot_list, cols = 4)
  dev.off()
}

## boxplot Figure 7D:
dfinal <- as.data.frame(cbind(Dg = Diag, MRD = Rem, R = Rel, Pat = Pati, Gene = Gene))
write.csv(dfinal, "data/MTchar/BoxplotDEG")
dfinal <- dfinal %>% filter(!is.na(Dg))
dfinal$Gene <- sapply(strsplit(dfinal$Gene, split = "_", fixed = TRUE), function(x) (x[1]))
melt_dfinal <- melt(dfinal, id = colnames(dfinal)[4:5])
melt_dfinal$value <- as.numeric(melt_dfinal$value)

Pat <- c("1216", "1886", "3904", "40389", "4978", "5750", "6323", "3822", "3853", "5143")
nas <- c(paste0("Patient ", 1:10))
level <-  as.factor(c(paste0("Patient ", 1:10)))
n <- c(1:10)

ns <- as.data.frame(cbind(Pat, nas, n))
ns$Pat <- as.integer(ns$Pat)
melt_dfinal <- left_join(melt_dfinal, ns, by = "Pat")

my_comparisons <- list(c("Dg", "R"))
melt_dfinal$variable <- as.character(melt_dfinal$variable)
melt_dfinal <- as.data.frame(melt_dfinal)

pl <- ggplot(melt_dfinal, aes(x = variable, y = value, color =  variable)) + geom_boxplot() + theme_classic() +
  facet_wrap(~factor(nas, levels = level)) +
  xlab("") + ylab("Percentage of cells") + theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif", label.y = 35) + ylim(0, 40)
