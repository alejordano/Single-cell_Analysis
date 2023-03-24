#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(patchwork))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(forcats))

path <- "/data/DEG/DEGinput/"
path2 <- "/data/DEG/DEGresults/"

## getting list of all DEGs to calculate FDR
lis <- list.files(path, pattern = "NormSubsetCluster.Rds")

for (i in seq_along(lis)){
  data <- readRDS(file = paste0(path, lis[i]))
  adata <- as.Seurat(data, counts = "counts", data = "logcounts",
          assay = "RNA", project = "SingleCellExperiment")
  perc <- ncol(adata) / 100
  expr <- as.matrix(GetAssayData(object = adata, slot = "data"))
  good <- rowSums(expr != 0)
  good <- which(good > perc)
  adata <- adata[good, ]

  patient <- sapply(strsplit(lis[i], split = "_", fixed = TRUE), function(x) (x[1]))
  tps <- c(paste0(patient, "_Dg"), paste0(patient, "_R"))
  y <- 0
  degexp <- c()
  data1 <- subset(adata, batch %in% tps)
  clstr <- unique(data1$clusters)

  for (y in seq_along(clstr)) {
    R <- c()
    d0 <- c()
    dat <- subset(data1, clusters == clstr[y])
    cond <- unique(dat$batch)
    print(paste0("working on ", patient, "_", clstr[y]))
    if (length(cond) == 2) {
      R <- ncol(subset(dat, batch == tps[2]))
      d0 <- ncol(subset(dat, batch == tps[1]))
      if ((R < 3) || (d0 < 3)) {
        print(paste0(patient, clstr[y], " had ", R, " samples at R and ", d0, " samples at Dg"))
      } else {
      Idents(dat) <- "batch"
      degexp <- FindMarkers(dat, ident.1 = tps[2], ident.2 = tps[1], verbose = FALSE, test.use = "LR", logfc.threshold = 0, min.pct = 0)
      write.csv(degexp, paste0(path2, patient, "_", clstr[y], "_Rd0LR"))
      }
    }
  }
}

list <- list.files(pattern = "Rd0LR")
tableUp <- c()

for (i in seq_along(list)) {
  patient <- c()
  cellType <- c()
  df <- c()
  pat <- c()
  ct <- c()
  dff <- c()

  #read DEG result file, calculate FDR and separate up regulated genes
  deg <- read.csv(paste0(path2, list[i]), header = TRUE, sep = ",")
  pval <- as.numeric(deg$p_val)
  fdr <- p.adjust(pval, method = "fdr", n = length(pval))
  deg <- cbind(deg, fdr)
  df <- subset(deg, fdr < 0.05 & avg_log2FC > 0.25 & (pct.1 > 0.25 | pct.2 > 0.25)) #test performed as R-0. UP regulated genes have R > 0.
  df2 <- subset(deg, pct.1 > 0.25 | pct.2 > 0.25) #need df2 for plot
  ct <- sapply(strsplit(list[i], split = "_", fixed = TRUE), function(x) (x[2]))
  pat <- sapply(strsplit(list[i], split = "_", fixed = TRUE), function(x) (x[1]))
  patient <- rep(pat, times = nrow(df))
  cluster <- rep(ct, times = nrow(df))
  dff <- cbind(patient, cluster, df)
  dff2 <- cbind(patient, cluster, df2)
  tableUp <- rbind(tableUp, dff)
  write.csv(tableUp, paste0(path2, "TableAllUpFDR"))
}

## Plot Figure 2A
tableUp$cluster <- as.numeric(tableUp$cluster)
pl <- ggplot(tableUp2 %>%
      filter(patient == "Patient1") %>%
      mutate(lab = fct_reorder(paste("Cluster", cluster), cluster)),
      aes(x = avg_logFC, y = -log10(fdr), col = (fdr <= 0.05 & abs(avg_logFC) >= 0.25))) +
      geom_point(show.legend = FALSE, size = 0.25) +
      scale_shape_manual(values = c(20, 19)) +
      scale_color_manual(values = c("lightblue", "blueviolet")) +
      theme_minimal() +
      theme(strip.background = element_rect(colour = "black")) +
      facet_wrap(~lab, nrow = 6) +
      ggtitle("Patient 1 - R x Dg") + xlab("log Fold Change") + ylab("-log10 FDR")
pl