suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(patchwork))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(clusterProfiler))
organism = "org.Hs.eg.db"
suppressPackageStartupMessages(require(organism, character.only = TRUE))

path <- "/data/DEG/DEGresults/"
path2 <- "/data/GOEA/"

tableUp <- read.csv(paste0(path, "RDg_tableAllUpFDR"), header = TRUE, sep = ",")

#background genes#
data <- readRDS(paste0("./", sample[i], "_NormSubsetCluster.Rds"))
adata <- as.Seurat(data, counts = "counts", data = "logcounts", assay = "RNA", project = "SingleCellExperiment")

resultado <- c()
sample <- unique(tableUp$patient)

for (i in seq_along(sample)) {
patient <- c()
df <- c()
pat <- c()
ct <- c()
dff <- c()
resultadoup <- c()
resultadoUp <- c()

df1 <- subset(tableUp, patient == sample[i])
ct <- unique(df1$cluster)

    for (y in seq_along(ct)) {
        print(paste0("working on ", sample[i], "_", ct[y]))
        df <- subset(df1, cluster == ct[y])
        genes <- unique(df$X)
        genes <- sort(genes, decreasing = FALSE)
        ribo.genes <- grep(pattern = "^RP[SL]", genes, value = TRUE)
        degs <- genes[!(genes %in% ribo.genes)]

        #background genes#
        adata2 <- subset(adata, clusters == ct[y] & batch != paste0(sample[i], "_MRD"))
        expr <- as.matrix(GetAssayData(object = adata2, slot = "data"))
        bg_names <- names(rowSums(expr)[rowSums(expr) > 0])
        if (length(genes) < 3) {
            print(paste0(sample[i], " ", ct[y], "had only ", length(genes), "DEGs"))
        } else {
        ENRUp <- enrichGO(degs,
                ont = "BP",
                keyType = "SYMBOL",
                minGSSize = 3,
                maxGSSize = 800,
                pvalueCutoff = 0.05,
                OrgDb = organism,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                universe = bg_names)
            if (!is.na(ENRUp[1, 1])) {
                go <- head(ENRUp$ID, n = 30)
                des <- head(ENRUp$Description, n = 30)
                ratio <- head(ENRUp$GeneRatio, n = 30)
                ratio2 <- sapply(strsplit(ratio, split = "/", fixed = TRUE), function(x) (x[1]))
                total <- sapply(strsplit(ratio, split = "/", fixed = TRUE), function(x) (x[2]))
                GeneRatio <- as.numeric(ratio2)
                genes <- head(ENRUp$geneID, n = 30)
                padj <- head(ENRUp$p.adjust, n = 30)
                resultadoup <- cbind(go, des, ratio, GeneRatio, genes, padj)
                resultadoup <- as.data.frame(resultadoup)
                #I think the 10 lines above could be reduced to:
                #resultadoup <- head(ENRUp@result, n = 30)
                patient <- rep(sample[i], times = nrow(resultadoup))
                cluster <- rep(ct[y], times = nrow(resultadoup))
                resultadof <- cbind(resultadoup, patient, cluster)
                resultado <- rbind(resultado, resultadof)
            }
        }
    }
}
write.csv(resultado, paste0(path2, "CharAllUpGOFilterRPs"))