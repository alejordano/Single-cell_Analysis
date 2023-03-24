library(bigSCale)
library(SeuratObject)
library(stringr)

path2 <- "~/Documents/Iacono/"
md <- readRDS(paste0(path2, "data/AllMerged_NormSubsetCluster.Rds"))

#make loop and let it
md$patient <- sapply(strsplit(md$batch, split = "_", fixed = TRUE), function(x) (x[1]))
md$timepoint <- sapply(strsplit(md$batch, split = "_", fixed = TRUE), function(x) (x[2]))
pat <- unique(md$patient)
for (i in seq_along(pat)) {
    md1 <- subset(md, patient == pat[i])
    cltrs <- unique(md1$clusters)
    for (y in seq_along(cltrs)){
        #relapse
        md2 <- subset(md1, clusters == cltrs[y])
        test <- table(md2$timepoint)["R"]
        if ((!is.na(test) && test > 50) == TRUE) {
            md2 <- subset(md1, clusters == cltrs[y] & timepoint == "R")
            print(paste0(pat[i], "_", cltrs[y], "_R ", "number of cells: ", ncol(md2)))
            expr.rlp <- md2$RNA@counts
            gene.namesmd <- as.data.frame(rownames(md2$RNA@counts))
            gene.namesmd <- as.matrix(gene.namesmd)
            results.rlp = compute.network(expr.data = expr.rlp, gene.names = gene.namesmd, clustering = "direct")
            write.csv(results.rlp$correlations, paste0("/data/aconforte/edgeAnalysis/Results/", pat[i], "_", cltrs[y], "_R_correlations"))
            write.csv(as.data.frame(results.rlp$centrality), paste0("/data/aconforte/edgeAnalysis/Results/", pat[i], "_", cltrs[y], "_R_centrality"))
            write_graph(results.rlp$graph, paste0("/data/aconforte/edgeAnalysis/Results/", pat[i], "_", cltrs[y], "_R_graph.gml"), format = "gml")
        }

        #diagnosis
        md3 <- subset(md1, clusters == cltrs[y])
        test2 <- table(md3$timepoint)["Dg"]
        if (!is.na(test2) && test2 > 50) {
            md3 <- subset(md1, clusters == cltrs[y] & timepoint == "Dg")
            print(paste0(pat[i], "_", cltrs[y], "_Dg ", "number of cells: ", ncol(md3)))
            expr.dg <- md3$RNA@counts
            gene.namesdg <- as.data.frame(rownames(md3$RNA@counts))
            gene.namesdg <- as.matrix(gene.namesdg)
            results.dg = compute.network(expr.data = expr.dg, gene.names = gene.namesdg, clustering = "direct")
            write.csv(results.rlp$correlations, paste0("/data/aconforte/edgeAnalysis/Results/", pat[i], "_", cltrs[y], "_Dg_correlations"))
            write.csv(as.data.frame(results.dg$centrality), paste0("/data/aconforte/edgeAnalysis/Results/", pat[i], "_", cltrs[y], "_Dg_centrality"))
            write_graph(results.dg$graph, paste0("/data/aconforte/edgeAnalysis/Results/", pat[i], "_", cltrs[y], "_Dg_graph.gml"), format = "gml")
        }
    }
}