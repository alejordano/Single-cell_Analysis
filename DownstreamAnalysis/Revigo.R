library(dplyr)
library(rrvgo)

path <- "/data/GOEA/"

lista <- list.files(path, pattern = "^CharAll")

for (w in seq_along(lista)) {
  table <- read.csv(lista[w], header = TRUE, sep = ",")
  table <- table %>% filter(padj < 0.05)
  pat <- unique(table$patient)
  resultadof <- c()
  for (i in seq_along(pat)) {
    tablef <- table %>% filter(patient == pat[i])
    clst <- unique(tablef$cluster)
    for (y in seq_along(clst)) {
      tablefc <- tablef %>% filter(cluster == clst[y])
      print(paste0("working on ", lista[w], pat[i], clst[y]))
      if (nrow(tablefc) > 2) {
        simMatrix <- calculateSimMatrix(tablefc$go, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
        scores <- setNames(-log10(tablefc$padj), tablefc$go)
        reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")
        patient <- rep(pat[i], n = length(reducedTerms))
        cluster1 <- rep(clst[y], n = length(reducedTerms))
        resultado <- cbind(reducedTerms, patient, cluster1)
        resultadof <- rbind(resultadof, resultado)
      }
    }
  }
  resultadofinal <- cbind(des = resultadof$term, parentTerm = resultadof$parentTerm,
                    parentSimScore = resultadof$parentSimScore, patient = resultadof$patient,
                    cluster = resultadof$cluster1, size = resultadof$size, score = resultadof$score)
  resultadofinal <- as.data.frame(resultadofinal)
  list <- unique(resultadofinal$parentTerm)
  tablefinal <- merge(table, resultadofinal, by = c("des", "patient", "cluster"))
  write.csv(tablefinal, paste0(path, "rrvgoResult_", lista[w]))
}