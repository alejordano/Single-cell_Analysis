# Differential expression analysis and Gene Ontology enrichment analysis

The same environment created for processing samples can be used here?

## Differential expression analysis
We used logistic regression on normalised counts data, as indicated by [best practices](https://pubmed.ncbi.nlm.nih.gov/31217225/). We considered only genes expressed by at least 1% of all cells.\
The parameters min.pct = 0 and minFC = 0 includes all genes in the output so we can calculate the FDR. The corrected p-value given by seurat is obtained by bonferroni method which may insert many false negatives.\

```
Rscript DifferentialExpression.R
```
Note: I compared multiple differential expression methods before choosing the Seurat - LR. 
In summary, Seurat gave me top DEG with high expression in one timepoint and with its expression in a high number of cells. The LR results are very similar to the wilcoxon (default option), but it also allows to include potential batch effects. As each sample was analysed separately, I didn't include the batch effect option.\
On the other hand, methods such as Limma and DESeq2 returned top DEG with low expression in few cells. 

## Gene Ontology enrichment analysis

ClusterProfiler packages was used to obtain the GOEA for genes up-regulated at relapse in each subpopulation.\
The presence of ribosomal proteins among the input list lead to all subpopulations having the same enriched pathways, all related to RNA, catabolic processes, translation and etc. The exclusion of ribosomal genes from the input list resulted in more diversified terms, kept the presence of few catabolic and translation related processes, and showed specific terms for each subpopulation.\

```
Rscript GOEA.R
```

## Are the pathways up-regulated related to survival?
We hypothesised that survival mechanisms mechanisms should follow two principles:\
1 - Cells with high expression of those pathways should be present at the first timepoint and become more abundant later.\
2 - The up-regulation of those pathways should be driven by immature cells.

The following codes tested these hypothesis:
```
Rscript DrugResistancePotential.R
Rscript CellTypesDrugResistance.R
```
