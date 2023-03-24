# Differential expression analysis and Gene Ontology enrichment analysis

The same environment created for processing samples can be used here?

## Differential expression analysis
This step was performed for each sample separately considering Relapse - Diagnosis and MRD - Diagnosis.\
We used logistic regression on normalised counts data, as indicated by [best practices](https://pubmed.ncbi.nlm.nih.gov/31217225/), in each subpopulation separately. We considered only genes expressed by at least 1% of all cells.\
The parameters min.pct = 0 and minFC = 0 includes all genes in the output so we can calculate the FDR. The corrected p-value given by seurat is obtained by bonferroni method which may insert many false negatives.\
The final list of DEGs included those with at least 20% change in expression (log2FC = 0.25), FDR < 0.05, percentage of cells expressing the gene is higher than 25% in at least one timepoint.  

```
Rscript DifferentialExpression.R
```
Note: I compared multiple differential expression methods before choosing the Seurat - LR. 
In summary, Seurat gave me top DEG with high expression in one timepoint and with its expression in a high number of cells. The LR results are very similar to the wilcoxon (default option), but it also allows to include potential batch effects. As each sample was analysed separately, I didn't include the batch effect option.\
On the other hand, methods such as Limma and DESeq2 returned top DEG with low expression in few cells. 

## Gene Ontology enrichment analysis

ClusterProfiler packages was used to obtain the GOEA for genes up-regulated at relapse in each subpopulation.\
The presence of ribosomal proteins among the input list lead to all subpopulations having the same enriched pathways, all related to RNA, catabolic processes, translation and etc. The exclusion of ribosomal genes from the input list resulted in more diversified terms, kept the presence of few catabolic and translation related processes, and showed specific terms for each subpopulation.\
We used the revigo R package to reduce the redundancy among multiple GO terms and facilitate the search for meaningful answers.
```
Rscript GOEA.R
Rscript Revigo.R
```
## Prioritizing GO terms with subpopulation weight
Some subpopulations are minor to relapse and, as our goal is to find important GO terms related to relapse, we inserted a weight to prioritize responses observed in supopulations with abundant relapse cells. The weight was determined as the percentage of cells at relapse in each subpopulation. If a GO term shows up enriched for multiple subpopulations, its weight would be the sum of the subpopulations weights.

We also performed DE analysis, GOEA and weight comparing MRD with diagnosis. In this case, we observed the up-regulated GO terms at MRD and if those would persist during relapse or were transiently up-regulated.

You can observe this analysis and plots in:
```
Rscript WeightsGOEA.R
```

## Are the pathways up-regulated at relapse the drug resistance mechanisms?
We hypothesised that drug resistant mechanisms should follow two principles:\
1 - Cells with high expression of those pathways should be present at diagnosis and become more abundant at relapse.\
2 - The up-regulation of those pathways should be driven by immature myeloid cells.

The following codes tested these hypothesis:
```
Rscript DrugREsistancePotential.R
Rscript CellTypesDrugResistance.R
```
