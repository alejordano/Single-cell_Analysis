# Processing the samples: Normalisation, Batch correction and Subclustering

The same environment created for pre-processing can be used for processing the samples.

## Normalisation 
First, samples from all three timepoint for each patient were concatenated. Then, the normalisation step was performed by [scran](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5112579/) R package following the [tutorial](https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/Case-study_Mouse-intestinal-epithelium_1906.ipynb) from current best practices [paper](https://www.embopress.org/doi/full/10.15252/msb.20188746).

```
Normalisation.ipynb
```

## Batch correction
There are several methods currently available for batch correction and none of them will be the best for all samples. For this reason, we used BatchBench tool to choose the best method for our normalised samples. 

[BatchBench](https://academic.oup.com/nar/article/49/7/e42/6125660?login=false) performs 7 batch_correction methods (MNN correct, BBKNN, seurat3, combat, harmony, limma and scanorama) and evaluate batch and cell type mixing through entropy and average silhouette width (asw) measures. 

Entropy is a measure of desorder, related to the probability that for each cell, its k nearest neighbors come from a different batch. Higher entropies indicated mixed batches. 

The asw measure is the average silhouette width of batches. For each point p, the average distance between p and all other points in the same batch is calculated (this is a measure of cohesion, call it A). Then, the average distance between p and all points in the nearest group is calculated (this is a measure of separation from the closest other cluster, call it B). The silhouette coefficient for p is defined as the difference between B and A divided by the greater of the two (max(A,B)) and scaled between -1-1. High batch ASW indicate well-separated batches. 

Seurat3 and Harmony showed the highest batch entropy, indicating they mixed the batches very well. The first one corrects the gene expression matrix, while the second just corrects the dimensionality reduction. Both can be used with different purposes and we choose to procees with Harmony. 
