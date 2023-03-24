# Network analysis

### Install the required R packages using:
```
#!/usr/bin/bash
while IFS=" " read -r package version; 
do 
Rscript -e "devtools::install_version('"$package"', version='"$version"')"; 
done < "RequirementNetwork.txt"
Rscript -e devtools::install_github("iaconogi/bigSCale2")
```

### Build co-expression networks based bigSCale2 tool [GitHub BigSCale2](https://github.com/iaconogi/bigSCale2)
[Reference Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1713-4)\
Build co-expression networks for each subpopulation at both diagnosis and relapse timepoints separately.\
All the networks built are available on the workstation "/data/aconforte/edgeAnalysis/Results/".

```
Rscript BuildingNetworks.R
```

### Quality control of the networks
The degree distribution in biological networks are expected to follow a power law, where few nodes have high degree and most nodes have low degree.\
To check the degree distribution of the networks built run:
```
Rscript QualityControl.R
``` 

As observed in Supplementary Figure 2, the power law is only observed in networks built for subpopulations with more than 100 cells in a timepoint. 

### Evaluate the rewiring of signal transduction from diagnosis to relapse

STEP 1 - Identify nodes (genes) with at least 25% increase in centrality measure (degree and betweenness):
```
Rscript GetNodeMeasures.R
Rscript GetNodesWithIncreasedCent.R
```

STEP 2 - Identify the difference in GO terms enriched by genes interacting with central nodes identified in the previous step:
```
Rscript GetHubsBnGOterms.R
Rscript DiffBP.R
```