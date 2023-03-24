# Obtaining count matrix for each 10x genomics sample with cell ranger.

## Create folders reference, data and script. 
```
mkdir reference
mkdir data
mkdir script
```

## Download the human reference, check for updated versions.  
```
cd reference
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz 
tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz 
```

## Download your data available on links listed in DownloadList:
```
./downloadData.sh
```

## Check names of the files:
The name of fastq files need to be "Name_S001_L1_R2_001.fastq.gz". If changes are needed, the following commands can be used:
```
rename s/old/new/ files 
rename s/Patient_Timepoint*/Patient_Timepoint/ *fastq.gz
rename PatientTimepoint Patient_Timepoint ./Patient/Patient_Timepoint* #command for Slurm
for f in *.fastq.gz; do mv -- "$f" "${f%001_*.fastq.gz}001.fastq.gz"; done #take the "part" out in *_001_part.fastq.gz
```

## Create a file for each samples and run cell ranger to obtain the count matrix:
```
cd script
ls *fastq.gz | cut -d "_" -f 1-2 | sort -u > ../script/samples_names.txt 
./cellranger_count.sh
```

## Quality Control:
The quality control of the count matrix considered the number of counts, number of genes and percentage of mitochondrial genes.\
To run the jupyter notebook in python, you need to install the libraries listes in requirements.txt.
```
pip install -r requirements.txt
```
The correspondent codes are also available for R.