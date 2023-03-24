
#!/bin/bash 
#SBATCH -J "cellranger" 
#SBATCH --output=/data/log/cellrangerAll.out 
#SBATCH -p normal 
#SBATCH -N 1 
#SBATCH -n 8 

export PATH=path_to_cellranger-6.0.0:$PATH 

transcriptome="/reference/refdata-gex-GRCh38-2020-A" 
s=(Patient1 Patient2 Patient3) 

for j in ${s[@]} 
do 
data_dir="/data/$j" 
sams="$(grep $j /data/$j/samples_names)" 

    for i in ${sams[@]} 
        do 
        cd $data_dir 
        cellranger count --id=$i --transcriptome=$transcriptome --fas 
        tqs=$data_dir --sample=$i 
    done 
done 