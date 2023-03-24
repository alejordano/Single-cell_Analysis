#!/bin/bash 
#SBATCH -J "download_data" 
#SBATCH --output= insert_your_path
#SBATCH -N 1 
#SBATCH -n 8 

cd /data
proc_ids=() 

while IFS= read -r line 
do 
echo $line  
nohup curl -OJ --cookie cookie $line >& /data/log/progress.out & 
proc_ids+=("$!") 
done < /data/DownloadList 

for i in ${proc_ids[@]}; do 
  wait $i 
done 
