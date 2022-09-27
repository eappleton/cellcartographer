#!/bin/bash
#SBATCH -c 4                   # Request 16 cores
#SBATCH -t 24:00:00             # Runtime in D-HH:MM format
#SBATCH -p priority             # Partition to run in
#SBATCH --mem-per-cpu=8G        # 16 GB memory needed (memory PER CORE)
#SBATCH --open-mode=append      # append adds to outfile, truncate deletes first
### In filenames, %j=jobid, %a=index in job array
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-type=END         # Mail when the job ends

#This script is for cross-validating TBA results with RNA-seq results
#cd "resources/data/CrossValidate/"$STAMP

#Load modules
module load gcc/6.2.0
module load R/3.4.1

#Find all cell types
CTs=()
arr=( $(ls *_gene_rank.tsv) )
for t in "${arr[@]}" ; do
	CT="${t%_gene*}"
	CTs+=($(echo ${CT}))	
done

echo "${CTs[@]}"

#For each cell type, run cross-validation python script and plot density functions for RPKM
for t in "${CTs[@]}" ; do
	export CT=$t
        python crossValidateTBARPKM.py $t"_gene_rank.tsv" $t"_rpkm.tsv" $t"_TFs.tsv" $t"_TFs_threshold" 1.000
	Rscript densityPlots.R
done

#Move files into other folders for order
mkdir TBA_results
mkdir RPKM_results
mkdir PNG
mkdir TF
mv *_gene_rank.tsv TBA_results
mv *TFs* TF
mv *_rpkm.tsv RPKM_results
mv *.png PNG
