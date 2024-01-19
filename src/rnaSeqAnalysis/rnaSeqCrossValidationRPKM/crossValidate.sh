#!/bin/bash

#This script is for cross-validating TBA results with RNA-seq results
#cd "resources/data/CrossValidate/"$STAMP

#module load gcc/6.2.0
#module load R/3.4.1

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
        python3 crossValidateTBARPKM.py $t"_gene_rank.tsv" $t"_rpkm.tsv" $t"_TFs.tsv" $t"_TFs_threshold" 1.000
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
