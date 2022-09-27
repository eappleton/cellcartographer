#!/bin/bash
#SBATCH -c 8                   # Request 16 cores
#SBATCH -t 3-00:00             # Runtime in D-HH:MM format
#SBATCH -p priority             # Partition to run in
#SBATCH --mem-per-cpu=12G        # 8 GB memory needed (memory PER CORE)
#SBATCH --open-mode=append      # append adds to outfile, truncate deletes first
### In filenames, %j=jobid, %a=index in job array
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-type=END         # Mail when the job ends

#This file defines some testing data to go from ENCODE keyword queries to tba results
set -e
cd ../../src/controller

### TEST - Perform TBA and RNAseq analysis for uploaded data for refined cell lines from TFOME2 paper ### 

export keyTermsTBA=';;;homo sapiens'
export keyTermsRNA=';;;homo sapiens'
export STAMP='clonal_rd3'
export numFiles=3
export nProc=8
./RNASeq_DataUpload.sh
