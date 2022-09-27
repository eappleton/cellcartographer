#!/bin/bash
#SBATCH -c 8                   # Request 16 cores
#SBATCH -t 0-12:00             # Runtime in D-HH:MM format
#SBATCH -p priority             # Partition to run in
#SBATCH --mem-per-cpu=16G        # 8 GB memory needed (memory PER CORE)
#SBATCH --open-mode=append      # append adds to outfile, truncate deletes first
### In filenames, %j=jobid, %a=index in job array
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-type=END         # Mail when the job ends

#This file defines some testing data to go from ENCODE keyword queries to tba results
set -e
cd ../../src/controller

### TEST - Perform TBA analysis for all uploaded ENCODE bed files for TFOME2 paper ### 

export numFiles=3
export nProc=8

export keyTermsTBA=';;;mus musculus'

export STAMP='astrocyte'
./TBA_DataUpload2.sh
