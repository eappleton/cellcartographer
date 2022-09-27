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

### TEST - Perform TBA and RNAseq analysis for uploaded data for refined cell lines from TFOME2 paper ### 
export STAMP='clonal_rd2'

set -e

#Move results back into the testing and scratch folders
#mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/fastq
mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/HOMER
mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/IDR

#mv "resources/data/fastq/"$STAMP /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/fastq
mv "resources/data/HOMER/"$STAMP '../../../uploaded_data/transcriptomics/'
#mv "resources/data/IDR/"$STAMP /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/IDR
mv "resources/data/RPKM/"$STAMP ../../../test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/RPKM
mv "resources/data/CrossValidate/"$STAMP ../../../test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/CrossValidate

