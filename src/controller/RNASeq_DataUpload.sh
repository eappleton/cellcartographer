#!/bin/bash
#This script goes all the way from inputing an ENCODE keywords query to tba results
#Expects the following environment variables:
### INPUTS ###
# $keyTerms = 'CT;seqTech;source;organism' -> terms to search ENCODE for data, used by searchENCODE.py and this bash script
# $numFiles = int -> number of files from query to download and analyze, used by getENCODECurl.sh (searchENCODE.py)
# $nProc = int -> number of nodes to use for run, used by runTBA.sh and bowtie2
# $out = "folder_name" -> name of folder to put TBA results, used by runTBA.sh

set -e
#timestamp() {  STAMP=`date +%Y_%m_%d_%H_%M_%S`; }
#timestamp
export STAMP=$STAMP

#Query data from ENCODE and download for tba
module load gcc/6.2.0 

#Move the data into the rnaSeqCrossValidationRPKM folder
mv '../../uploaded_data/transcriptomics/'$STAMP '../rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/fastq/'
mkdir -p "../rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/HOMER/"$STAMP
mkdir -p "../rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/RPKM/"$STAMP

#Define all variables for RNA-seq analysis scripts
cd ../rnaSeqAnalysis/rnaSeqCrossValidationRPKM
export out=$STAMP
export bedfolder=""

IFS=';' read -r -a array <<< "$keyTermsRNA"
organism=${array[3]}
if [ "homosapiens" = "$(echo -e "${organism}" | tr -d '[[:space:]]')" ]; then
        export species="human"
elif [ "musmusculus" = "$(echo -e "${organism}" | tr -d '[[:space:]]')" ]; then
        export species="mouse"
fi

#Runs RNA-seq analysis for each folder in the 'fastq' folder
#Runs PEAR -> (merged fastq reads) ->  Bowtie2 -> (sam files) ->  HOMER -> (peaks file) -> RPKM files 
./runRNASeqAnalysisENCODE_upload.sh

#Move results back into the testing and scratch folders
#mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/fastq
mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/HOMER
mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/IDR

mv "resources/data/fastq/"$STAMP ../../../uploaded_data/transcriptomics 
mv "resources/data/HOMER/"$STAMP /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/HOMER
mv "resources/data/RPKM/"$STAMP ../../../test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/RPKM
