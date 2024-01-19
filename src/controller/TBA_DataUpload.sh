#!/bin/bash

#This script goes all the way from inputing an ENCODE keywords query to tba results
#Expects the following environment variables:
### INPUTS ###
# $keyTerms = 'CT;seqTech;source;organism' -> terms to search ENCODE for data, used by searchENCODE.py and this bash script
# $numFiles = int -> number of files from query to download and analyze, used by getENCODECurl.sh (searchENCODE.py)
# $nProc = int -> number of nodes to use for run, used by runTBA.sh and bowtie2
# $out = "folder_name" -> name of folder to put TBA results, used by runTBA.sh

set -e
export STAMP=$STAMP

#Query data from ENCODE and download for tba

#module load gcc/6.2.0 

#Move the data into the tfAnalysis folder
mv '../../uploaded_data/epigenetics/'$STAMP '../../src/tfAnalysis/resources/data/fastq/'
mkdir -p "../tfAnalysis/resources/data/HOMER/"$STAMP
mkdir -p "../tfAnalysis/resources/data/IDR/"$STAMP
mkdir -p "../tfAnalysis/resources/data/bed/"$STAMP
mkdir -p "../tfAnalysis/resources/data/tba/"$STAMP
mkdir -p "../rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/CrossValidate/"$STAMP

#Define all variables for RNA-seq analysis scripts
cd ../tfAnalysis
export out=$STAMP
export bedfolder=""

#Get the organism and cell type query string
IFS=';' read -r -a array <<< "$keyTermsTBA"
organism=${array[3]}
if [ "homosapiens" = "$(echo -e "${organism}" | tr -d '[[:space:]]')" ]; then
        export species="human"
elif [ "musmusculus" = "$(echo -e "${organism}" | tr -d '[[:space:]]')" ]; then
        export species="mouse"
fi

#Runs TBA analysis for each folder in the 'fastq' folder
#Runs Fastq files ->  Bowtie2 -> (sam files) ->  HOMER -> (peaks files) ->  IDR -> (bed file) -> TBA -> (coefficients files) -> Processing -> (rank-ordered TF list)
./runAllTBA_upload.sh

#Move ENCODE search results to testing folder
cp -r "resources/data/tba/"$STAMP"/." "../rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/CrossValidate/"$STAMP
