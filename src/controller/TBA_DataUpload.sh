#!/bin/bash
#This script takes input BED files and makes motif score files only  
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

#Move the data into the tfAnalysis folder
mv '../../uploaded_data/epigenetics/'$STAMP '../../src/tfAnalysis/resources/data/bed/'
mkdir -p "../tfAnalysis/resources/data/tba/"$STAMP

#Define all variables for TBA scripts
cd ../tfAnalysis
export out=$STAMP

#Get the organism and cell type query string
IFS=';' read -r -a array <<< "$keyTermsTBA"
organism=${array[3]}
if [ "homosapiens" = "$(echo -e "${organism}" | tr -d '[[:space:]]')" ]; then
        export species="human"
elif [ "musmusculus" = "$(echo -e "${organism}" | tr -d '[[:space:]]')" ]; then
        export species="mouse"
fi

echo $species

#Runs TBA analysis for each folder in the 'fastq' folder
#Runs PEAR -> (merged fastq reads) ->  Bowtie2 -> (sam files) ->  HOMER -> (peaks file) ->  IDR -> (bed file) -> TBA -> (coefficients files) -> Processing -> (rank-ordered TF list)
./runTBA_bedfiles_upload.sh

#Move ENCODE search results to testing folder
#cp -r "resources/data/tba/"$STAMP"/." "../rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/CrossValidate/"$STAMP
mv "resources/data/tba/"$STAMP ../../test/tfAnalysis/results/tba
