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

##Move the data into the tfAnalysis folder
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

echo $species

#Runs TBA analysis for each folder in the 'fastq' folder
#Runs PEAR -> (merged fastq reads) ->  Bowtie2 -> (sam files) ->  HOMER -> (peaks file) ->  IDR -> (bed file) -> TBA -> (coefficients files) -> Processing -> (rank-ordered TF list)
./runAllTBA_upload.sh

#Move ENCODE search results to testing folder
cp -r "resources/data/tba/"$STAMP"/." "../rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/CrossValidate/"$STAMP

#Move results back into the testing and scratch folders#
#mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/tfAnalysis/results/fastq
mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/tfAnalysis/results/HOMER
mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/tfAnalysis/results/IDR
mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/tfAnalysis/results/bed

mv "resources/data/fastq/"$STAMP '../../uploaded_data/epigenetics/' 
mv "resources/data/HOMER/"$STAMP /n/scratch3/users/e/ema26/combinatorial_screen_software/test/tfAnalysis/results/HOMER
mv "resources/data/IDR/"$STAMP /n/scratch3/users/e/ema26/combinatorial_screen_software/test/tfAnalysis/results/IDR
mv "resources/data/bed/"$STAMP /n/scratch3/users/e/ema26/combinatorial_screen_software/test/tfAnalysis/results/bed
mv "resources/data/tba/"$STAMP ../../test/tfAnalysis/results/tba

#Move the data into the rnaSeqCrossValidationRPKM folder
mv '../../uploaded_data/transcriptomics/'$STAMP '../rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/fastq/'
mkdir -p "../rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/HOMER/"$STAMP
mkdir -p "../rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/RPKM/"$STAMP
#mkdir -p "../rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/IDR/"$STAMP

#Define all variables for RNA-seq analysis scripts
cd ../rnaSeqAnalysis/rnaSeqCrossValidationRPKM
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

#Run cross-validation script
cp -r "resources/data/RPKM/"$STAMP"/." "resources/data/CrossValidate/"$STAMP
cp crossValidate.sh "resources/data/CrossValidate/"$STAMP
cp crossValidateTBARPKM.py "resources/data/CrossValidate/"$STAMP
cp densityPlots.R "resources/data/CrossValidate/"$STAMP
cd "resources/data/CrossValidate/"$STAMP
./crossValidate.sh
cd ../../../..

#Move results back into the testing and scratch folders
#mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/fastq
mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/HOMER
mkdir -p /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/IDR

#mv "resources/data/fastq/"$STAMP /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/fastq
mv "resources/data/fastq/"$STAMP '../../../uploaded_data/transcriptomics/' 
mv "resources/data/HOMER/"$STAMP /n/scratch3/users/e/ema26/combinatorial_screen_software/test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/HOMER
mv "resources/data/RPKM/"$STAMP ../../../test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/RPKM
mv "resources/data/CrossValidate/"$STAMP ../../../test/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/results/CrossValidate
