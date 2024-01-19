#!/bin/bash

#This file defines some testing data to go from ENCODE keyword queries to tba results
set -e
cd ../../src/controller

### TEST - Perform TBA and RNAseq analysis for uploaded data for refined cell lines from TFOME2 paper ### 

export keyTermsTBA=';;;mus musculus'
export keyTermsRNA=';;;mus musculus'

export numFiles=3
export nProc=16
export STAMP='MG'

./TBARNASeqCrossValidation_DataUpload.sh
