#!/bin/bash

#This file defines some testing data to go from ENCODE keyword queries to tba results
set -e
cd ../../src/controller

### Perform RNAseq analysis for uploaded data ### 

export keyTermsRNA=';;;mus musculus'
export STAMP='MG'
export nProc=8

./RNASeq_DataUpload.sh
