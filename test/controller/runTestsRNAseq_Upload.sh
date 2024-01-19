#!/bin/bash

#This file defines some testing data to go from ENCODE keyword queries to tba results
set -e
cd ../../src/controller

### Perform RNAseq analysis for uploaded data ### 

export keyTermsRNA=';;;homo sapiens'
export STAMP='clonal_rd3'
export numFiles=3
export nProc=8

./RNASeq_DataUpload.sh
