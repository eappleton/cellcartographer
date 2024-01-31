#!/bin/bash

#This file defines some testing data to go from ENCODE keyword queries to tba results
set -e
cd ../../src/controller

### Perform TBA analysis for all uploaded fastq files ### 

export numFiles=3
export nProc=8

export keyTermsTBA=';;;mus musculus'

export STAMP='MG'

./TBA_DataUpload.sh
