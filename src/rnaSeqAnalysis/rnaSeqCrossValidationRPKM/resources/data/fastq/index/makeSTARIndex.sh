#!/bin/bash

#This script creates STAR indexes for alignments
module load gcc/6.2.0
module load star/2.7.9a

#export nProc=8

STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir hg38 \
--genomeFastaFiles /n/groups/church/eappleton/combinatorial_screen_software/src/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/fastq/index/hg38/hg38.fa \
--sjdbGTFfile /n/groups/church/eappleton/combinatorial_screen_software/src/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/fastq/index/hg38/hg38.ncbiRefSeq.gtf \
--sjdbOverhang 99

STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir mm10 \
--genomeFastaFiles /n/groups/church/eappleton/combinatorial_screen_software/src/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/fastq/index/mm10/mm10.fa \
--sjdbGTFfile /n/groups/church/eappleton/combinatorial_screen_software/src/rnaSeqAnalysis/rnaSeqCrossValidationRPKM/resources/data/fastq/index/mm10/mm10.ncbiRefSeq.gtf \
--sjdbOverhang 99
