#!/bin/bash

#This script runs HOMER, IDR, and TBA to go from a list of fastq files organized by cell type into folders, outputting lists of TFs most correlated to the CT

#module load gcc/6.2.0
#module load homer/4.9

set -e

#Set environment variables
#Pick species - either mouse or human right now
#species=human
if [ "human" = "$(echo -e "${species}" | tr -d '[[:space:]]')" ]; then
	export genome=hg38
elif [ "mouse" = "$(echo -e "${species}" | tr -d '[[:space:]]')" ]; then
	export genome=mm10
elif [ "elephant" = "$(echo -e "${species}" | tr -d '[[:space:]]')" ]; then
        export genome=mEleMax1
fi

#Run STAR to align fastqs to a genome
cd "resources/data/fastq/"$STAMP
for d in */ ; do
        cd $d

	#Run STAR on all files in each directory to get a sam file for each and move sam files into HOMER directories
        #module load bowtie2/2.2.9
	module load star/2.7.9a

        count0=`ls -1 *_R1.fastq 2>/dev/null | wc -l`
	if [ $count0 != 0 ]; then

                #Paired end reads alignment, followed by non-paired end reads alignment
                mkdir temp
                for f in *_R1.fastq; do
                        base="${f%_R1.*}"
                        homerTools trim -len 40 ${base}_R1.fastq
                        homerTools trim -len 20 ${base}_R2.fastq
                        #bowtie2 -x "../../index/"${genome} -1 ${base}_R1.fastq".trimmed" -2 ${base}_R2.fastq".trimmed" -p ${nProc} -S ${base}".sam"
                        STAR --runThreadN ${nProc} --genomeDir "../../index/"${genome} --readFilesIn ${base}_R1.fastq".trimmed" --outFileNamePrefix ${base}"_R1"
			STAR --runThreadN ${nProc} --genomeDir "../../index/"${genome} --readFilesIn ${base}_R2.fastq".trimmed" --outFileNamePrefix ${base}"_R2"

			###This time only, only consider the first read###
			#bowtie2 -x "../../index/"${genome} -U ${base}_R1.fastq".trimmed" -p ${nProc} -S ${base}".sam"

			mv ${base}_R1.fastq* temp
                        mv ${base}_R2.fastq* temp
                done

                cd temp; mv * ..; cd ..; rm -r temp;
        else
                #Non-paired end reads alignment
                for f in *.fastq; do
                        base="${f%.*}"
                        homerTools trim -len 40 ${f}
                        #bowtie2 -x "../../index/"${genome} -U ${f}".trimmed" -p ${nProc} -S ${base}".sam"
			STAR --runThreadN ${nProc} --genomeDir "../../index/"${genome} --readFilesIn ${f}".trimmed" --outFileNamePrefix ${f}
                done
        fi	
        
	#module unload bowtie2/2.2.9
	
	#Only make new folders for cell types with sam files and add to CTs variable
	count=`ls -1 *.sam 2>/dev/null | wc -l`
	if [ $count != 0 ]; then
		#d_name=$(echo $d | cut -d'/' -f 1)
		CTs+=($(echo ${d%%/}))
		mkdir ./../../../HOMER/$STAMP/$d
		mv *.sam ./../../../HOMER/$STAMP/$d
	fi
        cd ..
done

#Run HOMER to produce exon counts for cross-validation with TBA results
cd ../../HOMER/$STAMP
for d in */ ; do
        cd $d
	
	#For HOMER, a tag directory should be made for each replicate sam file
	dirbases=()
	for f in *.sam ; do
		base="${f%.*}"
		makeTagDirectory $base ${f} -format sam
		dirbases+=$(echo ${base}" ")
	done
	
	d_name=$(echo ${d%%/})

	#analyzeRepeats.pl rna ${genome} -count exons -condenseGenes -rpkm -d $dirbases > $d_name"_rpkm.tsv"
	if [ "mEleMax1" = "$(echo -e "${genome}" | tr -d '[[:space:]]')" ]; then
		analyzeRepeats.pl ../../genes.gtf none -count exons -condenseGenes -noadj -d $dirbases > $d_name"_RAWREADS.tsv"
		analyzeRepeats.pl ../../genes.gtf none -count exons -condenseGenes -rpkm -d $dirbases > $d_name"_RPKM.tsv"
		#analyzeRepeats.pl ../../GCF_024166365.1_mEleMax1_primary_haplotype_genomic_P53RTG.gtf none > $d_name"_rpkm.tsv"
	else
		analyzeRepeats.pl rna ${genome} -count genes -condenseGenes -noadj -strand both -d $dirbases > $d_name"_RAWREADS.tsv"
		analyzeRepeats.pl rna ${genome} -count genes -condenseGenes -rpkm -strand both -d $dirbases > $d_name"_RPKM.tsv"
	fi

	mv *.tsv ./../../../RPKM/$STAMP/
	cd ..
done
