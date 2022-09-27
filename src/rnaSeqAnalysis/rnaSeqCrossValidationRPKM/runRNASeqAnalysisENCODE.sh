#!/bin/bash
#This script runs HOMER, IDR, and TBA to go from a list of fastq files organized by cell type into folders, outputting lists of TFs most correlated to the CT
module load gcc/6.2.0
module load homer/4.9
set -e

#Set environment variables
#Pick species - either mouse or human right now
#species=human
if [ "human" = "$(echo -e "${species}" | tr -d '[[:space:]]')" ]; then
	export genome=hg38
elif [ "mouse" = "$(echo -e "${species}" | tr -d '[[:space:]]')" ]; then
	export genome=mm10
fi

#Run PEAR to merge reads if paired
cp runPEAR.sh "resources/data/fastq/"$STAMP
cd "resources/data/fastq/"$STAMP

#Find all cell types based upon folder names in the fastq folder
CTs=()

#Run PEAR to merge all paired end read files
for d in */ ; do
	cd $d
	cp ../runPEAR.sh ./
	./runPEAR.sh
	rm runPEAR.sh
	cd ..
done
rm runPEAR.sh

pwd
#Run Bowtie2 to align fastqs to a genome
for d in */ ; do
        cd $d

	#Run Bowtie2 on all files in each directory to get a sam file for each and move sam files into HOMER directories
        module load bowtie2/2.2.9
	for f in *.fastq; do
		base="${f%.*}"
       		homerTools trim -len 30 -3 13 ${f}
		#homerTools trim -len 30 ${f}
		bowtie2 -x "../../index/"${genome} -U ${f}".trimmed" -p ${nProc} -S ${base}".sam" 	
	done
	module unload bowtie2/2.2.9
	
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
#module load homer/4.9
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
	analyzeRepeats.pl rna ${genome} -count exons -condenseGenes -rpkm -d $dirbases > $d_name"_rpkm.tsv"

	mv *.tsv ./../../../RPKM/$STAMP/
	cd ..
done
