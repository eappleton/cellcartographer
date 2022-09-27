#!/bin/bash
#This script runs HOMER, IDR, and TBA to go from a list of fastq files organized by cell type into folders, outputting lists of TFs most correlated to the CT

#Set environment variables
#Pick species - either mouse or human right now
set -e
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

#Load HOMER
module load gcc/6.2.0
module load homer/4.9

#Run Bowtie2 to align fastqs to a genome
for d in */ ; do
        cd $d
	#cp ../../index/* .
	nProcBowtie=$((${nProc} - 2))

	#Run Bowtie2 on all files in each directory to get a sam file for each and move sam files into HOMER directories
        module load bowtie2/2.2.9
	for f in *.fastq; do
		base="${f%.*}"
       		homerTools trim -len 30 -3 13 ${f}
		bowtie2 -x "../../index/"${genome} -U ${f}".trimmed" -p ${nProcBowtie} -S ${base}".sam" 	
	done
	module unload bowtie2/2.2.9
	#rm hg38*
	#rm mm10*
	
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

#Run HOMER to convert sam files into tsv files for IDR
cd ../../HOMER/$STAMP
for d in */ ; do
        cd $d
	
	#For HOMER, a tag directory should be made for each replicate sam file
	for f in *.sam ; do
		base="${f%.*}"
		makeTagDirectory $base ${f} -format sam
		findPeaks $base -C 0 -L 0 -fdr 0.9 -o $base"_peaks.tsv"
	done
	mkdir ./../../../IDR/$STAMP/$d
	mv *.tsv ./../../../IDR/$STAMP/$d
	cd ..
done

#Run IDR to get bed files for the cell type
cd ../../../..
module load python/3.6.0
module load idr/2.0.2
source tbaVENVpy3_4/bin/activate
#module load idr/2.0.2
#cp run_idr_homerPeaks_withBed.py resources/data/IDR_$STAMP
cd resources/data/IDR/$STAMP

for d in */ ; do
	cd $d
	
	#Checking to make sure there are tsv peak files before running this command
	count=`ls -1 *.tsv 2>/dev/null | wc -l`
	if [ $count != 0 ] ; then
		
		#Run Jenhan's IDR script
		flist=$(ls -p | grep -v / | tr '\n' ' ')
		cp ../../../../../run_idr_homerPeaks_withBed.py .
		python3 run_idr_homerPeaks_withBed.py ${flist} ./../../../bed/$STAMP -threshold 0.01
        	#python3 run_idr_homerPeaks_withBed.py ${flist} ./../../../bed/$STAMP	
	fi

	#Select bed file with the largest size (i.e. greatest number of peaks in this case)
	cd ../../../bed/$STAMP

	maxSize=0
	bigFile=""
	for f in *.bed; do
		size=$(du -sb $f | awk '{ print $1 }')
		if (( $size > $maxSize )); then
			bigFile=$f
			maxSize=$size
		fi
	done
	d_name=$(echo $d | cut -d'/' -f 1)
	cp $bigFile ${d_name}".bed"

	cd ../../IDR/$STAMP
done

#Run TBA code to get list of TFs
cd ../../../..
if [ -n "$CTs" ]; then 
	for t in "${CTs[@]}" ; do
		export CT=$t
		./runTBA.sh
	done
fi
deactivate
