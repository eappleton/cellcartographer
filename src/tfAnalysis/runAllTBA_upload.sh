#!/bin/bash

#This script runs HOMER, IDR, and TBA to go from a list of fastq files organized by cell type into folders, outputting lists of TFs most correlated to the CT

#Set environment variables
#Pick species - either mouse or human right now
set -e

if [ "human" = "$(echo -e "${species}" | tr -d '[[:space:]]')" ]; then
	export genome=hg38
elif [ "mouse" = "$(echo -e "${species}" | tr -d '[[:space:]]')" ]; then
	export genome=mm10
fi

cd "resources/data/fastq/"$STAMP

#module load gcc/6.2.0
#module load homer/4.9

#Run Bowtie2 to align fastqs to a genome
for d in */ ; do
        cd $d

	#Run Bowtie2 on all files in each directory to get a sam file for each and move sam files into HOMER directories
        module load bowtie2/2.2.9

        count0=`ls -1 *_R1.fastq 2>/dev/null | wc -l`
        if [ $count0 != 0 ]; then
		
		#Paired end reads alignment, followed by non-paired end reads alignment
		mkdir temp
                for f in *_R1.fastq; do

                        base="${f%_R1.*}"
                        homerTools trim -len 40 ${base}_R1.fastq
                        homerTools trim -len 40 ${base}_R2.fastq
                        bowtie2 -x "../../index/"${genome} -1 ${base}_R1.fastq".trimmed" -2 ${base}_R2.fastq".trimmed" -p ${nProc} -S ${base}".sam"

			mv ${base}_R1.fastq* temp
			mv ${base}_R2.fastq* temp

                done

		cd temp; mv * ..; cd ..; rm -r temp;
        else
		
		#Non-paired end reads alignment
		for f in *.fastq; do
                	base="${f%.*}"
                	homerTools trim -len 10000 ${f}
                	bowtie2 -x "../../index/"${genome} -U ${f}".trimmed" -p ${nProc} -S ${base}".sam"
        	done
	fi

	mkdir ./../../../HOMER/$STAMP/$d
	mv *.sam ./../../../HOMER/$STAMP/$d

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

#module load python/3.6.0
#module load idr/2.0.2

source tbaVENV/bin/activate
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
	
		CTs+=($(echo ${f%%.*}))
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
		#export out=$t
		./runTBA.sh
	done
fi
deactivate
