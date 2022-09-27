#!/bin/bash
#This is a script for unpacking .bz2 or .gz fastq files and assembling reads with PEAR.

#First unpack .bz2 files into .fastq files for PEAR
count=`ls -1 *.bz2 2>/dev/null | wc -l`
if [ $count != 0 ]; then
	for i in *.bz2; do
		bunzip2 $i
	done
fi

#Unpack .gz files into .fastq files for PEAR
count_gz=`ls -1 *.gz 2>/dev/null | wc -l`
if [ $count_gz != 0 ]; then
        for i in *.gz; do
                gunzip $i
        done
fi

echo "RAN PEAR"

#Next we will run PEAR for each pair of F and R sequencing results
#Find the unique file prefixes first
#ALL_FASTQ=`ls | grep '.*R.\.fastq'`

function prefix 
{
	for file in `ls | grep '.*R.\.fastq'`; 
	do
		echo "${file%_*.fastq*}";
	done | sort -u
}

UNIQUE_PREFIXES=$(prefix)
echo $UNIQUE_PREFIXES

#For each unique file prefix, run PEAR 
#module load seq/pear/0.9.6

for p in $UNIQUE_PREFIXES
do
	#echo $p"_R1.fastq"
	if [ ! -f $p".fastq" ]; then		
		pear -f $p"_R1.fastq" -r $p"_R2.fastq" -o $p
		echo $p".fastq DID NOT EXIST"
		
		#Right now we aren't going to do anything with the unassembled, discarded or indexCount files
        	#We will clean up the files we do not care about here for now
        	rm $p"_R1.fastq.indexCount.tab"
        	rm $p".discarded.fastq"
        	rm $p".unassembled.forward.fastq"
        	rm $p".unassembled.reverse.fastq"
		mv $p".assembled.fastq" $p".fastq"
	else
		echo $p".fastq EXISTS"
	fi
done
