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
cd "resources/data/fastq/"$STAMP

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
                        ###This time only, only consider the first read###
			#bowtie2 -x "../../index/"${genome} -U ${base}_R1.fastq".trimmed" -p ${nProc} -S ${base}".sam"

			mv ${base}_R1.fastq* temp
                        mv ${base}_R2.fastq* temp
                done

                #Non-paired end reads alignment if they exist in the same folder as paired-end fastq files
                #for f in *.fastq; do
                #        base="${f%.*}"
                #        homerTools trim -len 30 ${f}
                #        bowtie2 -x "../../index/"${genome} -U ${f}".trimmed" -p ${nProc} -S ${base}".sam"
                #done
                cd temp; mv * ..; cd ..; rm -r temp;
        else
                #Non-paired end reads alignment
                for f in *.fastq; do
                        base="${f%.*}"
                        homerTools trim -len 40 ${f}
                        bowtie2 -x "../../index/"${genome} -U ${f}".trimmed" -p ${nProc} -S ${base}".sam"
                done
        fi	
        
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
	analyzeRepeats.pl rna ${genome} -count exons -condenseGenes -noadj -d $dirbases > $d_name"_rpkm.tsv"

	mv *.tsv ./../../../RPKM/$STAMP/
	cd ..
done
