#This script subsets the GEO search results and uses SRATools to grab fastq files
#Make abbreviated file from full search
cd $folder
head -n 10 $file > tmp.csv

#Loop through csv file to make 2 arrays - one of paired end reads, one for unpaired 
pairedFileArray=()
unpairedFileArray=()
awk -F "\"*,\"*" -v OFS=, '{print $1,$10,$16}' tmp.csv > tmpS.csv

{ read
while IFS=, read -r c1 c2 c3
do
	if [ "${c3}" == "SINGLE" ]; then
		unpairedFileArray+=("$c1")
	else
		pairedFileArray+=("$c1")
	fi
done
} < tmpS.csv

#Grab each file with SRA Toolkit
parallel --jobs 8 fastq-dump {} ::: ${unpairedFileArray[*]}
parallel --jobs 8 fastq-dump -I --split-files {} ::: ${pairedFileArray[*]}

#Cleanup
rm tmp.csv
rm tmpS.csv
cd ..
