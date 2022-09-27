#This script takes a file name with a list of ENCODE fastq data files and downloads the first $numFiles files
#Make a new folder to download the data and a temporary file to read to subset data and download a few

file="search_results.txt"
ctFile="search_results_ct.txt"
preFile="search_results_prefixes.txt"
extFile="search_results_extensions.txt"
sizeFile="search_results_size.txt"

#Loop through the ctFile to get each unique type and only grab ${numFiles} of that type, create temp files
#Doing a nested loop here might not be best, but it should be ok as long as the query results aren't too long
### CTs=( $(uniq $ctFile) )
declare -a CTs
CTs+=$(head -n 1 $ctFile)

readarray ctArrayAll < $ctFile
readarray sizeArrayAll < $sizeFile
readarray fileArrayAll < $file
readarray prefixArrayAll < $preFile
readarray fileExtArrayAll < $extFile

#Loop through CTs
for ct in "${CTs[@]}";
do
	count=0

	#Loop through the whole ctFile, append the first n lines to each of the growing files
	for (( i = 0 ; i < ${#ctArrayAll[@]} ; i=$i+1 ));
	do
	
		#If the cell type in the ctFile matches the cell type from CTs
		if [ $ct == ${ctArrayAll[${i}]} ]; then

			#Check that the size is > 100MB and less than ~10GB for TBA, but no limit for RNA-seq
			if [ -z ${RNA+x} ]; then
				if (( ${sizeArrayAll[${i}]} > 270000000 && ${sizeArrayAll[${i}]} < 2700000000 )); then

					#Only add $numFiles of replicates to the temporary files
					if (( $count < $numFiles )); then
						((count++))
						echo ${fileArrayAll[${i}]} >> tmp.txt
						echo ${ctArrayAll[${i}]} >> tmpCt.txt
						echo ${prefixArrayAll[${i}]} >> tmpPre.txt
						echo ${fileExtArrayAll[${i}]} >> tmpExt.txt
					fi
				fi
			else
				#Only add $numFiles of replicates to the temporary files
                                if (( $count < $numFiles )); then
                                        ((count++))
                                        echo ${fileArrayAll[${i}]} >> tmp.txt
                                        echo ${ctArrayAll[${i}]} >> tmpCt.txt
                                        echo ${prefixArrayAll[${i}]} >> tmpPre.txt
                                        echo ${fileExtArrayAll[${i}]} >> tmpExt.txt
                                fi
			fi
		fi		
	done

	#Files without replicates will crash in the IDR step, so we will remove types with only one fastq file
	if (( $count < 2 )); then
		tail -n 1 tmp.txt | wc -c | xargs -I {} truncate tmp.txt -s -{}
		tail -n 1 tmpCt.txt | wc -c | xargs -I {} truncate tmpCt.txt -s -{}
		tail -n 1 tmpPre.txt | wc -c | xargs -I {} truncate tmpPre.txt -s -{}
		tail -n 1 tmpExt.txt | wc -c | xargs -I {} truncate tmpExt.txt -s -{}
	fi
done

#Grab each file with curl
readarray ctArray < tmpCt.txt
readarray fileArray < tmp.txt
readarray prefixArray < tmpPre.txt
readarray fileExtArray < tmpExt.txt	 
parallel --jobs ${nProc} curl -s -O -L https://www.encodeproject.org{} ::: ${fileArray[*]}   

#Rename files with descriptive extensions
for (( i = 0 ; i < ${#prefixArray[@]} ; i=$i+1 ));
do
	prefix=${prefixArray[${i}]}
	extension=${fileExtArray[${i}]}
	ct=${ctArray[${i}]}
	dest=$(echo $ct"/"$prefix$extension | tr -d ' ')
	st=$(echo $prefix".fastq.gz" | tr -d ' ')
	if [ -f ./$st ]; then
    		mv ./$st ./$dest	
	fi
done

#Cleanup
rm tmp.txt
rm tmpCt.txt
rm tmpPre.txt
rm tmpExt.txt

#Removes empty data directories from downstream analysis crashes
find . -type d -empty -delete
