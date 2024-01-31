#This simple script downloads data from a list of SRA or ENCODE entries

#Grab ENCODE data
readarray fileArray < encode_files_EXAMPLE.txt
parallel --jobs 8 curl -s -O -L https://www.encodeproject.org{} ::: ${fileArray[*]} 

#Grab SRA data
#readarray pairedFileArray < SRA_files_EXAMPLE.txt
#parallel --jobs 8 fasterq-dump {} ::: ${pairedFileArray[*]}
