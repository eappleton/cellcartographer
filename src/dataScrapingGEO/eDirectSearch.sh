#This is the the script to run the Entrez queries of GEO (SRA) via E-Direct tools
#The returned CSV file needs to be cleaned up
file=$dir"_results.csv"
mkdir $dir
awk 'NF' search.csv > $dir"/"$file
cd $dir
header=$(head -n 1 $file)
sed '/^Run/ d' $file > $file".tmp"
sed -i 1i$header $file".tmp"
mv $file".tmp" $file
cd ..
rm search.csv
