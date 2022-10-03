### This script does analysis for barcode appearces for all fastq files in a data directory ###
#mkdir fastq
#mkdir barcode

# Run Alex Ng's Perl script to extract barcodes from raw fastq reads
# USAGE EX: 'perl 20190104-extract-TFbarcode-from-RNAseq.pl 2019-01-02-Evan-test.fastq' #
#cd fastq
#for file in *.fastq; do
#	perl ../extractBC.pl $file
#	newFileName="${file%.*}"
#	mv ${file}.csv ${newFileName}.csv
#done
#mv *.csv ../barcode
#cd ../barcode

cd barcode
# Use a python script to count barcodes using Hamming/Levinstein distance and plot results in R
for file in *.csv; do
	newFileName="${file%.*}"
	python ../tallyBC.py -BC_file $file -outPath ${newFileName}_BC_counts.csv
done
mkdir results
mv *BC_counts* results
mv *.txt results
cd ..

# Use an R script to plot the data
Rscript analyze_BC_results.R
