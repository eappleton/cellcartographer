### This python script is used to take a list of barcodes from sequencing files and compare distance to list of known barcodes to tally apearance frequency using both Hamming and Levenstein distance

__author__ = "Evan Appleton"
__email__ = "evan.appleton@gmail.com"

### imports ###
import sys
import os
import argparse
import pandas as pd
import distance as dis

### functions ###

### main operations ###
if __name__ == "__main__":
	execPath = os.path.dirname(os.path.realpath(__file__))

	# build argument parser
	parser = argparse.ArgumentParser(description='Given some TBA and RPKM results, markup the TBA file with RPKM')
	parser.add_argument("-BC_file",
		help="CSV file with raw barcodes extracted from fastq file",
		default="", type=str)
	parser.add_argument("-BC_lib",
		help="Barcode library mapped to gene names",
		default="../PBAN_barcodes.csv", type=str)
	parser.add_argument("-CTmap",
		help="Map of which genes are expected in which cell type screen",
		default="../CT_BC_map.csv", type=str)
	parser.add_argument("-fastqCTmap",
                help="Map of fastq file names to cell type",
                default="../fastq_CT_map.csv", type=str)
	parser.add_argument("-outPath",
		help="Name of output file",
		default="", type=str)

	# parse arguments
	args = parser.parse_args()
	BCs = args.BC_file
	lib = args.BC_lib
	CT_BC = args.CTmap
	fastq_CT = args.fastqCTmap
	out = args.outPath

	# Read in BC library to make dictionary
	BCdict = {}
	BCct = {}
	TFct = {}
	BCdict['unknown'] = 'unknown'
	BCct['unknown'] = 0
	TFct['unknown'] = 0
	f = open(lib, 'rU')
	BClib = f.readlines()[1:]

	for line in BClib:
		TF=line.replace('\n',"").split(',')[0]
		BC=line.replace('\n',"").split(',')[1]
		BCdict[BC] = TF
		BCct[BC] = 0
		TFct[TF] = 0

	#print(BCdict)

	# Upload raw BC data, loop through to get Hamming or Levenstein distance from each BC in the BC library
	BCdf = pd.read_csv(BCs, header=0, skip_blank_lines=True)
	BClist = list(BCdf.ix[:,0])
	unknown_BC_file = os.path.splitext(BCs)[0] + '_UNKNOWN_BC.txt'

	for BC in BClist:
		minDistance = 3
		match = 'unknown'
		for BCl in BCct.keys():					
			d = 20

			#Get distance
			if len(BC) == len(BCl):
				d = dis.hamming(BC, BCl)
			else:
				d = int(dis.levenshtein(BC, BCl))

			if d < minDistance:
				minDistance = d
				match = BCl
		
		#Save unmatched barcodes to file for further analysis
		if match == 'unknown':
			with open (unknown_BC_file, 'a') as f:
				f.write(BC + '\n')

		ct = BCct[match]
		ct = ct + 1
		TF = BCdict[match]
		TFct[TF] = ct
		BCct[match] = ct

	# Determine which CT this file was so that we can subset the TFs plotted
	fqCTdf = pd.read_csv(fastq_CT)
	fqCTdict = fqCTdf.set_index('FASTQ').T.to_dict('list')	
	fqFname = os.path.splitext(BCs)[0] + '.fastq'
	CT = fqCTdict[fqFname][0]
	ctTFdf = pd.read_csv(CT_BC)

	# Save the counts to a csv file
	df = pd.DataFrame.from_dict(TFct, orient="index")
	df.columns = ['BC_COUNTS']
	all_out = os.path.splitext(out)[0] + '_ALL.csv'
	
	df.to_csv(all_out, sep=',')
	df = pd.read_csv(all_out) #This saves for all barcodes in case samples were swapped somehow
	subDF = df[df['Unnamed: 0'].isin(ctTFdf[CT])]
	subDF.to_csv(out, sep=',', index=False)
	df.to_csv(all_out, sep=',', index=False)
