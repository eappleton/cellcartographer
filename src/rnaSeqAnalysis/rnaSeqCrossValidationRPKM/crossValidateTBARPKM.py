#! /usr/bin/env python
################################################################################
'''
Given some TBA and RPKM results, markup the TBA file with RPKM
'''

### header ###
__author__ = "Evan Appleton"
__email__ = "evan_appleton@hms.harvard.edu"

### imports ###
import sys
import os
import argparse

### functions ###


### main operations ###
if __name__ == "__main__":
	execPath = os.path.dirname(os.path.realpath(__file__))
    
	# build argument parser
	parser = argparse.ArgumentParser(description='Given some TBA and RPKM results, markup the TBA file with RPKM')
	parser.add_argument("TBAFile",
		help="TBA file name",
		default="", type=str)
	parser.add_argument("RPKMFile",
        	help="RPKM file name",
        	default="", type=str)
	parser.add_argument("outfile",
		help="Marked up TBA gene rank file name",
		default="gene_rank_TBA_RPKM.tsv", type=str)
	parser.add_argument("outfile_threshold",
                help="Marked up TBA gene rank file with threshold filter name",
                default="gene_rank_TBA_RPKM_threshold", type=str)
	parser.add_argument("threshold",
                help="RPKM threshold",
                default="0.000", type=float) 

	#parse arguments
	args = parser.parse_args()
	TBAFile = args.TBAFile
	RPKMFile = args.RPKMFile
	outfile = args.outfile 
	threshold = args.threshold
	outfile_thr = args.outfile_threshold
	outfile_thr = (outfile_thr + '_' + str(threshold) + '_rpkm.tsv')

	#Read in RPKM file, line by line and make dictionary of gene name to RPKM
	geneRPKM = {}
	with open(RPKMFile) as f:
		next(f)
		for line in f:
	
			#Parse the line into array split by \t
			values = line.split("\t")
			geneList = values[7]
			RPKM = values[8].replace('\n',"")

			#Split the 7th column into an array, add all genes to dictionary with value in the 8th column if it's a protein coding gene
			l = geneList.split("|")
			geneType = l[len(l)-1]

			if geneType == 'protein-coding':
				genes = l[0:len(l)-3]
				for gene in genes:
					GENE = gene.upper()
					geneRPKM[GENE] = RPKM


	#Read in TBA file, line by line and look up the genes in the RPKM dictionary, and build to the outfile
	out = open(outfile, 'w')
	out_thr = open(outfile_thr, 'w')

	#Write header to outfile
	with open(TBAFile) as f:
		header = f.readline()
		out.write(header.replace('\n','\t')+'RPKM\n')
		out_thr.write(header.replace('\n','\t')+'RPKM\n')

	#Loop through TBA file
	with open(TBAFile) as f2:
		next(f2)
		for line in f2:
			
			#Parse the line into array split by \t
			values = line.split("\t")
			values = [v.replace('\n','') for v in values]
			genes = values[0].split("|")

			#If this gene is in the dictionary, add it to the outfile
			for gene in genes:
				GENE = gene.upper()	
				if GENE in geneRPKM.keys():
					out.write(gene+'\t'+values[1]+'\t'+values[2]+'\t'+geneRPKM.get(GENE)+'\n')
					#out.write(line.replace('\n','\t')+geneRPKM.get(gene)+'\n')
					#break
					
					#If greater than the threshold, add to filtered outfile
					if float(geneRPKM.get(GENE)) > threshold:
						out_thr.write(gene+'\t'+values[1]+'\t'+values[2]+'\t'+geneRPKM.get(GENE)+'\n')		 
				
				else:
					out.write(gene+'\t'+values[1]+'\t'+values[2].replace('\n','')+'\t'+'0.000'+'\n')
					#out.write(line.replace('\n','\t')+'0.000'+'\n')			
					#break

	out.close()
	out_thr.close()

########################################################################
