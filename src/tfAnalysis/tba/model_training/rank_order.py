#!/usr/bin/env python
"""
Given a TBA coefficients file and TBA significance file produce a rank_ordered list of TFs to screen for differentiation
"""

### imports ###
import argparse
import numpy as np
import os    
import operator
    
if __name__ == '__main__':

	# How to parse arguments
	parser = argparse.ArgumentParser(description='Given a TBA coefficients file or a TBA significance file, maps the motif names to gene names')
	parser.add_argument("coeff_path", help="path to a TBA coefficients file", type=str)
	parser.add_argument("sig_path", help="path to a TBA significance file", type=str)
	parser.add_argument("output_path", help="path to output file", type=str)
    
	# Parse arguments
	args = parser.parse_args()

	coeff_path = args.coeff_path
	sig_path = args.sig_path
	output_path = args.output_path

	# Read in the significance and coefficients files 
	with open(coeff_path) as f:
		coeffs = f.readlines()[1:]
	with open(sig_path) as f:
                sig = f.readlines()[1:]

	# Create dictionaries with average of coeffs and discard the negative valued genes
	coeff_avg_dict = {}
	for line in coeffs:
		tf=line.replace('\n',"").split('\t')[1]
		all_coeffs_str=line.replace('\n',"").split('\t')[2:]
		all_coeffs=[float(i) for i in all_coeffs_str]
		avg=sum(all_coeffs)/len(all_coeffs)

		#Only consider positive values for now 
		if avg > 0:
			coeff_avg_dict[tf]=avg

	# Create dictionaries with average of coeffs and discard the negative valued genes
	sig_avg_dict = {}
	for line in sig:
		tf=line.replace('\n',"").split('\t')[1]
		all_sig_str=line.replace('\n',"").split('\t')[2:]
		all_sig=[float(i) for i in all_sig_str]
		avg=sum(all_sig)/len(all_sig)
		#sig_avg_dict[tf]=avg

		#Only consider positive values, i.e. non-zero significance 
		if avg > 0:
			if avg < 0.00316227766:
				sig_avg_dict[tf]=avg

	sorted_sig_avg_dict = sorted(sig_avg_dict.items(), key=operator.itemgetter(1))

	# Write results to output path
	out_file = open(output_path, 'w')
	out_file.write('GENE\tAVG_SIG\tAVG_COEF\n') #Header

	for idx, val in enumerate(sorted_sig_avg_dict):
		if val[0] in coeff_avg_dict:
			out_line=str(val[0])+'\t'+str(val[1])+'\t'+str(coeff_avg_dict[val[0]])+'\n'
			out_file.write(out_line)
	
	out_file.close()
