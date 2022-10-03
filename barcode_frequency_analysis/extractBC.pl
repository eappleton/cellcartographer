#!/usr/bin/perl

use strict;
use warnings;

#####
# Extract a sequence in relation to a common sequence in a FASTQ file
#####

## Author: Alex Ng
## Version 2.4
## Date: 2019 Jan 4

#####
# History of updates
#####

# 2.4 - for use on raw RNA-seq. Looks for common sequence on both sides of barcode (both strands). Note: it requires a perfect match. An improvement would be to use fuzzy matching or Hamming distance-based matching, but I suspect it will dramatically increase runtime. CSV file only contains the barcode sequence.
# v2.3 - for use on a FASTQ file
# v2.2 - CSV files has three columns now: cell barcode, gene name, and molecular barcode
# v2.1 - also outputs a CSV file with two columns: cell barcode and gene name. This file may be used to subtract counts from the 10x-generated count matrix
# v2.0 - instead of outputting the entire line, it converts the output to FASTQ, keeping the sequence name, sequence, quality, and moves the 10x Chromium CB tag (correct cellular barcode) to the sequence name

# Example command:
# Input: FASTQ file from a 10x scRNAseq TF-amplified run
# Note: it will only find one instance of it per line (Cell Ranger SAM files are one sequence per line)
# Output:
# A CSV file with one column: the TF barcode (20bp)
# Example run:
# perl 20180724-extract-TFbarcode-from-scRNAseq.pl little.fastq

#####
# Main code
#####

# Hard-coded AttL sites to look for
# Approximately 25-mers of the common region between TF and polyA on bcPBAN, forward strand
my $fwdseq1 = "CCAAGCACCTGCTACATAGC";

# Approximately 25-mers of the common region between TF and polyA on bcPBAN, reverse strand (other side of barcode, still 5' to barcode, on reverse strand)
my $revseq1 = "GGCACAGTCGAGGGGTACCG";

# Open FASTQ file
my $inputfile = $ARGV[0];
open(INPUTFILE,$inputfile);

# Open output files
open(OUTPUTCSV,">$inputfile.csv") or die "Couldn't open $inputfile.csv: $!";

# Use every fourth line, starting on line 2 (to grab the sequencing read, and not the header, "+" or Q30 score)
my $useline = -1;

while(my $thisline = <INPUTFILE>)
{
	# Remove returns
	chomp($thisline);

	# Uses simple assumption that this FASTQ file starts with @ header, and 2nd line (and every 4th onwards) contains the sequencing read.
	#### If it's a line that has a name, print the name and go to next line
	# Extract the line after the header line
#	if ($thisline =~ /^@/) {

	# Use a line counter to determine which line to grab
	if($useline % 4) # notice $useline is set to -1
	{	
		# not divisible
		$useline++;
		next;
	} else {
		# divisible
		$useline++;
		
		# If the common sequence exists in this read, continue, else skip
		if ($thisline =~ /$fwdseq1/) {
		
			# Find position where barcode is, by looking for known sequence
			my $positionknownsequence = index($thisline,$fwdseq1);
			my $lengthknownsequence = length($fwdseq1);
			my $barcode = substr($thisline,$positionknownsequence+$lengthknownsequence,$lengthknownsequence);
	
		print OUTPUTCSV "$barcode\n";
		}

		# Also consider the scenario that the sequencing read is on the bottom strand
		if ($thisline =~ /$revseq1/) {
	
		# Find position where barcode is, by looking for known sequence
		my $positionknownsequence = index($thisline,$revseq1);
		my $lengthknownsequence = length($revseq1);
		my $barcode = substr($thisline,$positionknownsequence+$lengthknownsequence,$lengthknownsequence);
		
		# IMPORTANT: reverse complements it so the barcode is output as top strand
		$barcode = reverse $barcode;
		$barcode =~ tr/ATGCatgc/TACGtacg/;
		print OUTPUTCSV "$barcode\n";
	
		}

	}
}	

close INPUTFILE;
#close OUTPUTFASTQ;
close OUTPUTCSV;

