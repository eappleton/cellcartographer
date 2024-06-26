#!/usr/bin/env python
"""
Given a set of genomic coordinates in BED format:
chr start end
...

Extracts the genomic sequences
"""

### imports ###
import sys
import os
import inspect

def read_bed_file(input_path):
    '''
    reads a bed file and returns the genomic coordinates
    '''
    with open(input_path) as f:
        data = f.readlines()
    coordinates = []
    if data[0].strip()[0] == '#':
        data = data[1:]
    for line in data:
        tokens = line.strip().split()
        chrom = tokens[0]
        start = tokens[1]
        end = tokens[2]
        coordinates.append((chrom,start, end))
    return coordinates

def extract_sequence(coordinates, genome, out_file_path):
    '''
    Given a list of genomic coordinates, extracts sequences
    inputs: [(chrom1, start1, end1), ..., (chromN, startN, endN)]
    outputs: [seq1, seq2, ...seqN]
    '''
    script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    if genome == 'mm10':
        genome_path = script_path + '/mm10/'
        chromosomes = ['chr1' , 'chr2' , 'chr3' , 'chr4' , 'chr5' , 
                        'chr6' , 'chr7' , 'chr8' , 'chr9' , 'chr10', 
                        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
                        'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
    elif genome == 'hg38':
        genome_path = script_path + '/hg38/'
        chromosomes = ['chr1' , 'chr2' , 'chr3' , 'chr4' , 'chr5' , 
                        'chr6' , 'chr7' , 'chr8' , 'chr9' , 'chr10', 
                        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
                        'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 
                        'chr21', 'chr22', 'chrX']
    chrom_size_dict = {}
    chrom_seq_dict = {}

    print('reading genome', genome)
    for chrom in chromosomes:
        with open(genome_path + chrom + '.fa') as f:
            data = f.readlines()
        seq = ''.join(x.upper().strip() for x in data[1:])
        size = len(seq)
        chrom_size_dict[chrom] = size
        chrom_seq_dict[chrom] = seq 
    
    out_file = open(out_file_path, 'w')
    for coord in coordinates:
        chrom = coord[0]
        # chrom_seq dict is 0 indexed, genome coords are 1 indexed
        start = int(coord[1]) - 1
        end = int(coord[2]) - 1
        if chrom in chrom_seq_dict:
            seq = chrom_seq_dict[chrom][start:end] 
            id_line = '>' + str(coord[0]) + ':' +str(coord[1]) + '-' + str(coord[2]) + '\n'
            out_file.write(id_line)
            out_file.write(seq + '\n')
    out_file.close()

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage')
        print('extract_sequences.py <bed file> <genome> <output_file_path>')
        sys.exit(0)
    else:
        bed_path = sys.argv[1]
        genome = sys.argv[2]
        output_path = sys.argv[3]

        coordinates = read_bed_file(bed_path)
        extract_sequence(coordinates, genome, output_path)
    
