#! /usr/bin/env python
################################################################################
'''
Given a list of peak files - merges them
'''

### header ###
__author__ = "Jenhan Tao"
__license__ = "BSD"
__email__ = "jenhantao@gmail.com"


### imports ###
import sys
import os
import pandas as pd
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import seaborn as sns
from itertools import combinations
from itertools import product
from collections import Counter

### functions ###
def split_rpkm_file(rpkmFilePath, conditions, outPath):
    '''
    Reads a Homer rpkm file as a Pandas DataFrame  
    inputs: path to the Homer peak file, condition represented by each column,
        and output path
    outputs: writes a series of narrow peak files for each column, 
        returns a data frame giving gene information and a list containing
        tuples of replicates
    '''
    condition_frequency_dict = {x:1 for x in set(conditions)}
    condition_tuple_dict = {x:[] for x in set(conditions)}
    rpkm_frame = pd.read_csv(rpkmFilePath, sep = "\t")
    columns = list(rpkm_frame.columns.values)
    for i in range(len(conditions)):
        currentCondition = conditions[i]
        newColName = currentCondition + '_' + \
            str(condition_frequency_dict[currentCondition])
        condition_tuple_dict[currentCondition].append(newColName)    
        condition_frequency_dict[currentCondition] += 1
        columns[8 + i] = newColName

    condition_tuple_dict = {x:tuple(condition_tuple_dict[x]) for x in \
        condition_tuple_dict.keys()}
    columns[0] = 'refseq'
    rpkm_frame.columns = columns
    rpkm_frame['score'] = 0
    narrowPeakColumns = ['chr', 'start', 'end', 'refseq', 'score', 'strand'] 
    for col in columns[8:]:
        currentColumns = narrowPeakColumns + [col]
        narrowPeak_frame = rpkm_frame[currentColumns]
        narrowPeak_frame['pValue'] = -1
        narrowPeak_frame['qValue'] = -1
        narrowPeak_frame['peak'] = [int(x) for x in ((narrowPeak_frame['end'] \
            - narrowPeak_frame['start'])/2).values]
        narrowPeak_frame.to_csv(outPath + '/' + col + '.narrowPeak', 
            sep = '\t', index = False, header = False)
    return rpkm_frame[['refseq', 
                'Annotation/Divergence', 
                'chr', 
                'start', 
                'end',
                'strand']], condition_tuple_dict
    
        

def merge_idrFrames(idrOut_frames, outPath):
    # merge results from the replicates together; 
    # calculate the union and the intersection
    intersection_genes = set(idrOut_frames[0]['refseq'].values)
    union_genes = set(idrOut_frames[0]['refseq'].values)
    for of in idrOut_frames[1:]:
       intersection_genes = \
           intersection_genes.intersection(set(of['refseq'].values))
       union_genes = \
           union_genes.union(set(of['refseq'].values))
    
    intersection_frame = \
        idrOut_frames[0][idrOut_frames[0]['refseq'].isin(intersection_genes)]
    union_frame = \
        idrOut_frames[0][idrOut_frames[0]['refseq'].isin(union_genes)]
    
    # merge data frames
    count =0
    for of in idrOut_frames[1:]:
        count+=1
        intersection_frame = intersection_frame.merge(of, on=['refseq', 
            'chr', 
            'start', 
            'end', 
            'strand',
            'annotation'])
        intersection_frame = \
            intersection_frame.set_index(intersection_frame['refseq'].values)

        union_frame = union_frame.merge(of, how='outer', on=['refseq', 
            'chr', 
            'start', 
            'end', 
            'strand',
            'annotation'])

        union_frame.dropna(subset=['refseq'],inplace=True)

        union_frame = union_frame.set_index(union_frame['refseq'].values)

        union_frame = union_frame.fillna('0')

    # concatenate idrScore and count columns
    # concatenate idrScore columns

    # for the intersection
    idrScore_intersection_indices = []
    idrScore_union_indices = []
    intersection_columns = intersection_frame.columns.values
    for i in range(len(intersection_columns)):
        if 'idrScore' in intersection_columns[i]:
            idrScore_intersection_indices.append(i)
    idrScoreColumn_intersection = intersection_frame.ix[:,idrScore_intersection_indices[0]].values
    for ind in idrScore_intersection_indices[1:]:
        idrScoreColumn_intersection = idrScoreColumn_intersection + ', ' \
            + intersection_frame.ix[:,ind].values
    intersection_frame['idrScore'] = idrScoreColumn_intersection

    intersection_frame = intersection_frame[[
        'chr', 
        'start', 
        'end',
        'refseq',
        'annotation',
        'idrScore']]
    intersection_frame.to_csv(outPath + '_idrIntersection.tsv', index=False, sep="\t")

    # for the union
    idrScore_union_indices = []
    idrScore_union_indices = []
    union_columns = union_frame.columns.values
    for i in range(len(union_columns)):
        if 'idrScore' in union_columns[i]:
            idrScore_union_indices.append(i)
    idrScoreColumn_union = union_frame.ix[:,idrScore_union_indices[0]].values
    for ind in idrScore_union_indices[1:]:
        idrScoreColumn_union = idrScoreColumn_union + ', ' \
            + union_frame.ix[:,ind].values
    union_frame['idrScore'] = idrScoreColumn_union

    union_frame = union_frame[[
        'chr', 
        'start', 
        'end',
        'refseq',
        'annotation',
        'idrScore']]
    union_frame.to_csv(outPath + '_idrUnion.tsv', index=True, sep="\t")

    return set(list(intersection_frame['refseq'].values)), set(list(union_frame['refseq'].values))

def run_idr(sample1, sample2, outPath):
    resultPath = outPath + '/' + sample1 + '_' + sample2+ '_idr.out'
    if not os.path.isfile(resultPath):
        os.system('idr --samples ' +
                  outPath + '/' + sample1 + '.narrowPeak ' +
                  outPath + '/' + sample2 + '.narrowPeak ' + 
                  '--output-file ' + resultPath
                  + ' --plot'  
                  + ' --idr-threshold ' + str(threshold)
                  + ' --peak-list ' +outPath + '/' + sample2 + '.narrowPeak '
                 )

    with open(resultPath) as f:
        lines = f.readlines()
    chromosomes = []
    starts = []
    ends = []
    scores = []
    strands = []
    counts = [] # rpkms or tag counts
    refseqs = []
    annotations = []
    for line in lines:
        tokens = line.strip().split("\t")
        chromosomes.append(tokens[0])
        starts.append(tokens[1])
        ends.append(tokens[2])
        gene = position_gene_dict[(tokens[0], 
            int(tokens[1]), 
            int(tokens[2]), 
            tokens[5])]
        refseqs.append(gene[0])
        annotations.append(gene[1])
        scores.append(tokens[4])
        strands.append(tokens[5])

    data = {'refseq':refseqs,
            'annotation':annotations,
            'chr':chromosomes,
            'start':starts,
            'end':ends,
            'strand':strands,
            'idrScore':scores
            } 
    outFrame = pd.DataFrame(data)
    outFrame = outFrame.set_index(outFrame['refseq'])
    outFrame = outFrame[[
            'chr',
            'start',
            'refseq',
            'annotation',
            'end',
            'strand',
            'idrScore']]
    outFrame = outFrame.drop_duplicates(['chr',
        'start',
        'end', 
        'strand'])
    outFrame.to_csv(resultPath.replace('.out','.tsv'), 
        index=False, sep='\t')
    return outFrame

if __name__ == "__main__":
    execPath = os.path.dirname(os.path.realpath(__file__))
    # build argument parser
    parser = argparse.ArgumentParser(description='Given a Homer rpkm file'+
        ' and a list of experimental conditions corresponding to each column'+
        ' returns a list of peaks that passes the IDR threshold')
    parser.add_argument("rpkm",
        help="path to a Homer rpkm file")
    parser.add_argument("outputPath",
        help="directory where output files should be written",
        default="~/", type=str)
    parser.add_argument('conditions',
        help='list of experimental conditions corresponding to each column of' +
        ' RPKM file', nargs='+')
    parser.add_argument("-threshold",
        help="idr threshold to use",
        default = "0.05", type=float)
    parser.add_argument("-print", 
        help="just print commands", 
        default = False, action = "store_true")
    # parse arguments
    args = parser.parse_args()

    rpkmFilePath = args.rpkm
    outPath = args.outputPath
    conditions = args.conditions
    threshold = args.threshold
    justPrint = False
    if args.print:
        justPrint= True
    
    if not os.path.exists(outPath):
        os.makedirs(outPath)
    with open(rpkmFilePath) as f:
        data = f.readlines()
    columns = data[0].strip().split("\t")[8:]
   
    if len(columns) != len(conditions):
        print("The number of conditions (" + str(len(conditions)) + ")" + \
            " must match the number of columns in the rpkm file (" + \
            str(len(columns)) + ")" )
        sys.exit(1)
    
    print("Output files will be written to:", outPath)
    print("Using the following IDR threshold:", threshold)
    print("Performing IDR analysis using the following columns and conditions:")
    for i in range(len(columns)):
        print(columns[i], "\t", conditions[i])
    
    # split rpkm file into narrow peak files
    gene_frame, condition_replicate_dict = \
        split_rpkm_file(rpkmFilePath, conditions, outPath)
    gene_positions = zip(gene_frame['chr'], 
        gene_frame['start'], 
        gene_frame['end'], 
        gene_frame['strand'])
    gene_refseqAnnotation = zip(gene_frame['refseq'], gene_frame['Annotation/Divergence'])
    position_gene_dict = dict(zip(gene_positions, gene_refseqAnnotation))
    
    # run IDR between all pairs of replicates
    intersection_refseqs = None
    union_refseqs = set()

    for condition in condition_replicate_dict.keys():
        replicate_tuple = condition_replicate_dict[condition]
        idrOut_frames = []
        for i in range(len(replicate_tuple) - 1):
            for j in range(i + 1, len(replicate_tuple)):
                sample1 = replicate_tuple[i]
                sample2 = replicate_tuple[j]

                if justPrint:
                    resultPath = outPath + '/' + sample1 + '_' + sample2+ '_idr.out'
                    print('idr --samples ' +
                              outPath + '/' + sample1 + '.narrowPeak ' +
                              outPath + '/' + sample2 + '.narrowPeak ' + 
                              '--output-file ' + resultPath
                              + ' --plot'  
                              + ' --idr-threshold ' + str(threshold)
                              + ' --peak-list ' +outPath + '/' + sample2 + '.narrowPeak &'
                             )
                else:

                    outFrame = run_idr(sample1, sample2, outPath)
                    if outFrame.shape[0] > 0:
                        idrOut_frames.append(outFrame)
        if not justPrint:
        
            if len(idrOut_frames) > 0:
                # report IDR scores
                intersection_refs, union_refs = merge_idrFrames(idrOut_frames,
                                                         outPath + '/' + condition)
            else:
                intersection_refs = set()
                union_refs = set()

            if intersection_refseqs != None:
                intersection_refseqs = intersection_refseqs.intersection(intersection_refs)
            else:
                if intersection_refs != None:
                    intersection_refseqs = intersection_refs.copy()

            union_refseqs = union_refseqs.union(union_refs)

    # produce rpkm spread sheet genes that pass through IDR only
    if not justPrint:
        rpkm_frame = pd.read_csv(rpkmFilePath, sep='\t')
        rpkm_frame[rpkm_frame.ix[:,0].isin(intersection_refseqs)].to_csv(outPath + '/idr_intersection_rpkm.tsv', sep='\t', index=False)
        rpkm_frame[rpkm_frame.ix[:,0].isin(union_refseqs)].to_csv(outPath + '/idr_union_rpkm.tsv', sep='\t', index=False)
    

    
        
########################################################################
