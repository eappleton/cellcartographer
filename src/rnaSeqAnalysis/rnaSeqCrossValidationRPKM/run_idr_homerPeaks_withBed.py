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
from matplotlib_venn import venn2

### functions ###
def makeOverlapHeatmap(groupSummaryPath, outputPath):
    '''
    Creates a heat mapimage given a path to a group summary file
    inputs: a path to a group summary file
    outputs: output directory path
    '''
    # create output directory if it doesn't exit
    if not os.path.isdir(outputPath) or not os.path.exists(outputPath):
        os.mkdir(outputPath)

    # read in input
    with open(groupSummaryPath) as f:
        data = f.readlines()
    # read in list of factors from first line of summary file
    factors = data[0].strip().split("\t")[4:]
    factorList = [] # list of all factors that appear in group summary
    factorIndexDict = {} # key: factor, value: order in which factors appear; starts from 0
    factorPeakDict = {} # key: factor, value: array of peak ids; which will be converted to a set of IDs
    counter = 0
    for factor in factors:
        factorList.append(factor)
        factorIndexDict[factor] = counter
        factorPeakDict[factor] = []
        counter += 1
    # read in data
    for line in data[1:]:
        tokens = line.strip().split("\t")
        factorsInPeak = tokens[1].split(",")
        peakID = tokens[2]

        # find which peaks belong to each factor
        for factor in factorsInPeak:
            if factor in factorPeakDict:
                factorPeakDict[factor].append(peakID)
            else:
                factorPeakDict[factor] = [peakID]
        
    # make sets out of the array of peak ids
    matrix = np.zeros((len(factors), len(factors)))
    for key in factorPeakDict.keys():
        factorPeakDict[key] = set(factorPeakDict[key])
        index = factorIndexDict[key]
        matrix[index][index] = len(factorPeakDict[key])
    factorCombos = []
    for i in range(len(factors) - 1):
        for j in range(i+1, len(factors)):
            factorCombos.append([factors[i], factors[j]])
        
    # create peak files for all triplets
    for factorCombination in factorCombos:
        factor1 = factorCombination[0]
        factor2 = factorCombination[1]
        index1 = factorIndexDict[factor1]
        index2 = factorIndexDict[factor2]
        intersection = len(factorPeakDict[factor1].intersection(factorPeakDict[factor2]))
        matrix[index1][index2] = intersection
        matrix[index2][index1] = intersection
#    matrix = np.log2(matrix)
#    matrix[np.isneginf(matrix)] = 0
    sns.heatmap(matrix,xticklabels=factors, yticklabels=factors)
    plt.title('Number of Intersecting Peaks')
    plt.savefig(outputPath + '/intersectionHeatmap.pdf',bbox_inches='tight')

def makeDoubleVennDiagrams(groupSummaryPath, outputPath, factorCombos):
    '''
    Creates a series of Venn diagram images that shows the overlap between 
        replicates and idr peaksgiven a group summary file
    inputs: a path to a group summary file
    outputs: output directory path
    '''
    # create output directory if it doesn't exit
    if not os.path.isdir(outputPath) or not os.path.exists(outputPath):
        os.mkdir(outputPath)

    # read in input
    with open(groupSummaryPath) as f:
        data = f.readlines()
    # read in list of factors from first line of summary file
    factors = data[0].strip().split("\t")[4:]
    factorList = [] # list of all factors that appear in group summary
    factorIndexDict = {} # key: factor, value: order in which factors appear; starts from 0
    counter = 0
    for factor in factors:
        factorList.append(factor)
        factorIndexDict[factor] = counter
        counter += 1
    factorPeakDict = {} # key: factor, value: array of peak ids; which will be converted to a set of IDs
    # read in data
    for line in data[1:]:
        tokens = line.strip().split("\t")
        factorsInPeak = tokens[1].split(",")
        peakID = tokens[2]

        # find which peaks belong to each factor
        for factor in factorsInPeak:
            if factor in factorPeakDict:
                factorPeakDict[factor].append(peakID)
            else:
                factorPeakDict[factor] = [peakID]
    # make sets out of the array of peak ids
    for key in factorPeakDict.keys():
        factorPeakDict[key] = set(factorPeakDict[key])

    # create peak files for all triplets
    for factorCombination in factorCombos:
        vennDiagramSets = []
        for factor in factorCombination:
            if factor in factorPeakDict:
                vennDiagramSets.append(factorPeakDict[factor])
            else:
                vennDiagramSets.append(set())
        vd = venn2(vennDiagramSets, set_labels = tuple(factorCombination))

        plt.savefig(outputPath+"/"+"_AND_".join(factorCombination)+".pdf", bbox_inches='tight')
        plt.close()

def read_homerPeak(peakFilePath):
    '''
    Reads a Homer peak file as a Pandas DataFrame  
    inputs: path to the Homer peak filen
    outputs: returns a Pandas data frame
    '''
    # determine how many lines to skip
    with open(peakFilePath) as f:
        data = f.readlines()
    numRowsToSkip = 0
    for line in data:
        if line[0] == '#':
            numRowsToSkip +=1
        else:
            break
    numRowsToSkip -= 1
    peakFrame = pd.read_csv(peakFilePath,sep='\t', skiprows=numRowsToSkip)
    peakFrame.name = peakFilePath.split("/")[-1].split(".")[0]
    peakFrame['Score'] = 0
    #peakFrame.columns = [x.split('(')[0].strip() for x in peakFrame.columns.values]
    return peakFrame

def convert_homerPeak_to_narrowPeak(homerPeakFrame, scoreColumn):
    # for running this script on its own output
    if not scoreColumn in homerPeakFrame.columns:
        scoreColumn = 'count'
    
    narrowPeakFrame= homerPeakFrame[['chr',
                                    'start',
                                    'end',
                                    '#PeakID',
                                    'Score',
                                    'strand',
                                    scoreColumn]]
    narrowPeakFrame['pValue'] = -1
    narrowPeakFrame['qValue'] = -1
    narrowPeakFrame['summit'] = [int(x) for x in ((narrowPeakFrame['end'] - narrowPeakFrame['start'])/2).values]
    #narrowPeakFrame['summit'] = -1
    narrowPeakFrame.name = homerPeakFrame.name
    return narrowPeakFrame

def convert_narrowPeak_to_bed(inputPath, outputPath, threshold):
    '''
    converts narrow peak files to Homer files
    '''
    idr_threshold_score = min([int(-125*np.log2(threshold)), 1000])

    with open(inputPath) as f:
        data = f.readlines()
    if len(data) > 0:
        narrowPeak_frame = pd.read_csv(inputPath, sep='\t')
        narrowPeak_frame = narrowPeak_frame.iloc[:,:12]
        narrowPeak_frame.columns = ['chrom', 'chromStart', 'chromEnd', '#PeakID', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak', 'local IDR', 'global IDR']
    else:
        narrowPeak_frame = pd.DataFrame(columns = ['chrom', 'chromStart', 'chromEnd', '#PeakID', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak', 'local IDR', 'global IDR'])

    narrowPeak_frame = narrowPeak_frame[['chrom', 'chromStart', 'chromEnd','#PeakID', 'score', 'strand']]
    narrowPeak_frame = narrowPeak_frame[narrowPeak_frame['score'] >= threshold]
    narrowPeak_frame.columns = ['#chr', 'start', 'end','PeakID', 'idrScore', 'strand']

    if '.' in narrowPeak_frame['PeakID'].values:
        rand = np.random.randint(0,1000)
        narrowPeak_frame['PeakID'] = [str(rand) + '_'+ str(x) for x in np.array(range(narrowPeak_frame.shape[0]))]
    if '.' in narrowPeak_frame['strand'].values:
        narrowPeak_frame['strand'] = '+'
    narrowPeak_frame.to_csv(outputPath, index= False, sep='\t')

def convert_narrowPeak_to_homerPeak(inputPath, outputPath, threshold):
    ''' 
    converts narrow peak files to Homer files
    '''
    idr_threshold_score = min([int(-125*np.log2(threshold)), 1000])

    with open(inputPath) as f:
        data = f.readlines()
    if len(data) > 0:
        narrowPeak_frame = pd.read_csv(inputPath, sep='\t')
        narrowPeak_frame = narrowPeak_frame.iloc[:,:12]
        narrowPeak_frame.columns = ['chrom', 'chromStart', 'chromEnd', '#PeakID', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak', 'local IDR', 'global IDR']
    else:
        narrowPeak_frame = pd.DataFrame(columns = ['chrom', 'chromStart', 'chromEnd', '#PeakID', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak', 'local IDR', 'global IDR'])

    narrowPeak_frame = narrowPeak_frame[['#PeakID','chrom', 'chromStart', 'chromEnd','strand', 'score', 'signalValue']]
    narrowPeak_frame = narrowPeak_frame[narrowPeak_frame['score'] >= threshold]
    narrowPeak_frame.columns = ['#PeakID','chr', 'start', 'end','strand', 'idrScore', 'count']

    if '.' in narrowPeak_frame['#PeakID'].values:
        rand = np.random.randint(0,1000)
        narrowPeak_frame['#PeakID'] = [str(rand) + '_'+ str(x) for x in np.array(range(narrowPeak_frame.shape[0]))]
    if '.' in narrowPeak_frame['strand'].values:
        narrowPeak_frame['strand'] = '+'
    narrowPeak_frame.to_csv(outputPath, index= False, sep='\t')


if __name__ == "__main__":
    execPath = os.path.dirname(os.path.realpath(__file__))
    # build argument parser
    parser = argparse.ArgumentParser(description='Given a set of Homer peak '+
        'files that are replicates, performs IDR analysis on all of them')
    parser.add_argument("samples",
        help="space separated listed of Homer peak files", nargs='+')
    parser.add_argument("outputPath",
        help="directory where output files should be written",
        default="~/", type=str)
    parser.add_argument("-threshold",
        help="idr threshold to use",
        default = "0.05", type=float)
    parser.add_argument("-scoreColumn",
        help="column to use for ranking peaks",
        default = ['findPeaks','Score'],
        type=str,
        nargs='+')
    parser.add_argument("-print", 
        help="just print commands", 
        default = False, action = "store_true")
    # parse arguments
    args = parser.parse_args()

    samples = args.samples
    outPath = args.outputPath
    threshold = args.threshold
    scoreColumn = ' '.join(args.scoreColumn)
    justPrint = False
    if args.print:
        justPrint = True
    
    if not os.path.exists(outPath):
        os.makedirs(outPath)

    print("Performing IDR analysis on the following samples:", ", ".join(samples))
    print("Output files will be written to:", outPath)
    print("Using the following IDR threshold:", threshold)
    print("Peaks will be ranked using:", scoreColumn)
    

    # read in peak files
    peak_frames = []
    toMerge = []
    for s in samples:
        pf = read_homerPeak(s)
        peak_frames.append(pf)
        toMerge.append(s)
    print("Other available scoreColumns:", peak_frames[0].columns.values[5:])

    # convert peak files to narrowPeak format for use with IDR
    for pf in peak_frames:
        npf = convert_homerPeak_to_narrowPeak(pf, scoreColumn)
        npf.to_csv(outPath + '/' + pf.name + '.narrowPeak', 
                   sep = '\t', 
                   index = False,
                   header = False)
        
    # run IDR between all pairs of replicates
    for i in range(len(peak_frames) - 1):
        for j in range(i + 1, len(peak_frames)):
            sample1 = peak_frames[i].name 
            sample2 = peak_frames[j].name
            resultPath = outPath + '/' + sample1 + '_' + sample2 + '_idr.out' 
            
            print('idr --samples ' +
                      outPath + '/' + sample1 + '.narrowPeak ' +
                      outPath + '/' + sample2 + '.narrowPeak ' + 
                      '--output-file ' + resultPath
                      + ' --plot'  
                      + ' --idr-threshold ' + str(threshold) + ' &'
                     )
            if not justPrint:
                if not os.path.isfile(resultPath):
                    os.system('idr --samples ' +
                              outPath + '/' + sample1 + '.narrowPeak ' +
                              outPath + '/' + sample2 + '.narrowPeak ' + 
                              '--output-file ' + resultPath
                              + ' --plot'  
                              + ' --idr-threshold ' + str(threshold)
                             )
                toMerge.append(resultPath.replace('out','tsv'))

            if not justPrint:
                # convert IDR output to Homer format
                convert_narrowPeak_to_bed(resultPath, 
                    resultPath.replace('out','bed'),
                    threshold) 
                convert_narrowPeak_to_homerPeak(resultPath,
                    resultPath.replace('out','tsv'),
                    threshold)
                



        
            
    
        
########################################################################
