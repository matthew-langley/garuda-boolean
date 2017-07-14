# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:08:51 2017

@author: Matthew
"""

from sklearn.cluster import KMeans
import pandas as pd
import numpy as np
import os
import argparse


def discretize(inputFile, outputFile=None, inputFileType='csv', k=2, 
               invertOrder=False, perGene=True):
    """
    Converts a set of expression data (from microarray, qRT-PCR, etc.) into a 
    discretized data set, in which each expression value is classified into one 
    of a finite set of levels. Useful for binarizing data (ex. ON, OFF) or 
    preparing for (fuzzy) logic analysis (ex. HIGH, MED, LOW).
    
    Args:
        inputFile (str)
            Filepath for original expression data. In general, input file 
            should be in tabular form with rows:samples and columns:genes. 
            Expects first row to be gene names/IDs. Expects first column to be 
            sample names
        outputFile (str, optional)
            Filepath for discretized output. If not provided, data will only
            be returned to console and not saved to file
        inputFileType (str, optional)
            One of 'csv'
        k (int, optional)
            Used only for k-means. The number of discrete groups to generate
        invertOrder (bool, optional)
            If true, highest values in the original input file will be 
            assigned to the lowest-rank group (i.e. OFF, LOW, ...);  useful 
            when starting with Ct values, where greater Ct indicates lower 
            gene expression
        perGene (bool, optional)
            If true, discretization will be done on a per-gene basis; i.e. 
            expression values for Gene A will be discretized separately from 
            expression values for Gene B. Use if a consistent threshold 
            between groups (OFF vs. ON, for example) cannot be defined across 
            all genes.
            
    Returns:
        pd.DataFrame
            Table in same layout as input data, but with expression values 
            listed in discrete form. Discrete levels take the form of integers,
            with 0 representing the lowest expression level.
    """
    
    # Read input from file
    if inputFileType == 'csv':
        df = pd.read_csv(inputFile, index_col=0)
    else:
        raise ValueError('Unsupported file type: %s' %(inputFileType))
    
    # Perform discretization using kmeans
    if perGene is True:
        xList = []
        # results = pd.DataFrame()
        for gene in df.columns.tolist():
            x = __discretize_with_kmeans(df[gene], k, invert=invertOrder)
            xList.append(x)
        results = pd.concat(xList, axis=1)
    else:
        stack = df.stack()
        x = __discretize_with_kmeans(stack, k, invert=invertOrder)
        results = x.unstack()    
    
    if outputFile is not None:
        results.to_csv(outputFile)
        
    return results
            

def __discretize_with_kmeans(df, k, invert=False):
    
    # Do k-means clustering
    values = df.as_matrix().reshape((-1,1))
    km = KMeans(n_clusters=k)
    km.fit(values)
    labels0 = km.labels_
    labels = labels0.copy()
    centers = km.cluster_centers_.flatten()
    
    # Put groups in the right order
    if invert is True:
        sortorder = np.argsort(centers)[::-1]
    else:
        sortorder = np.argsort(centers)
    i = 0
    for j in sortorder:
        labels[labels0 == j] = i
        i += 1
        
    # Return the groups, indexed by the same as original data
    return pd.Series(labels, index=df.index, name=df.name)
    
    
if __name__ == '__main__':
    
    inputFile = ''
    outputFile = ''
    k = 2
    invertOrder = False
    perGene = False
    
    # Configure command line arguments
    parser = argparse.ArgumentParser(description='Discretize gene expression data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input',
                        help='Continuous gene expression data to be binarized, csv format, columns: genes, rows: samples')
    parser.add_argument('output',
                        nargs='?',
                        help='(optional) Destination file for discretized values')
    parser.add_argument('--k',
                        dest='k',
                        type=int,
                        default=2,
                        help='Number of levels')
    parser.add_argument('--invertOrder',
                        dest='invertOrder',
                        action='store_true',
                        default=False,
                        help='Invert ranking of discrete levels, i.e. high values --> 0, low values --> 1. Useful for binarizing Ct values from qRT-PCR')
    parser.add_argument('--perGene',
                        dest='perGene',
                        action='store_true',
                        default=False,
                        help='Set thresholds specific to each gene')                    

    # Switch between listening to command-line arguments or hard-coded 
    # arguments depending on whether running in IDE or from cmd.
    if any([name.startswith('SPYDER') for name in os.environ]):
        myArgs = 'test/single_cell_ct.csv --k 2 --invertOrder'
        args = parser.parse_args(myArgs.split())
    else:
        args = parser.parse_args()    
    
    if args.output is None:
        if '\\' in args.input:
            args.input = args.input.replace('\\', '/')
        args.output = '../Output/' + args.input.rsplit('/', 1)[-1].rsplit('.', 1)[0] + '_discretized.csv'
    
    results = discretize(inputFile=args.input, outputFile=args.output, k=args.k,
                         invertOrder=args.invertOrder, perGene=args.perGene)