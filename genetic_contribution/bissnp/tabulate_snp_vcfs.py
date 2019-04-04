#!/usr/bin/env python3

"""
> tabulate_snp_vcfs.py <

Compiles SNPs called by Bis-SNP from multiple files, then dumps genotype
information from all files into a giant tsv table.
"""
import argparse
import csv
import re

import numpy as np

import natural_sort

def convert_gt_to_int(gt_string):
    """
    This function converts GT calls from GATK into ints:
        0/0 (homozygous reference) --> 0
        0/1 (heterozygous reference) --> 1
        1/1, 1/2, ... (homozygous non-reference) --> 2
    
    Most homozygous non-reference bases are 1/1 (i.e. homozygous for first
    gt allele) anyway.
    """
    gt_ints = [int(x) for x in gt_string.split('/')]
    sum_gt_ints = sum(gt_ints)
    
    if sum_gt_ints > 2:
        sum_gt_ints = 2
    
    return sum_gt_ints

parser = argparse.ArgumentParser(description="""
Compiles SNPs called by Bis-SNP from multiple files, then dumps genotype
information from all files into a giant tsv table.""")

parser.add_argument('vcfs', metavar='vcf_files',
                    type=argparse.FileType('r'), nargs='+',
                    help='VCFs containing coverage information.')

args = parser.parse_args()

vcf_filenames = natural_sort.natural_sort([x.name for x in args.vcfs])

# read styl scaffold lengths to create appropriately-sized NumPy arrays
scaf_lens = {}

tsv_reader = csv.reader(open(
    '../../raw_data/pdae_genome.v1.scaffold_lengths.tsv'), delimiter='\t')
for row in tsv_reader:
    if not row: continue
    
    scaf_lens[row[0]] = int(row[1])

# chuck gt info into a dict containing many NumPy arrays
#   gt_info[file][scaf] = convert_gt_to_int(gt_bases)
gt_info = {}

for v in vcf_filenames:
    gt_info[v] = {}
    for s in scaf_lens:
        gt_info[v][s] = np.zeros(scaf_lens[s], dtype='int8')
    
    tsv_reader = csv.reader(open(v), delimiter='\t')
    for row in tsv_reader:
        if not row: continue
        if row[0][0] == '#': continue       # ignore commented lines
        
        scaf = row[0]
        pos = int(row[1]) - 1               # vcf files are 1-based
        gt_string = row[9].split(':')[0]
        gt_info[v][scaf][pos] = convert_gt_to_int(gt_string)

# print stuff out
print ('scaf', 'pos', *natural_sort.natural_sort(gt_info), sep='\t')
for s in natural_sort.natural_sort(scaf_lens):
    for n in range(scaf_lens[s]):
        pos_covs = [gt_info[x][s][n] for x in natural_sort.natural_sort(gt_info)]
        
        # do not print positions that are 0 coverage in all files
        if not any(pos_covs): continue
        
        # remember to convert 0-based positions back to 1-based positions
        results = [s, n + 1] + pos_covs
        print (*results, sep='\t')
