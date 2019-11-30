#!/usr/bin/env python3

"""
> filter_adult_sperm_snps.py <

Script looks at "tabulated_genotypes.cov10-100.tsv", and for each position,
checks whether the SNP lies in an adult or sperm sample (i.e. A1-A8, S1-S8).

This is because the methylation-SNP analysis can only be done on the subset
of 16 samples, which has a balanced design (8 Abu Dhabi vs. 8 Fujairah samples;
each location has 4 adults and 4 sperm samples).
"""
import csv

# start parsing unfiltered table
tsv_reader = csv.reader(open('tabulated_genotypes.cov10-100.annot.tsv'),
                        delimiter='\t')

# print header
print (*next(tsv_reader), sep='\t')

for row in tsv_reader:
    if not row: continue
    
    a_genotypes = [int(x) for x in row[2:10]]
    s_genotypes = [int(x) for x in row[17:25]]
    
    # positions are discarded if genotypes of all adult/sperm samples from
    # F and AD are the same.
    if len(set(a_genotypes + s_genotypes)) == 1: continue
    
    # print stuff that passes this filter
    print (*row, sep='\t')
