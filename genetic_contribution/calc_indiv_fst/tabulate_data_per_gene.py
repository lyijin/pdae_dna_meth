#!/usr/bin/env python3

"""
> tabulate_data_per_gene.py <

Parse genetic and epigenetic data into per-gene bins, each row contains one bin
with all relevant data for future calculations.
"""
import csv

import pandas as pd

import natural_sort

def fst_ratio_of_avg(list_fst_n, list_fst_d):
    """
    Calculates and returns the ratio of averages of the Fsts in a specific
    window.
    """
    n = len(list_fst_d[list_fst_d > 0])
    if not n:
        return [0, 0]
    else:
        ratio_of_avg = sum(list_fst_n) / sum(list_fst_d)
        return [n, ratio_of_avg]

# read meth data
meth_data = pd.read_table(
    '../../diff_meth_genes/bias_density_medians/all.median_meths.tsv',
    index_col=0)

# retain genes with >= 5 methylated positions
meth_data = meth_data[meth_data['meth pos'] >= 5]

# subselect meth pos, adult and sperm samples
meth_data = meth_data.filter(regex='^meth pos|^A|^S')

# rename 'meth pos' to 'n_meth'
meth_data = meth_data.rename({'meth pos': 'n_meth'}, axis='columns')

# retain genes that are methylated in the subselected samples
meth_data = meth_data[meth_data.filter(regex='^A|^S').max(axis=1) > 0]

# use regexes to subselect groups of data
fuj = meth_data.filter(regex='-F$')
ad = meth_data.filter(regex='-AD$')

# calculate delta meth (as AD - F) and mean methylation of [fuj, ad]
meth_data['delta_meth'] = ad.mean(axis=1) - fuj.mean(axis=1)
meth_data['mean_meth'] = (fuj.mean(axis=1) + ad.mean(axis=1)) / 2

# read fst data
hudson_fsts = pd.read_table('hudson_fsts.tsv', index_col=0)

# for each gene, calculate ratio of averages
fst_data = {}
for g in set(hudson_fsts['gene']):
    fst_data[g] = {}
    
    # compute fst over entire gene
    hudson_gene = hudson_fsts[hudson_fsts['gene'] == g]
    fst_data[g]['n_snp_gene'], fst_data[g]['fst_gene'] = \
        fst_ratio_of_avg(hudson_gene['N'], hudson_gene['D'])
    
    # compute fst of exonic regions within gene
    hudson_exon = hudson_gene[hudson_gene['ei'] == 'exon']
    fst_data[g]['n_snp_exon'], fst_data[g]['fst_exon'] = \
        fst_ratio_of_avg(hudson_exon['N'], hudson_exon['D'])

fst_data = pd.DataFrame(fst_data).T
fst_data = fst_data[['n_snp_gene', 'fst_gene', 'n_snp_exon', 'fst_exon']]
fst_data = fst_data.astype({'n_snp_gene': int, 'fst_gene': float,
                            'n_snp_exon': int, 'fst_exon': float})

# combine methylation and fst data together--but first drop unnecessary
# columns in the methylation data
meth_data = meth_data.filter(regex='meth')
combined_data = meth_data.merge(fst_data, left_index=True, right_index=True)

# save combined data as file
combined_data.to_csv('delta_meth_vs_fst.gene.tsv', sep='\t', float_format='%.5f')
