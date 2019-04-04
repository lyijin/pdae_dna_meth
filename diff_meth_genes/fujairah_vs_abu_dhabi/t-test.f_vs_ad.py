#!/usr/bin/env python3

"""
> t-test.f_vs_ad.py <

Perform paired t-tests on a per-gene basis to detect significant changes in
median methylation levels due to location (Fujairah vs. Abu Dhabi).

Only considers sperm and adult samples--the inclusion of planula samples in
Abu Dhabi causes an upward bias in the median methylation.
"""
import pandas as pd
import scipy.stats

import correct_p_values

# read data
data = pd.read_table('../bias_density_medians/all.median_meths.tsv', index_col=0)

# retain genes with >= 5 methylated positions
data = data[data['meth pos'] >= 5]
data = data.drop('meth pos', axis=1)

# subselect only adult and sperm samples
data = data.filter(regex='^A|^S')

# retain genes that are methylated in the subselected samples
data = data[data.max(axis=1) > 0]

# use regexes to subselect groups of data
fuj = data.filter(regex='-F$')
ad = data.filter(regex='-AD$')

# while we're at it, check which groups have a delta meth of > 15% or < -15%
data['delta_15pc'] = abs(fuj.mean(axis=1) - ad.mean(axis=1)) > 15

# then perform t-test
results = scipy.stats.ttest_ind(fuj.transpose(), ad.transpose())
data['t-test'] = results.pvalue

# correct for multiple testing
x = pd.Series(correct_p_values.correct_p_values(data['t-test']),
              name='corr.t-test')
data = pd.concat([data, x], axis=1, sort=False)

# save output as two files: one for debug purposes, one for plotting purposes
data.to_csv('t-test.output.tsv', sep='\t')

# save genes of interest sans details, just gene names
genes_of_interest = data[data['corr.t-test'] < 0.05][data['delta_15pc'] == True]
genes_of_interest.to_csv('t-test.goi.tsv', columns=[], header=False)
