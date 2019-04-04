#!/usr/bin/env python3

"""
> t-test.a_vs_s.py <

Perform paired t-tests on a per-gene basis to detect significant changes in
median methylation levels between adult and sperm samples (A1-A8 vs. S1-S8).

Does not consider larval samples, as larval samples (L1-L6) does not have the
same genetic background as A1-A8 or S1-S8. The larval methylation levels are
much more similar to their parents (A7/S7/A8/S8).
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
adult = data.filter(regex='^A')
sperm = data.filter(regex='^S')

# while we're at it, check which groups have a delta meth of > 15% or < -15%
data['delta_15pc'] = abs(adult.mean(axis=1) - sperm.mean(axis=1)) > 15

# then perform t-test
results = scipy.stats.ttest_ind(adult.transpose(), sperm.transpose())
data['t-test'] = results.pvalue

# correct for multiple testing
x = pd.Series(correct_p_values.correct_p_values(data['t-test']),
              name='corr.t-test')
data = pd.concat([data, x], axis=1, sort=False)

# save output as two files: one for debug purposes, one for plotting purposes
data.to_csv('t-test.a_vs_s.output.tsv', sep='\t')

# save genes of interest sans details, just gene names
genes_of_interest = data[data['corr.t-test'] < 0.05][data['delta_15pc'] == True]
genes_of_interest.to_csv('t-test.a_vs_s.goi.tsv', columns=[], header=False)
