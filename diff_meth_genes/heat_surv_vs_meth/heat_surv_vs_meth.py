#!/usr/bin/env python3

"""
> heat_surv_vs_meth.py <

Perform a linear correlation a per-gene basis to check which genes show a
strong correlation between median methylation levels and larval survival rates.

Code assumes Python >= 3.7, where dicts have preserved order.
"""
import pandas as pd
import scipy.stats

# hardcode survival index values derived in excel spreadsheet
surv_index = {'S1-F': 0.789439953,
              'S3-F': 0.998209258,
              'S4-F': 0.22654263,
              'S5-AD': 0.954202839,
              'S6-AD': 1.477754999,
              'S8-AD': 1.300719112}
surv_keys = list(surv_index.keys())
surv_values = list(surv_index.values())

# read data
data = pd.read_table('../bias_density_medians/all.median_meths.tsv', index_col=0)

# retain genes with >= 5 methylated positions
data = data[data['meth pos'] >= 5]

# retain the six columns with survival data
data = data[surv_keys]

# retain genes that are methylated in all six samples
data = data[data.min(axis=1) > 0]

# carry out correlation of meth values and surv values
r = data.apply(lambda x: scipy.stats.pearsonr(x, surv_values)[0], axis=1)
p = data.apply(lambda x: scipy.stats.pearsonr(x, surv_values)[1], axis=1)
data['pearson_r'] = r
data['p_val'] = p

# save output as file
data.to_csv('heat_surv_vs_meth.tsv', sep='\t')
