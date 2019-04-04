#!/usr/bin/env python3

"""
> glm.mean_meths.py <

Perform a GLM analysis on a per-gene basis to detect significant changes in
mean methylation levels that can be attributed to developmental stage
(adult, sperm, larval) or geographical origin (Fujairah, Abu Dhabi), or
through the interaction of both variables!
"""
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

import correct_p_values

# read data
data = pd.read_table('../bias_density_medians/all.median_meths.tsv', index_col=0)

# retain genes with >= 5 methylated positions and is methylated
data = data[data['meth pos'] >= 5]
data = data.drop('meth pos', axis=1)
data = data[data.max(axis=1) > 0]

# drop E8-AD: not much one can do, statistically, with n=1
data = data.drop('E8-AD', axis=1)

# convert each row into a dataframe that can be used for subsequent GLM
pval_table = pd.DataFrame()
for gene, row in data.iterrows():
    data_per_gene = pd.DataFrame({'location': [], 'devt_stage': [], 
                                  'mean_meth': []})
    for r in row.iteritems():
        loc = r[0].split('-')[1]
        devt = r[0][0]
        meth = float(r[1])
        temp = pd.DataFrame({'location': loc, 'devt_stage': devt,
                             'mean_meth': meth}, index=[0])
        data_per_gene = pd.concat([data_per_gene, temp], ignore_index=True)
        
    model = smf.glm('mean_meth ~ devt_stage * location', data=data_per_gene,
                    family=sm.families.Gaussian())
    results = model.fit()
    
    # hardcode parsing of results' p values
    pvals = results.pvalues.values
    temp = pd.DataFrame({'devt_larval': pvals[1],
                         'devt_sperm': pvals[2],
                         'location': pvals[3],
                         'devt_sperm:location': pvals[5]}, index=[gene])
    
    pval_table = pd.concat([pval_table, temp])

# correct for multiple testing
x = pd.Series(correct_p_values.correct_p_values(pval_table['devt_larval']),
              name='corr.devt_larval')
pval_table = pd.concat([pval_table, x], axis=1, sort=False)

x = pd.Series(correct_p_values.correct_p_values(pval_table['devt_sperm']),
              name='corr.devt_sperm')
pval_table = pd.concat([pval_table, x], axis=1, sort=False)

x = pd.Series(correct_p_values.correct_p_values(pval_table['location']),
              name='corr.location')
pval_table = pd.concat([pval_table, x], axis=1, sort=False)

x = pd.Series(correct_p_values.correct_p_values(pval_table['devt_sperm:location']),
              name='corr.devt_sperm:location')
pval_table = pd.concat([pval_table, x], axis=1, sort=False)

# save output as file
pval_table.to_csv('glm_output.tsv', sep='\t')
