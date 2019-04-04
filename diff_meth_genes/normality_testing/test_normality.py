#!/usr/bin/env python3

"""
> test_normality.py <

Perform a GLM analysis on a per-gene basis to detect significant changes in
mean methylation levels that can be attributed to developmental stage
(adult, sperm, larval) or geographical origin (Fujairah, Abu Dhabi), or
through the interaction of both variables!
"""
import pandas as pd
import scipy.stats

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
                                  'median_meth': []})
    for r in row.iteritems():
        loc = r[0].split('-')[1]
        devt = r[0][0]
        meth = float(r[1])
        temp = pd.DataFrame({'location': loc, 'devt_stage': devt,
                             'median_meth': meth}, index=[0])
        data_per_gene = pd.concat([data_per_gene, temp], ignore_index=True)
    
    # perform normality/equal variance tests
    a_meths = data_per_gene[data_per_gene['devt_stage'] == 'A']['median_meth']
    s_meths = data_per_gene[data_per_gene['devt_stage'] == 'S']['median_meth']
    l_meths = data_per_gene[data_per_gene['devt_stage'] == 'L']['median_meth']
    ad_meths = data_per_gene[data_per_gene['location'] == 'AD'][data_per_gene['devt_stage'] != 'L']['median_meth']
    f_meths = data_per_gene[data_per_gene['location'] == 'F']['median_meth']
    
    levene_devt_p = scipy.stats.levene(a_meths, s_meths, l_meths)[1]
    levene_location_p = scipy.stats.levene(ad_meths, f_meths)[1]
    
    shapiro_a_p = scipy.stats.shapiro(a_meths)[1]
    shapiro_s_p = scipy.stats.shapiro(s_meths)[1]
    shapiro_l_p = scipy.stats.shapiro(l_meths)[1]
    shapiro_ad_p = scipy.stats.shapiro(ad_meths)[1]
    shapiro_f_p = scipy.stats.shapiro(f_meths)[1]
    
    # hardcode parsing of results' p values
    temp = pd.DataFrame({'levene_devt': levene_devt_p,
                         'levene_location': levene_location_p,
                         'shapiro_a': shapiro_a_p,
                         'shapiro_s': shapiro_s_p,
                         'shapiro_l': shapiro_l_p,
                         'shapiro_ad': shapiro_ad_p,
                         'shapiro_f': shapiro_f_p}, index=[gene])
    
    pval_table = pd.concat([pval_table, temp])

# correct for multiple testing
for cols in pval_table.columns:
    x = pd.Series(correct_p_values.correct_p_values(pval_table[cols]),
                  name='corr.' + cols)
    pval_table = pd.concat([pval_table, x], axis=1, sort=False)

# save output as file
pval_table.to_csv('normality_tests.tsv', sep='\t')
