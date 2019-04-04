#!/usr/bin/env python3

"""
> check_pattern_passed_down.a_vs_s.py <

Some genes have methylation patterns that correlate well to larval survival
under heat stress. Are these positions inherited better from the respective
adult samples, when compared to genome-wide rates?

In numbers--if overall inheritance of A1-F --> S1-F has a tau of 0.65, if 
we subselected methylated positions in those interesting genes, would these
positions produce a tau > 0.65 (where tau = 1 means perfect concordance)?
"""
import numpy as np
import pandas as pd
import scipy.stats

# read survival data, and grab genes with changes in meth that strongly
# correlates with survival data (|r| > 0.8) 
heat_surv_data = pd.read_table('heat_surv_vs_meth.tsv', index_col=0)
all_genes = heat_surv_data.index
high_corr_genes = heat_surv_data[abs(heat_surv_data['pearson_r']) > 0.8].index

# read genomic context of methylated positions, then rename
# 'PdaeGeneX' --> 'PdaeX'
meth_pos_context_data = pd.read_table('../../A1-F.cov', usecols=[6],
                                      header=None, names=['gene'])
meth_pos_context_data['gene'] = \
    meth_pos_context_data['gene'].apply(lambda x: x.replace('PdaeGene', 'Pdae'))

# read table of meth values, subselect adult and sperm samples, then append
# gene context into the table
meth_data = pd.read_table('../../all.filt.pct.tsv', usecols=range(2,25))
meth_data = meth_data.filter(regex='^[A|S]')
meth_data['gene'] = meth_pos_context_data['gene']

# subselect meth_data rows corresponding to all genes used in earlier analysis
all_genes_data = meth_data[meth_data['gene'].isin(all_genes)]

# subselect meth_data rows that correspond to high_corr_genes
high_corr_genes_data = meth_data[meth_data['gene'].isin(high_corr_genes)]

# code is a bit inelegant, but meh, gets the job done.
# slice out individual colonies, then check kendall-tau for all meth pos vs.
# kendall-tau for meth pos in genes with high correlation to survival index
with open('check_pattern_passed_down.a_vs_s.tsv', 'w') as f:    
    # print header
    print ('Comparison', 'Tau (genome-wide)', 'Tau (all genes in analysis)',
           'Tau (genes with strong corr)', sep='\t', file=f)
    tau1 = []
    tau2 = []
    tau3 = []
    for n in range(1, 9):
        temp1 = meth_data.filter(regex=f'{n}-')
        temp2 = all_genes_data.filter(regex=f'{n}-')
        temp3 = high_corr_genes_data.filter(regex=f'{n}-')
        
        # crash code intentionally if either temps do not have 2 columns
        assert len(temp1.columns) == len(temp2.columns) == len(temp3.columns) == 2, \
            'One of the temp variables do not have two columns!'
        
        tau1.append(scipy.stats.kendalltau(temp1.iloc[:, 0], temp1.iloc[:, 1])[0])
        tau2.append(scipy.stats.kendalltau(temp2.iloc[:, 0], temp2.iloc[:, 1])[0])
        tau3.append(scipy.stats.kendalltau(temp3.iloc[:, 0], temp3.iloc[:, 1])[0])
        print (' vs. '.join(temp1.columns), round(tau1[-1], 5),
               round(tau2[-1], 5), round(tau3[-1], 5), sep='\t', file=f)
    
    # print means, SE and paired t-test p values
    print ('Mean', round(np.mean(tau1), 5), round(np.mean(tau2), 5),
           round(np.mean(tau3), 5), sep='\t', file=f)
    print ('Standard error', round(scipy.stats.sem(tau1), 5),
           round(scipy.stats.sem(tau2), 5), round(scipy.stats.sem(tau3), 5),
           sep='\t', file=f)
    print ('Paired t-test, genome-wide vs. genes with strong corr',
           scipy.stats.ttest_rel(tau1, tau3)[1], sep='\t', file=f)
    print ('Paired t-test, all genes vs. genes with strong corr',
           round(scipy.stats.ttest_rel(tau2, tau3)[1], 5), sep='\t', file=f)
    
    # print a few interesting stats
    print ('Methylated positions', len(meth_data), len(all_genes_data),
           len(high_corr_genes_data), sep='\t', file=f)
    print ('# genes', 'N/A', len(all_genes_data['gene'].unique()),
           len(high_corr_genes_data['gene'].unique()), sep='\t', file=f)
