#!/usr/bin/env python3

"""
> plot_pairplot.fst_vs_meth.py <

Plots a pairplot using seaborn to show the correlation between all genetic
variables vs. all methylation variables.
"""
import math
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
import seaborn as sns

def pearson_r_text(x, y, **kws):
    # convert representation of p value to scientific if below 0.0001
    # for scientific: round number to nearest negative power of 10
    r, p_val = scipy.stats.pearsonr(x, y)
    if p_val < 0.0001:
        p_val = 10 ** (math.ceil(math.log10(p_val)))
        p_string = f'p < {p_val:.0e}'
    else:
        p_string = f'p = {p_val:.4f}'
    
    ax = plt.gca()
    if kws['label'] == 'Fujairah':
        ax.annotate(f'r = {r:.2f}\n{p_string}',
                    xy=(.95, .05), xycoords=ax.transAxes,
                    verticalalignment='bottom', horizontalalignment='right')
    else:
        ax.annotate(f'r = {r:.2f}\n{p_string}',
                    xy=(.95, .95), xycoords=ax.transAxes,
                    verticalalignment='top', horizontalalignment='right')
        

# read data, then filter out rows with < 10 SNPs within gene
# (i.e. valid rows have >= 5 meth positions and >= 10 SNPs in gene)
data = pd.read_table('delta_meth_vs_fst.gene.tsv', index_col=0)
data = data[data['n_snp_gene'] >= 10]

# create column to store which location has higher meth
data['higher_meth'] = data['delta_meth'].apply(
    lambda x: 'Fujairah' if x < 0 else 'Abu Dhabi')

# plot the pairplot!
sns.set(style='ticks', font_scale=1)
fig, ax = plt.subplots(figsize=(2, 8))

g = sns.pairplot(data, hue='higher_meth',
                 palette={'Fujairah': '#6baed6', 'Abu Dhabi': '#fcae91'},
                 kind='reg',
                 x_vars=['n_snp_gene', 'fst_gene', 'n_snp_exon', 'fst_exon'],
                 y_vars=['delta_meth'],
                 plot_kws={'line_kws':{'color':'#262626', 'alpha': 0.7},
                           'scatter_kws': {'s': 5, 'alpha': 0.2,
                                           'edgecolor': 'none'}})
g.map(pearson_r_text)

sns.despine(offset=2, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'pairplot.fst_vs_meth.loc_split.pdf'
fig.savefig(output_filename, bbox_inches='tight')
