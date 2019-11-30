#!/usr/bin/env python3

"""
> plot_pairplot.allelic_diff_vs_meth.py <

Plots a pairplot using seaborn to show the correlation between all genetic
variables vs. all methylation variables; but instead of using Fst, use 
mean allelic frequency differences.
"""
import csv
import math
import statistics

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
    ax.annotate(f'r = {r:.2f}\n{p_string}',
                xy=(.95, .95), xycoords=ax.transAxes, 
                verticalalignment='top', horizontalalignment='right')

# read data, then filter out rows with < 10 SNPs within gene
# (i.e. valid rows have >= 5 meth positions and >= 10 SNPs in gene)
data = pd.read_table('delta_meth_vs_fst.gene.tsv', index_col=0)
data = data[data['n_snp_gene'] >= 10]

# convert delta_meth to absolute values (the non-abs values will be plotted 
# by a different script)
data['delta_meth'] = abs(data['delta_meth'])

# calculate per-gene mean allelic differences
gene_allelic_diff = {}
exon_allelic_diff = {}
tsv_reader = csv.reader(open('../bissnp/tabulated_genotypes.cov10-100.AS.annot.tsv'), 
                        delimiter='\t')
for row in tsv_reader:
    # ignore intergenic positions
    if 'PdaeGene' not in row[25]: continue
    
    # shorten annot
    gene = row[25].replace('PdaeGene', 'Pdae')
    if gene not in gene_allelic_diff:
        gene_allelic_diff[gene] = []
        exon_allelic_diff[gene] = []
    
    # sum allelic differences in A1-F to A4-F; S1-F to S4-F
    f_allelic_dist = sum([int(row[x]) for x in [2, 3, 4, 5, 17, 18, 19, 20]])
    
    # likewise for AD
    ad_allelic_dist = sum([int(row[x]) for x in [6, 7, 8, 9, 10, 11, 12, 13]])
    
    # calculate per-position mean allelic diff
    pos_allelic_dist = abs(f_allelic_dist - ad_allelic_dist) / 8
    gene_allelic_diff[gene].append(pos_allelic_dist)
    if 'Exon_' in row[29]:
        exon_allelic_diff[gene].append(pos_allelic_dist)

# go through dict, and calculate mean per-gene allelic diffs
gene_allelic_diff = {x:statistics.mean(y) for x, y in gene_allelic_diff.items()}
exon_allelic_diff = {x:(statistics.mean(y) if len(y) > 0 else 0) for x, y in exon_allelic_diff.items()}

# join these dicts into the original data dataframe
data = pd.concat([data, pd.Series(gene_allelic_diff, name='ad_gene')], axis=1, join='inner')
data = pd.concat([data, pd.Series(exon_allelic_diff, name='ad_exon')], axis=1, join='inner')

# plot the pairplot!
sns.set(style='ticks', font_scale=1)
fig, ax = plt.subplots(figsize=(4, 8))

g = sns.pairplot(data, kind='reg',
                 x_vars=['n_snp_gene', 'ad_gene', 'n_snp_exon', 'ad_exon'],
                 y_vars=['n_meth', 'delta_meth'],
                 plot_kws={'line_kws':{'color':'#262626', 'alpha': 0.7}, 
                           'scatter_kws': {'color': '#31a354', 's': 5, 
                                           'alpha': 0.2, 'edgecolor': 'none'}})
g.map(pearson_r_text)

sns.despine(offset=2, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'pairplot.allelic_diff_vs_meth.pdf'
fig.savefig(output_filename, bbox_inches='tight')
