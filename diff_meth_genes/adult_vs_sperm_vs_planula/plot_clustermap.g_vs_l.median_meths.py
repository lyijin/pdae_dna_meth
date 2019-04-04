#!/usr/bin/env python3

"""
> plot_clustermap.g_vs_l.median_meths.py <

Plots a clustered per-gene median methylation heatmap to illustrate the
difference in methylation in gametic (S7, E8, S8) vs. larval samples.
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# read genes of interest
genes_of_interest = open('t-test.g_vs_l.goi.tsv').read().strip().split('\n')

# read data
data = pd.read_table('../bias_density_medians/all.median_meths.tsv', index_col=0)

# subselect only appropriate sperm and larval samples
data = data.filter(regex='^E8|^S7|^S8|^L')

# pick out subset containing genes of interest
data = data[data.index.isin(genes_of_interest)]

# plot the heatmap!
fig, ax = plt.subplots()
#sns.set(font_scale=1.2)

cm = sns.clustermap(data,
                    method='ward', z_score=0,
                    figsize=(6,8),
                    vmin=-2, vmax=2,
                    cmap='RdBu_r', annot=False, linewidths=0,
                    yticklabels=False,
                    cbar_kws={'ticks': [-2, -1, 0, 1, 2]})

# save figure
fig.tight_layout()
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'clustermap.g_vs_l.median_meths.pdf'
fig.savefig(output_filename, bbox_inches='tight')

# also save the values of the heatmap in a plain text file
cm.data2d.to_csv('clustermap.g_vs_l.median_meths.tsv', sep='\t')
