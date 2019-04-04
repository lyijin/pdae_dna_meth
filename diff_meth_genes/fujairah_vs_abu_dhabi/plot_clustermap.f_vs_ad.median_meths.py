#!/usr/bin/env python3

"""
> plot_clustermap.f_vs_ad.median_meths.py <

Plots a clustered per-gene median methylation heatmap to illustrate the
difference in methylation in Fujairah vs. Abu Dhabi samples.
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# read genes of interest
genes_of_interest = open('t-test.goi.tsv').read().strip().split('\n')

# read data
data = pd.read_table('../bias_density_medians/all.median_meths.tsv', index_col=0)

# subselect only adult and sperm samples
data = data.filter(regex='^A|^S')

# pick out subset containing genes of interest
data = data[data.index.isin(genes_of_interest)]

# plot the heatmap!
fig, ax = plt.subplots()
#sns.set(font_scale=1.2)

cg = sns.clustermap(data,
                    method='ward', z_score=0,
                    figsize=(5,6),
                    vmin=-2, vmax=2,
                    cmap='RdBu_r', annot=False, linewidth=0,
                    yticklabels=False,
                    cbar_kws={'ticks': [-2, -1, 0, 1, 2]})

# save figure
fig.tight_layout()
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'clustermap.f_vs_ad.median_meths.pdf'
fig.savefig(output_filename, bbox_inches='tight')

# also save the values of the heatmap in a plain text file
cg.data2d.to_csv('clustermap.f_vs_ad.median_meths.tsv', sep='\t')
