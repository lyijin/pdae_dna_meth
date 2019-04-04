#!/usr/bin/env python3

"""
> plot_clustermap.a_vs_s.median_meths.py <

Plots a clustered per-gene median methylation heatmap to illustrate the
difference in methylation in adults vs. sperm samples.
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# read genes of interest
genes_of_interest = open('t-test.a_vs_s.goi.tsv').read().strip().split('\n')

# read data
data = pd.read_table('../bias_density_medians/all.median_meths.tsv', index_col=0)

# subselect only adult and sperm samples
data = data.filter(regex='^A|^S')

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
output_filename = 'clustermap.a_vs_s.median_meths.pdf'
fig.savefig(output_filename, bbox_inches='tight')

# also save the values of the heatmap in a plain text file
cm.data2d.to_csv('clustermap.a_vs_s.median_meths.tsv', sep='\t')
