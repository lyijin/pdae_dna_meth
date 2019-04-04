#!/usr/bin/env python3

"""
> plot_correlation.meth.py <

Plots a correlation heatmap using seaborn.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# read in methylation %s: [0-100]
data = pd.read_table('../../all.filt.pct.tsv', usecols=range(2,25))
data_corr = data.corr(method='kendall')

# save values in the correlation matrix as a file
data_corr.to_csv('corr_plot.meth.tsv', sep='\t')

mask = np.zeros((len(data_corr), len(data_corr)), 'int8')
np.fill_diagonal(mask, 1)

# plot the heatmap!
fig, ax = plt.subplots(figsize=(6, 6))
sns.set(font_scale=1.3)

cg = sns.clustermap(data_corr, method='complete',
                    vmin=0.40, vmax=0.70, mask=mask,
                    cmap='YlGn', annot=False, linewidth=0.5,
                    cbar_kws={'ticks': [0.4, 0.5, 0.6, 0.7]})

# save figure
fig.tight_layout()
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'corr_plot.meth.pdf'
fig.savefig(output_filename, bbox_inches='tight')
