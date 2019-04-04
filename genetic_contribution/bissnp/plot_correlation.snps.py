#!/usr/bin/env python3

"""
> plot_correlation.snps.py <

Plots a heatmap using seaborn to show pairwise cityblock distance (i.e. overall 
genetic distance) between all samples.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.spatial.distance
import seaborn as sns

# remove gene names, and filter out rows containing all 0s
data = pd.read_table('tabulated_genotypes.cov10-100.tsv', usecols=range(2,25))

# calculate pairwise city-block distances
data_corr = scipy.spatial.distance.squareform(
    scipy.spatial.distance.pdist(data.transpose(), 'cityblock'))
data_corr = pd.DataFrame(data_corr, 
    index=data.columns, columns=data.columns)

# dump values into a table for supp data
data_corr.to_csv('correlation.snps.tsv', sep='\t')

mask = np.zeros((len(data_corr), len(data_corr)), 'int8')
np.fill_diagonal(mask, 1)

# plot the heatmap!
fig, ax = plt.subplots(figsize=(6, 6))
sns.set(font_scale=1.2)

cg = sns.clustermap(data_corr, mask=mask,
                    cmap='YlGn_r', annot=False, linewidth=0.5)
plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

# save figure
fig.tight_layout()
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'correlation.snps.pdf'
fig.savefig(output_filename, bbox_inches='tight')
