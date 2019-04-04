#!/usr/bin/env python3

"""
> plot_2d_pca.py <

Plot PCA of all methylation levels from all samples.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

def assign_samples(label):
    """
    Assign sample types based on the label.
    """
    if label[-2:] == '-F':
        sample_type = 'Fujairah '
    else:
        sample_type = 'Abu Dhabi '
    
    if label[0] == 'A':
        sample_type += 'Adult'
    elif label[0] == 'E':
        sample_type += 'Egg'
    elif label[0] == 'S':
        sample_type += 'Sperm'
    elif label[0] == 'L':
        sample_type += 'Larva'
    
    return sample_type

# read in methylation %s: [0-100]
data = pd.read_table('../../all.filt.pct.tsv', usecols=range(2,25))

# PCA stuff
pca = PCA(n_components=3)
pca_array = pca.fit_transform(data.transpose())
pca_evr = [round(x * 100, 1) for x in pca.explained_variance_ratio_]
    
# plot PCA
fig, ax_array = plt.subplots(1, 3, figsize=(18, 6))
point_labels = data.columns.tolist()

labels = [assign_samples(x) for x in point_labels]

colours = {'Fujairah Adult': '#08519c', 'Fujairah Sperm': '#6baed6',
           'Abu Dhabi Adult': '#a50f15', 'Abu Dhabi Sperm': '#fcae91',
           'Abu Dhabi Egg': '#fb6a4a', 'Abu Dhabi Larva': '#de2d26'}
markers = {'Fujairah Adult': 's', 'Fujairah Sperm': 'o',
           'Abu Dhabi Adult': 's', 'Abu Dhabi Sperm': 'o',
           'Abu Dhabi Egg': '*', 'Abu Dhabi Larva': '^'}

pc1_values = [x[0] for x in pca_array]
offset = 0.02 * (max(pc1_values) - min(pc1_values))

for m, ax in enumerate(ax_array):
    x_pc = [0, 0, 1]
    y_pc = [1, 2, 2]
    for l in sorted(set(labels)):
        x = []
        y = []
        for n in range(len(pca_array)):
            if labels[n] == l:
                x.append(pca_array[n][x_pc[m]])
                y.append(pca_array[n][y_pc[m]])
                ax.text(pca_array[n][x_pc[m]] + offset, pca_array[n][y_pc[m]],
                        point_labels[n][:2], ha='left', va='center')
        
        ax.scatter(x, y, c=colours[l], marker=markers[l], label=l)
    
    # change ordering of legend
    old_handles, old_labels = ax.get_legend_handles_labels()
    new_order = [4, 5, 0, 3, 1, 2]
    ax.legend([old_handles[x] for x in new_order],
              [old_labels[x] for x in new_order],
              bbox_to_anchor=(0,1.02,1,0.2),
              loc='lower left',
              mode='expand',
              borderaxespad=0,
              ncol=3)
    
    ax.set_xlabel('PC{} ({}%)'.format(x_pc[m] + 1, pca_evr[x_pc[m]]))
    ax.set_ylabel('PC{} ({}%)'.format(y_pc[m] + 1, pca_evr[y_pc[m]]))

# save figure
plt.tight_layout()
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'all_strains.2d_pca.pdf'
fig.savefig(output_filename, bbox_inches='tight')
