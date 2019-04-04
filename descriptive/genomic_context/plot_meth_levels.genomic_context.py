#!/usr/bin/env python3

"""
> plot_meth_levels.genomic_context.py <

Uses seaborn to plot a histogram of the methylation to contrast methylation
across treatments.
"""
import csv

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def define_genomic_feature(descriptor):
    """
    From descriptor, parse and return whether position is 'intergenic',
    'exon' or 'intron'.
    """
    if 'Exon_' in descriptor:
        return 'exon'
    elif 'Intron_' in descriptor:
        return 'intron'
    else:
        return 'intergenic'

data = pd.read_table('../../merged_final_covs/all.filt.annot.merged.cov', header=None,
                     usecols=[3, 6, 10],
                     names=['meth_pct', 'gene', 'ei'])

data['temp'] = data['gene'] + data['ei']
data['temp'] = data['temp'].fillna(value='intergenic')
data['genomic_feature'] = data['temp'].apply(define_genomic_feature)

data = data.drop('temp', axis=1)
data = data.drop('gene', axis=1)
data = data.drop('ei', axis=1)

# seaborn
sns.set(style='ticks')

fig, ax = plt.subplots(figsize=(5, 2.5))
np_bins = np.linspace(0, 100, 51)
color = {'exon': '#2166ac', 'intron': '#67a9cf', 'intergenic': '#f4a582'}
for f in ['exon', 'intron', 'intergenic']:
    sns.distplot(data[data['genomic_feature'] == f]['meth_pct'],
                 bins=np_bins,
                 hist_kws={'alpha': 0.3, 'edgecolor': '#999999'},
                 kde_kws={'bw': 1},
                 color=color[f],
                 label=f.title())

plt.legend()
plt.xlim(0, 100)
plt.ylim(0, 0.07)
sns.despine(offset=10, trim=True)
ax.set_xlabel('Methylation level (%)')
ax.set_ylabel('Relative frequency')

plt.legend(loc='upper center', ncol=3, frameon=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'meth_levels.genomic_context.pdf'
fig.savefig(output_filename, bbox_inches='tight')
