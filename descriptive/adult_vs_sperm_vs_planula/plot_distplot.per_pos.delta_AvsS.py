#!/usr/bin/env python3

"""
> plot_distplot.per_pos.delta_AvsS.py <

Uses seaborn to plot a histogram of the methylation differences between adult
and sperm samples.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# read data
data = pd.read_table('../../all.filt.pct.tsv')
data['a_mean'] = data.filter(regex='^A').mean(axis=1)
data['s_mean'] = data.filter(regex='^S').mean(axis=1)

# calculate delta only if both methylation levels are > 0
data['delta_svsa'] = data.apply(lambda x: x['s_mean'] - x['a_mean']
                                    if all([x['s_mean'], x['a_mean']])
                                    else np.nan, axis=1)

data = data['delta_svsa'].dropna()

# seaborn stuff!
sns.set_style('white')
sns.set_style('ticks')

fig, ax = plt.subplots(figsize=(4, 3))
np_bins = np.linspace(-20, 20, 41)
sns.distplot(data,
             bins=np_bins,
             hist_kws={'alpha': 0.3, 'edgecolor': '#999999'},
             kde_kws={'bw': 1})

# write down mean/median values somewhere
ax.text(np.median(data), 0.075, f'Median: {np.median(data):.1f}%')
ax.text(np.mean(data), 0.07, f'Mean: {np.mean(data):.1f}%')

ax.set_xlim([-20.20, 20.20])
ax.set_xlabel('Methylation level Sperm - Adult')
ax.set_ylabel('Relative frequency')

sns.despine(offset=5, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'distplot.per_pos.delta_AvsS.pdf'
fig.savefig(output_filename, bbox_inches='tight')
