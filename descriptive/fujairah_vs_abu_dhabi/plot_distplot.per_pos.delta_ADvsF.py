#!/usr/bin/env python3

"""
> plot_distplot.per_pos.delta_ADvsF.py <

Uses seaborn to plot a histogram of the methylation differences between
Fujairah and Abu Dhabi samples.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# read data
data = pd.read_table('../../all.filt.pct.tsv')
data = data.filter(regex='^A|^S')
data['f_mean'] = data.filter(regex='-F').mean(axis=1)
data['ad_mean'] = data.filter(regex='-AD').mean(axis=1)

# calculate delta only if both methylation levels are > 0
data['delta_advsf'] = data.apply(lambda x: x['ad_mean'] - x['f_mean']
                                     if all([x['ad_mean'], x['f_mean']])
                                     else np.nan, axis=1)

data = data['delta_advsf'].dropna()

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
ax.set_xlabel('Methylation level Abu Dhabi - Fujairah')
ax.set_ylabel('Relative frequency')

sns.despine(offset=5, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'distplot.per_pos.delta_ADvsF.pdf'
fig.savefig(output_filename, bbox_inches='tight')
