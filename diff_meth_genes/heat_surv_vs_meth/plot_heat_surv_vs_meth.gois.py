#!/usr/bin/env python3

"""
> plot_heat_surv_vs_meth.gois.py <

For the genes of interest, plot a lmplot to illustrate the correlation of
survival data with methylation levels.
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# hardcoded survival index values
surv_index = {'S1-F': 0.789439953,
              'S3-F': 0.998209258,
              'S4-F': 0.22654263,
              'S5-AD': 0.954202839,
              'S6-AD': 1.477754999,
              'S8-AD': 1.300719112}

# read from output file from previous script
data = pd.read_table('heat_surv_vs_meth.tsv')

# retain genes of interest
gois = ['Pdae491', 'Pdae479', 'Pdae2659', 'Pdae16130']
data = data[data['gene'].isin(gois)]
data = data.loc[:, 'gene':'S8-AD']

# melt data for plotting, then add survival data and sample location
data = data.melt(id_vars='gene', var_name='sample', value_name='meth_pct')
data['surv_index'] = data['sample'].apply(lambda x: surv_index[x])
data['location'] = data['sample'].apply(lambda x: x.split('-')[-1])

# plot the lmplot!
sns.set(style='ticks', font_scale=1)
fig, ax = plt.subplots(figsize=(7, 7))

sns.lmplot(x='meth_pct', y='surv_index', data=data, col='gene',
           col_wrap=2, height=3,
           sharex=False, sharey=True,
           ci=None,
           scatter_kws={'color': ['#6baed6'] * 3 + ['#fcae91'] * 3},
           line_kws={'color': '#000000', 'alpha': 0.7})

sns.despine(offset=5, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'heat_surv_vs_meth.gois.pdf'
fig.savefig(output_filename, bbox_inches='tight')
