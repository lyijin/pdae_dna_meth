#!/usr/bin/env python3

"""
> plot_ridgeline.genomic_context.loc.py <

Uses seaborn to plot a ridgeline plot of the methylation to contrast methylation
across locations.

Borrows code from https://seaborn.pydata.org/examples/kde_ridgeplot.html.
"""
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

# define and use a simple function to label the plot in axes coordinates
def label(x, color, label):
    ax = plt.gca()
    ax.text(.025, .025, label, fontweight='bold', color=color,
            ha='right', va='center', transform=ax.transAxes)

combined_data = pd.DataFrame({'meth_pct': [], 'genomic_feature': [],
                              'loc': []})
for n, loc in enumerate(['Fujairah', 'Abu_Dhabi']):
    data = pd.read_table(
        '../../merged_final_covs/all-{}.filt.annot.merged.cov'.format(loc),
        header=None,
        usecols=[3, 6, 10],
        names=['meth_pct', 'gene', 'ei'])
    
    data['temp'] = data['gene'] + data['ei']
    data['temp'] = data['temp'].fillna(value='intergenic')
    data['genomic_feature'] = data['temp'].apply(define_genomic_feature)
    
    data = data.drop('temp', axis=1)
    data = data.drop('gene', axis=1)
    data = data.drop('ei', axis=1)
    
    data['loc'] = loc
    
    combined_data = combined_data.append(data)

# seaborn stuff (requires seaborn > 0.9.0)
sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

# initialize the FacetGrid object
g = sns.FacetGrid(combined_data, row='loc', col='genomic_feature', hue='loc',
                  height=1.2, aspect=3, palette=['#67a9cf', '#ef8a62'])

# draw the densities in a few steps
g.map(sns.kdeplot, 'meth_pct', clip_on=False, shade=True, alpha=1, lw=1.5, bw=4)
g.map(sns.kdeplot, 'meth_pct', clip_on=False, color="w", lw=2, bw=4)
g.map(plt.axhline, y=0, lw=2, clip_on=False)

# label the lines with the `label` def
g.map(label, 'meth_pct')

# set the subplots to overlap
g.fig.subplots_adjust(wspace=0.1, hspace=-.55)

# remove axes details that don't play well with overlap
g.set(xlim=(-15, 115))
g.set_titles('{col_name}')
g.set(yticks=[])
g.despine(bottom=True, left=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'ridgeline.genomic_context.loc.pdf'
fig.savefig(output_filename, bbox_inches='tight')
