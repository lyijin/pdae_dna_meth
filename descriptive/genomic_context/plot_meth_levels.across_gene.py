#!/usr/bin/env python3

"""
> plot_meth_levels.across_gene.py <

Perhaps plotting line plots are clearer than heatmaps...

Plots methylation %, from left to right of the plot:
1. 4 kb upstream
2. exon 1
3. intron 1
4. exon 2
5. intron 2
6. exon 3
7. nothing
8. exon -3
9. intron -2
10. exon -2
11. intron -1
12. exon -1
13. 4kb downstream

NOTE: one needs to consider that different introns/exons have different lengths,
which will bias the density of methylated CpGs in those features. This is why
each genomic feature had to have its mean length calculated (with a different
script) prior to the plotting of this graph.
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

EXON_LENGTH = {1: 286, 2: 320, 3: 225, -3: 203, -2: 270, -1: 380}
INTRON_LENGTH = {1: 1971, 2: 1744, -2: 1598, -1: 1728}
START = -4000
NOTHING_START = EXON_LENGTH[1] + INTRON_LENGTH[1] + EXON_LENGTH[2] + \
                INTRON_LENGTH[2] + EXON_LENGTH[3]
NOTHING_END = NOTHING_START + 1000
END = NOTHING_END + EXON_LENGTH[-3] + INTRON_LENGTH[-2] + EXON_LENGTH[-2] + \
      INTRON_LENGTH[-1] + EXON_LENGTH[-1] + 4000

def get_starting_loc(feature):
    """
    Hardcode starting coordinates for each of the 13 starting regions listed
    above.
    """
    if feature == 'upstream':
        return START
    elif feature == 'exon 1':
        return 0
    elif feature == 'intron 1':
        return EXON_LENGTH[1]
    elif feature == 'exon 2':
        return EXON_LENGTH[1] + INTRON_LENGTH[1]
    elif feature == 'intron 2':
        return EXON_LENGTH[1] + INTRON_LENGTH[1] + EXON_LENGTH[2]
    elif feature == 'exon 3':
        return EXON_LENGTH[1] + INTRON_LENGTH[1] + EXON_LENGTH[2] + INTRON_LENGTH[2]
    elif feature == 'nothing':
        return NOTHING_START
    elif feature == 'exon -3':
        return NOTHING_END
    elif feature == 'intron -2':
        return NOTHING_END + EXON_LENGTH[-3]
    elif feature == 'exon -2':
        return NOTHING_END + EXON_LENGTH[-3] + INTRON_LENGTH[-2]
    elif feature == 'intron -1':
        return NOTHING_END + EXON_LENGTH[-3] + INTRON_LENGTH[-2] + EXON_LENGTH[-2]
    elif feature == 'exon -1':
        return NOTHING_END + EXON_LENGTH[-3] + INTRON_LENGTH[-2] + \
               EXON_LENGTH[-2] + INTRON_LENGTH[-1]
    elif feature == 'downstream':
        return NOTHING_END + EXON_LENGTH[-3] + INTRON_LENGTH[-2] + \
               EXON_LENGTH[-2] + INTRON_LENGTH[-1] + EXON_LENGTH[-1]
    else:
        raise ValueError('feature has an unexpected value!')

# read the cov file
data = pd.read_table('../../merged_final_covs/all.filt.annot.merged.cov', header=None,
                     names=['scaf', 'start_pos', 'end_pos', 'meth_pct', 
                            'meth', 'unmeth', 'gene', 'gene_relpos',
                            'gene_5p', 'gene_3p', 'ei', 'ei_relpos', 
                            'ei_5p', 'ei_3p'])

# remove lines with 'no_info'
data = data[data['gene'] != 'no_info'][data['ei'] != 'no_info']

# inefficient code but... sod it. manually filter cov files for desired contexts.
# grab exonic contexts: exons 1-3, then exons (-3)-(-1)
exons = []
for n in range(1, 4):
    tmp = data[data['ei'].str.match('Exon_{}_'.format(n))]
    float_relpos = pd.to_numeric(tmp['ei_relpos'])
    
    tmp2 = pd.DataFrame()
    tmp2['meth_pct'] = tmp['meth_pct']
    tmp2['recalibrated_x'] = get_starting_loc('exon {}'.format(n)) + \
                             float_relpos * EXON_LENGTH[n]
    exons.append(tmp2)

for n in range(-3, 0):
    tmp = data[data['ei'].str.match('Exon_\d+_{}$'.format(n))]
    float_relpos = pd.to_numeric(tmp['ei_relpos'])
    
    tmp2 = pd.DataFrame()
    tmp2['meth_pct'] = tmp['meth_pct']
    tmp2['recalibrated_x'] = get_starting_loc('exon {}'.format(n)) + \
                             float_relpos * EXON_LENGTH[n]
    exons.append(tmp2)

exons = pd.concat(exons)

# grab intronic contexts: introns 1-2, then introns (-2), (-1)
introns = []
for n in range(1, 3):
    tmp = data[data['ei'].str.match('Intron_{}_-\d+$'.format(n))]
    float_relpos = pd.to_numeric(tmp['ei_relpos'])
    
    tmp2 = pd.DataFrame()
    tmp2['meth_pct'] = tmp['meth_pct']
    tmp2['recalibrated_x'] = get_starting_loc('intron {}'.format(n)) + \
                             float_relpos * INTRON_LENGTH[n]
    introns.append(tmp2)

for n in range(-2, 0):
    tmp = data[data['ei'].str.match('Intron_\d+_{}$'.format(n))]
    float_relpos = pd.to_numeric(tmp['ei_relpos'])
    
    tmp2 = pd.DataFrame()
    tmp2['meth_pct'] = tmp['meth_pct']
    tmp2['recalibrated_x'] = get_starting_loc('intron {}'.format(n)) + \
                             float_relpos * INTRON_LENGTH[n]
    introns.append(tmp2)

introns = pd.concat(introns)

# grab upstream region
upstreams = []
tmp = data[data['ei'].str.match('upstream_crick\|')]
tmp2 = pd.DataFrame()
tmp2['meth_pct'] = tmp['meth_pct']
tmp2['recalibrated_x'] = -tmp['gene_5p'] - 1
upstreams.append(tmp2)

tmp = data[data['ei_relpos'].str.match('downstream_watson\|')]
tmp2 = pd.DataFrame()
tmp2['meth_pct'] = tmp['meth_pct']
tmp2['recalibrated_x'] = tmp['gene_3p']
upstreams.append(tmp2)

upstreams = pd.concat(upstreams)
upstreams = upstreams[upstreams['recalibrated_x'] >= START]

# grab downstream region
downstreams = []
tmp = data[data['ei'].str.match('upstream_watson\|')]
tmp2 = pd.DataFrame()
tmp2['meth_pct'] = tmp['meth_pct']
tmp2['recalibrated_x'] = get_starting_loc('downstream') + tmp['gene_5p'] + 1
downstreams.append(tmp2)

tmp = data[data['ei_relpos'].str.match('downstream_crick\|')]
tmp2 = pd.DataFrame()
tmp2['meth_pct'] = tmp['meth_pct']
tmp2['recalibrated_x'] = get_starting_loc('downstream') - tmp['gene_3p']
downstreams.append(tmp2)

downstreams = pd.concat(downstreams)
downstreams = downstreams[downstreams['recalibrated_x'] <= END]

# combine everything together
combined = pd.concat([upstreams, exons, introns, downstreams])

# plot the heatmap!
sns.set_style('white')
sns.set_style('ticks')

fig, ax = plt.subplots(figsize=(9, 3))
sns.distplot(combined['recalibrated_x'], bins=range(START, END + 50, 50),
             kde_kws={'color': '#e34a33', 'bw': 100}, 
             hist_kws={'color': '#fdbb84', 'edgecolor': '#999999'})

plt.xlabel('Genomic context')
plt.ylabel('Relative frequency')

plt.xlim(START, END)
#plt.ylim(0, 100)
xtick_labels = ['4 kb\nupstream', 'exon 1', 'intron 1', 'exon 2', 'intron 2',
                'exon 3', '', 'exon -3', 'intron -2', 'exon -2', 'intron -1',
                'exon -1', '4 kb\ndownstream']
xt = [-2000,
      get_starting_loc('exon 1') + 0.5 * EXON_LENGTH[1],
      get_starting_loc('intron 1') + 0.5 * INTRON_LENGTH[1],
      get_starting_loc('exon 2') + 0.5 * EXON_LENGTH[2],
      get_starting_loc('intron 2') + 0.5 * INTRON_LENGTH[2],
      get_starting_loc('exon 3') + 0.5 * EXON_LENGTH[3],
      NOTHING_START + 500,
      get_starting_loc('exon -3') + 0.5 * EXON_LENGTH[-3],
      get_starting_loc('intron -2') + 0.5 * INTRON_LENGTH[-2],
      get_starting_loc('exon -2') + 0.5 * EXON_LENGTH[-2],
      get_starting_loc('intron -1') + 0.5 * INTRON_LENGTH[-1],
      get_starting_loc('exon -1') + 0.5 * EXON_LENGTH[-1],
      END - 2000]
plt.xticks(xt, xtick_labels, rotation=45, ha='right')

# divide plot using solid/broken lines
plt.axvline(x=0, lw=1, color='k')
plt.axvline(x=get_starting_loc('downstream'), lw=1, color='k')
plt.axvline(x=get_starting_loc('nothing'), lw=1, ls='--', color='k')
for n in range(-3, 4):
    if not n: continue
    plt.axvline(x=get_starting_loc('exon {}'.format(n)),
                lw=1, ls='--', color='k')

for n in range(-2, 3):
    if not n: continue
    plt.axvline(x=get_starting_loc('intron {}'.format(n)),
                lw=1, ls='--', color='k')

sns.despine(offset=10, bottom=True, trim=True)

#fig.tight_layout()

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'meth_levels.across_gene.pdf'
fig.savefig(output_filename, bbox_inches='tight')
