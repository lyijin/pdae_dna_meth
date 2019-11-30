===================================================
Scripts to calculate F\ :sub:`ST` values on a per-gene basis
===================================================

Script order is fairly straightforward here:

1. ``calc_hudson_fst.py`` implements the calculation of F\ :sub:`ST` values. It produces an intermediate file called ``hudson_fsts.tsv`` which is ~200 MB, hence it's not included in this repository. This intermediate file is used only in the next script, which is...

2. ... ``tabulate_data_per_gene.py``, which calculates per-gene F\ :sub:`ST` values and tabulates them with delta methylation values. The output file ``delta_meth_vs_fst_gene.tsv`` is provided in compressed form here. This file is used in all three plotting scripts.

Plotting scripts:

1. ``plot_pairplot.fst_vs_meth.py``: original plot, tries to make the case that genetic changes and epigenetic changes are weakly correlated (Fig 2b).

2. Supervisor: what if different populations had different trends? Me: ``plot_pairplot.fst_vs_meth.loc_split.py`` (Supplementary Fig. S3).

3. Reviewer: what if F\ :sub:`ST` isn't the best way to gauge genetic difference--how about allelic difference? Me: ``plot_pairplot.allelic_diff_vs_meth.py`` (unpublished, but shows that the weak positive correlation is still there (and r still ~0.1)).
