==================================================
Scripts to analyse differentially methylated genes
==================================================

Subfolders will be described in the order that makes the most logical sense. Note that the outputs are chained--#1, #2 and #3 has to be carried out in order, while #4, #5 and #6 can be done in any order.

Brief description of subfolder contents
---------------------------------------
1. ``bias_density_medians/`` first calculates the CpG bias, methylation density and median methylation on a per-gene basis. Of these three values, in this project, I cared most about median methylation. The per-sample values (compressed as ``per_sample_bdm.tar.gz``) was parsed using ``tabulate_tsvs.py`` to produce  ``all.median_meths.tsv`` (a compressed version is in the subfolder). This file was used in the scripts in subsequent folders.

2. ``normality_testing/`` checks whether median methylation values were normally distributed. They mostly were, so GLMs were applied with a Gaussian distribution in mind, and it justifies the use of t-tests and ANOVAs.

3. ``glm_showing_independence/`` checks whether genic median methylation levels were affected by a significant interaction between developmental stage and sample origin. TL;DR: not really. Methylation levels of most genes were independently affected by developmental stage and sample origins.

4. ``fujairah_vs_abu_dhabi/`` runs a t-test to find genes that are differentially methylated due to the environment. The plotting script in this folder produces Fig. 3b.

5. ``heat_surv_vs_meth/`` correlates larval heat survival to gene methylation levels--genes with strong correlations were functionally enriched for stress- and growth-related annotations. The plotting script in this folder produces Fig. 3d.

6. ``adult_vs_sperm_vs_planula/`` runs two t-tests to identify genes that are differentially methylated in adult vs. sperm and E8/S7/S8 vs. larva. The plotting scripts in this folder produce Figs. 4b, 4c.
