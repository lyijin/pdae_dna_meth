==============================================================
Scripts to analyse genetic contribution to epigenetic patterns
==============================================================

Only two folders here, but these are... not so straightforward. Scripts in these subfolders have to be carried out in order. #2 depends on the results of #1.

Brief description of subfolder contents
---------------------------------------
1. ``bissnp/`` identified SNPs from WGBS data (from Bismark). The folder contains commands used to run Bis-SNP, and key intermediate files parsed from the raw outputs. This folder contains another README to walk the reader through the pipeline, and how I used the SNP data for filtering/plotting purposes. Plotting script in this folder produces Fig. 2a.

2. ``calc_indiv_fst/`` calculates F\ :sub:`ST` values on a per-gene basis, then correlates relative genotypic differences with relative epigenotypic differences. The plotting script in this folder produces Fig. 2b and Supplementary Fig. S3.
