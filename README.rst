================================================================
Intergenerational epigenetic inheritance in reef-building corals
================================================================

.. image:: https://zenodo.org/badge/122317576.svg
   :target: https://zenodo.org/badge/latestdoi/122317576

This has been published in *Nature Climate Change*--huzzah! Link to paper: https://www.nature.com/articles/s41558-019-0687-2

Disclaimer
----------
This repository contain scripts that carry out the key analyses for this project. As I have had a few more years of analysing methylation data, I'd like to think I'm more organised now (as compared to the other coral methylation project I was on, https://github.com/lyijin/spis_dna_meth. Experience really does make a difference!). However, there is still a slim chance that I might've missed uploading some code--alert me via email, and I'll upload them here.

To conserve disk space (and avoid GitHub clobbering me about using excessive space), I have gzipped large ``*.tsv`` (tab-separated values) files. Check out the scripts in the same folder to see how the files are produced, or used.

Methylation pipeline
--------------------
The methylation pipeline used in this project is available elsewhere at https://github.com/lyijin/working_with_dna_meth, which provides a fuller description (and theoretical considerations) of the pipeline used, and the scripts written to operate on Bismark's output.

Bismark ends by producing ``*.cov`` files, which were extracted from the corresponding deduplicated ``*.bam`` files. As the samples in this project are not genetically identical, there is an additional SNP identification step that was carried out individually on these ``*.bam`` files. This is described in further detail in the README located in ``genetic_contribution/bissnp``.

After SNP identification, instead of using the ``filter_pos.four_criteria.py`` in the ``working_with_dna_meth`` repository, use ``filter_meth_pos/filter_pos.five_criteria.py`` in this repository as it contains an additional criteria: the blacklisting of methylation positions that was identified as a potential SNP in ANY of the 23 ``*.bam`` files. This might sound super conservative, but in practice, it removed only 96k positions (out of 1.51m positions, ~6.4%), which isn't too bad.

Post-filtering, ``annotate_bismark_cov.py`` added per-position annotations for each methylated position, which I then ``gzip``ped to produce the ``*.cov.gz`` files you see at the root folder. For those who are interested in the analyses instead of the data filtering steps, these files are THE key intermediate files that get used in many other scripts to produce tables/graphs.

Brief description of folder contents
------------------------------------
To regenerate the key files used in the plotting scripts, you'll have to:

1. Decompress the ``*.cov.gz`` files in this folder.

2. ``all.filt.pct.tsv`` is produced by running the command 
   
   ``tabulate_tsvs.py *.cov -k 0 1 -c 3 -v > all.filt.pct.tsv``
   (you can wget/curl ``tabulate_tsvs.py`` from https://raw.githubusercontent.com/lyijin/common/master/tabulate_tsvs.py)
   
   then
   
   ``sed -i 's/^\t\t/scaffold\tpos\t/' all.filt.pct.tsv && sed -i 's/\.cov//g' all.filt.pct.tsv``
   
   ``head all.filt.pct.tsv`` should look like::

     scaffold        pos     A1-F    A2-F    A3-F    A4-F    A5-AD   A6-AD   A7-AD   A8-AD   E8-AD   L1-AD   L2-AD   L3-AD   L4-AD   L5-AD   L6-AD   S1-F    S2-F    S3-F    S4-F    S5-AD   S6-AD   S7-AD   S8-AD
     scaffold1|size5511861   3227    39.0244 21.4286 11.7647 5.2632  0.0     24.1379 35.2941 25.0    31.4286 20.0    16.2162 15.3846 20.5128 15.2174 32.0    52.6316 29.0323 13.7255 18.75   16.6667 28.5714 25.641  24.4444
     scaffold1|size5511861   3228    20.5882 15.3846 7.1429  0.0     14.2857 28.0    25.0    9.0909  26.3158 24.0    31.4286 18.1818 31.4815 23.3333 18.1818 37.7778 16.6667 10.4167 38.4615 3.8462  28.5714 15.1515 21.6216
     scaffold1|size5511861   3231    32.5    28.2609 11.1111 5.2632  0.0     22.5806 33.3333 33.3333 29.7297 25.0    24.3243 14.6341 20.0    16.6667 28.0    51.2821 31.4286 10.0    18.75   19.2308 35.1351 28.5714 22.2222
     etc.

3. Navigate to the folder ``merged_final_covs/``. This folder contains a shell script that merges these ``*.cov`` files via semantically meaningful ways (e.g. "all Fujairah samples", "all sperm samples", etc.) to produce another bunch of ``*.cov`` files. More info provided in the folder itself.

With that out of the way, the fun part starts. Each folder described below has another README within the folder, which describes the figures plotted from the code within the nested subfolders.

``descriptive/`` is the folder containing scripts for... descriptive plots of the methylation data.

``genetic_contribution/`` analyses the effect of genotype on epigenotype (are the changes at epigenotype level dependent on the changes in genotype?).

``diff_meth_genes/`` first justifies why we separately analysed the effect of development and environment, then goes into analysing and plotting stuff for both variables respectively.

More in-depth explanations are nestled within the subfolders.
