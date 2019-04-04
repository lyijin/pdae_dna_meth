#!/bin/bash
# 
# > merge_covs.sh <
# 
# Merge bismark covs into semantically meaningful files.

# merges relevant files
merge_bismark_cov.py ../A*.cov > all-Adult.filt.annot.merged.cov
merge_bismark_cov.py ../S*.cov > all-Sperm.filt.annot.merged.cov
merge_bismark_cov.py ../L*.cov > all-Larval.filt.annot.merged.cov

merge_bismark_cov.py ../A?-F.cov ../S?-F.cov > all-Fujairah.filt.annot.merged.cov
merge_bismark_cov.py ../A?-AD.cov ../S?-AD.cov > all-Abu_Dhabi.filt.annot.merged.cov

merge_bismark_cov.py ../*.cov > all.filt.annot.merged.cov

# port annotations over
cut -f 7- ../A1-F.cov > tmp
for a in *.cov; do paste $a tmp > tmp2 && mv -f tmp2 $a; done
rm -f tmp