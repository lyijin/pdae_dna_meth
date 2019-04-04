#!/usr/bin/env python3

"""
> filter_pos.five_criteria.py <

Script checks the "compiled_coverage_meth_unmeth.pre-filtered.tsv" file, and 
filters for positions that have median coverage >= 10 (criteria I) and present
in ALL samples (i.e. unmethylated + methylated reads > 0) (criteria II). To
make sure that lowly methylated loci are truly methylated, at least one of the 
treatments needs to have >= 1 methylated read in all replicates (criteria III).

The script reads in positions considered bona fide ("filter_miscalled_Cs.py" 
ran on "all.pre-filtered.merged.cov"), then filters out positions that are not
in this file. (criteria IV).

The script then checks the table of all potential SNPs (in the Bis-SNP folder)
to remove methylated positions that might be a potential SNP.

Afterwards, the script reads in meth_extract_covs/*.cov files, filters out 
unwanted positions, and saves the filtered files as ${a/cov/filt.cov} files.
"""
import csv
import glob
import gzip
import statistics
import sys

import numpy as np

def lazy_int(string):
    string = string.strip()
    return int(string) if string else 0

def calc_cols(input_list):
    return sum([1 for x in input_list if x])

# bit topsy-turvy, but start with criteria IV first!
print ('Checking positions that are considered bona fide (criteria IV)...', 
       file=sys.stderr)

# bona_fide_pos[scaf][pos]: 1 == bona fide; 0 == not bona fide
bona_fide_pos = {}
tsv_reader = csv.reader(gzip.open('all.bona_fide_meth_pos.cov.gz', 'rt'),
                        delimiter='\t')

for row in tsv_reader:
    scaf = row[0]
    pos = int(row[1])
    
    if scaf not in bona_fide_pos:
        scaf_size = int(scaf.split('size')[-1])
        # scaf_size, due to gap filling etc, tend to be MERELY SUGGESTIVE of
        # the actual length. the scaffolds are usually longer than labelled.
        # to counter this: create an array twice the expected length
        bona_fide_pos[scaf] = np.zeros(scaf_size * 2, dtype=np.int8)
    
    bona_fide_pos[scaf][pos] = 1

# then deal with criteria V!
print ('Blacklisting positions that are potentially SNPs (criteria V)...', 
       file=sys.stderr)
tsv_reader = csv.reader(gzip.open(
    '../genetic_contribution/bissnp/tabulated_genotypes.tsv.gz', 'rt'),
    delimiter='\t')

next(tsv_reader)
for row in tsv_reader:
    scaf = row[0]
    pos = int(row[1])
    
    # blacklist positions by multiplying the value of the position with -1,
    # turning bona-fide-but-potentially-SNP positions into -1
    if scaf in bona_fide_pos:
        bona_fide_pos[scaf][pos] *= -1

# criteria I, II, III: check compiled, pre-filtered, meth & unmeth reads
print ('Checking full compiled table...', file=sys.stderr)
tsv_reader = csv.reader(gzip.open(
    'compiled_coverage.pre-filt.meth_unmeth.tsv.gz', 'rt'), 
    delimiter='\t')

next(tsv_reader)
for row in tsv_reader:
    if not row: continue
    
    per_sample_cov = [lazy_int(i) + lazy_int(j)
                      for i, j in zip(row[2::2], row[3::2])]
    
    # exclude egg sample with very few reads (index #13)
    del per_sample_cov[13]
    
    # criteria II: coverage >= 5 for all treatments
    if not min(per_sample_cov) >= 5: continue
    
    # criteria I: median coverage >= 10
    if not statistics.median(per_sample_cov) >= 10: continue
    
    # criteria III: all replicates in at least one treatment with multiple
    # samples have methylation, i.e., any of
    # - Fujairah adult
    # - Fujairah sperm
    # - Abu Dhabi adult
    # - Abu Dhabi sperm
    # - Abu Dhabi larvae (E7 x S8)
    # - Abu Dhabi larvae (E8 x S7)
    per_sample_meth = [lazy_int(i) for i in row[2::2]]
    
    # exclude egg sample with very few reads (index #13)
    del per_sample_meth[13]
    
    # for this dataset, there are 6 treatments, each with 3 or 4 replicates.
    # in the per_sample_meth array, the treatments have indices of
    # 0,2,4,6; 1,3,5,7; 8,10,12,14; 9,11,13,16; 17,18,19; 20,21,22
    # following lines are hacky, but whatever works!
    one_treat_all_reps_meth = False
    if all([per_sample_meth[0], per_sample_meth[2], per_sample_meth[4],
            per_sample_meth[6]]): one_treat_all_reps_meth = True
    elif all([per_sample_meth[1], per_sample_meth[3], per_sample_meth[5],
              per_sample_meth[7]]): one_treat_all_reps_meth = True
    elif all([per_sample_meth[8], per_sample_meth[10], per_sample_meth[12],
              per_sample_meth[14]]): one_treat_all_reps_meth = True
    elif all([per_sample_meth[9], per_sample_meth[11], per_sample_meth[13],
              per_sample_meth[16]]): one_treat_all_reps_meth = True
    elif all(per_sample_meth[17:20]): one_treat_all_reps_meth = True
    elif all(per_sample_meth[21:23]): one_treat_all_reps_meth = True
    
    if not one_treat_all_reps_meth: continue
    
    scaf = row[0]
    pos = int(row[1])
    
    if scaf in bona_fide_pos:
        # adds 3 because it passes 3 additional criteria here
        bona_fide_pos[scaf][pos] += 3

# start filtering for correct positions
unfilt_files = glob.glob('../meth_extract_covs/*.cov.gz')

# exclude egg sample with very few reads ('85-Abu_Dhabi-Egg')
unfilt_files = [x for x in unfilt_files if '85-Abu_Dhabi-Egg' not in x]

for u in sorted(unfilt_files):
    print ('Processing {}...'.format(u), file=sys.stderr)
    
    tsv_reader = csv.reader(gzip.open(u, 'rt'), delimiter='\t')
    output_file = u.split('/')[-1].replace('cov.gz', 'filt.cov')
    with open(output_file, 'w') as o:
        for row in tsv_reader:
            if not row: continue
            
            scaf = row[0]
            pos = int(row[1])
            
            if scaf in bona_fide_pos:
                # pick out positions that pass criteria I-IV, and avoids the
                # -ve multiplication because of criteria V
                if bona_fide_pos[scaf][pos] == 4:
                    print ('\t'.join(row), file=o)
