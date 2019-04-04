#!/usr/bin/env python3

"""
> filter_high_cov_pos.py <

Script takes in `tabulated_genotype.tsv`, and based on `tabulated_depths.tsv`,
prints out positions with high coverage (minimum coverage >= 10), so that
"0/0" calls from Bis-SNP are more likely to mean homozygous reference base,
rather than unknown due to lack of coverage.

Script assumes that both files have headers, positions are ordered in the same
manner, with the same number of lines.
"""
import csv

def try_int(x):
    """
    Overrides int('') to return 0, instead of crashing.
    """
    try:
        return int(x)
    except ValueError:
        return 0

# start parsing unfiltered table
tsv_genotypes = csv.reader(open('tabulated_genotypes.tsv'), delimiter='\t')
tsv_depths = csv.reader(open('tabulated_depths.tsv'), delimiter='\t')

# print header
header = next(tsv_genotypes)
print (*header, sep='\t')
next(tsv_depths)

for row in tsv_depths:
    # read in genotypes
    covs = [try_int(x) for x in row[2:]]
    min_cov = min(covs)
    max_cov = max(covs)
    
    if min_cov >= 10 and max_cov <= 100:
        print (*next(tsv_genotypes), sep='\t')
    else:
        next(tsv_genotypes)
