#!/usr/bin/env python3

"""
> calc_hudson_fst.py <

Based on the output of BisSNP analysis, calculate the numerators and
denominators of Hudson's Fst (as formulated by Bhatia et. al, 2013) 
for every consistent position.

Equation 10 from the publication was used to calculate the numerators and 
denominators. The reason why both are calculated separately is because the same
publication recommends the calculation of Fst in a window containing multiple 
SNPs as the ratio of averages (i.e. sum N / sum D) instead of average of ratios
(i.e. mean (N/D)).

Hudson Fst values are saved directly into `hudson_fsts.tsv`.
"""
import csv

def calc_numerator(p1, p2, n1, n2):
    """
    Refer to equation 10 of the publication!
    
    If negative values are calculated, convert that to 0 as negative Fst 
    makes no sense.
    """
    a = (p1 - p2) ** 2
    b = (p1 * (1 - p1)) / (n1 - 1)
    c = (p2 * (1 - p2)) / (n2 - 1)
    return max(a - b - c, 0)

def calc_denominator(p1, p2):
    """
    Refer to equation 10 of the publication!
    """
    return p1 * (1 - p2) + p2 * (1 - p1)

# start parsing unfiltered table
tsv_reader = csv.reader(open('../bissnp/tabulated_genotypes.cov10-100.AS.annot.tsv'), 
                        delimiter='\t')

# print header
with open('hudson_fsts.tsv', 'w') as output_file:
    old_header = next(tsv_reader)
    # 'ei' stores exon/intron information
    print (*old_header[:2], 'gene', 'ei', 'N', 'D', sep='\t', file=output_file)
    
    for row in tsv_reader:
        if not row: continue
        
        f_genotypes = [int(x) for x in row[2:6] + row[17:21]]
        ad_genotypes = [int(x) for x in row[6:10] + row[21:25]]
        
        # rename PdaeGeneX to PdaeX
        gene = row[25].replace('Gene', '')
        
        # if genic, check whether exon or intron
        if 'Pdae' in row[25]:
            if 'Exon_' in row[29]:
                ei = 'exon'
            elif 'Intron_' in row[29]:
                ei = 'intron'
        else:
            ei = 'intergenic'
        
        # refresher for genotype: 0 = homozygous reference base
        #                         1 = heterozygous
        #                         2 = homozygous non-reference base
        # hence the proportion of non-ref base = sum(genotype) / (2 * n)
        n1 = len(f_genotypes)
        n2 = len(ad_genotypes)
        p1 = sum(f_genotypes) / (2 * n1)
        p2 = sum(ad_genotypes) / (2 * n2)
        
        N = calc_numerator(p1, p2, n1, n2)
        D = calc_denominator(p1, p2)
        
        print (*row[:2], gene, ei, round(N, 4), round(D, 4), sep='\t',
               file=output_file)
