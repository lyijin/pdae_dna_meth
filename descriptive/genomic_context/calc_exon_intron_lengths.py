#!/usr/bin/env python3

"""
> calc_exon_intron_lengths.py <

Parse a GFF3 file and calculate the lengths of individual exons and introns.
"""
import argparse
import statistics

import parse_gff3

parser = argparse.ArgumentParser(description="""
Parse a GFF3 file and calculate the lengths of individual exons and introns.""")

parser.add_argument('genome_gff3', metavar='gff3_file',
                    type=argparse.FileType('r'),
                    help='genome GFF3 annotation.')
args = parser.parse_args()

# read coordinates of genes and exons from .gff3 file.
scaffold_gff3 = parse_gff3.parse_gff3(args.genome_gff3, 'exon')
# as genes might contain overlapping isoforms, the longest isoform is chosen,
# if multiples exist.
scaffold_gff3 = parse_gff3.pick_longest_mRNA(scaffold_gff3)
# make sure features in all mRNAs are sorted properly (for exon numbering).
scaffold_gff3 = parse_gff3.sort_features(scaffold_gff3)

# create two dictionaries to store intron/exon lengths
exon_lengths = {}
intron_lengths = {}

for s in scaffold_gff3:
    for gene_id in scaffold_gff3[s]:
        # sole_mRNA is a mRNA object, with the assumption that *.mRNAs
        # returns a dict containing {mRNA_name: mRNA_object}
        sole_mRNA = list(scaffold_gff3[s][gene_id].mRNAs.values())[0]
        
        exon_coords = sole_mRNA.details['exon']
        intron_coords = [(x[1], y[0]) for x, y in zip(exon_coords[:-1],
                                                      exon_coords[1:])]
        
        for n, el in enumerate(exon_coords):
            # calculate exon number (e) and inverse exon number (inv_e)
            e = n + 1
            if e not in exon_lengths: exon_lengths[e] = []
            exon_lengths[e].append(abs(el[1] - el[0]))
            
            inv_e = n - len(exon_coords)
            if inv_e not in exon_lengths: exon_lengths[inv_e] = []
            exon_lengths[inv_e].append(abs(el[1] - el[0]))
        
        for n, il in enumerate(intron_coords):
            # same pattern for introns
            i = n + 1
            if i not in intron_lengths: intron_lengths[i] = []
            intron_lengths[i].append(abs(il[1] - il[0]))
            
            inv_i = n - len(intron_coords)
            if inv_i not in intron_lengths: intron_lengths[inv_i] = []
            intron_lengths[inv_i].append(abs(il[1] - il[0]))

# print statistics out
for n in range (-10, 11):
    if not n: continue
    print ('Exon {}: {} bp'.format(n, round(statistics.mean(exon_lengths[n]))))
    
for n in range (-10, 11):
    if not n: continue
    print ('Intron {}: {} bp'.format(n, round(statistics.mean(intron_lengths[n]))))
