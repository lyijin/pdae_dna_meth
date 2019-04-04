#!/usr/bin/env python3

"""
> calc_bias_density_medians.py <

Multi-functional script calculates, on a per-gene basis, three values:
1. the CpG bias (CpG O/E) ratios for individual genes in the genome
2. methylation density (dense vs. sparse methylation)
3. median methylation %

Requires a genome (to calculate expected CpG fraction), a GFF3 file with gene
information (to calculate observed CpG fraction) and an annotated cov file
(to extract # of methylated CpGs in each gene).
"""
import argparse
import collections
import csv
import itertools
import statistics

import natural_sort
import parse_fasta
import parse_gff3

def reverse_complement(seq):
    seq = seq.replace('U', 'T')

    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)

    seq = seq[::-1].translate(translation_table)

    return seq
    
def get_gene_sequence(sequence, coords):
    startpos, endpos = coords
    
    revcomp_seq = False
    if startpos > endpos:
        # set reverse complement flag!
        startpos, endpos = endpos, startpos
        revcomp_seq = True

    # slice sequence
    sliced_seq = sequence[startpos:endpos]
    if revcomp_seq: sliced_seq = reverse_complement(sliced_seq)
    
    return sliced_seq

parser = argparse.ArgumentParser(description="""
Multi-functional script calculates, on a per-gene basis, three values:
1. the CpG bias (CpG O/E) ratios for individual genes in the genome
2. methylation density (dense vs. sparse methylation)
3. median methylation %.""")

parser.add_argument('genome', metavar='FASTA file',
                    type=argparse.FileType('r'),
                    help='FASTA file of genome.')
parser.add_argument('gff', metavar='GFF3 file',
                    type=argparse.FileType('r'),
                    help='GFF3 file of genome.')
parser.add_argument('cov', metavar='annotated cov file',
                    type=argparse.FileType('r'),
                    help='cov file with methylated positions.')

args = parser.parse_args()

# parse genome for genomic sequences
genome_sequence = parse_fasta.get_all_sequences(args.genome, 'fasta')

# parse gff3 file
scaffold_gff3 = parse_gff3.parse_gff3(args.gff, 'gene')

# parse annotated cov file
cov_data = {}
tsv_reader = csv.reader(args.cov, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    
    gene_id = row[6]
    meth_pct = float(row[3])
    
    if gene_id not in cov_data:
        cov_data[gene_id] = []
    
    cov_data[gene_id].append(meth_pct)

# print header line
print ('gene', 'meth pos', 'observed CpG', 'expected CpG', 'CpG bias',
       'meth density', 'median meth %', sep='\t')

# note on CpG(expected) calculation:
#   E = p(C) . p(G) is vulnerable to ambiguous dinucleotides ('CD', 'CN', ...).
#   instead, calculate p(C) = p(CA) + p(CC) + p(CG) + p(CT).

# perform per-gene calculations!
valid_dinuc =  [''.join(x) for x in itertools.product('ACGT', repeat=2)]
cx_dinuc = ['C' + x for x in 'ACGT']
gx_dinuc = ['G' + x for x in 'ACGT']

results = {}
for scaf in scaffold_gff3:
    for gene in scaffold_gff3[scaf]:
        gene_sequence = get_gene_sequence(genome_sequence[scaf],
                                          scaffold_gff3[scaf][gene].coords)
        gene_dinuc = collections.Counter([''.join(x) for x in \
            zip(gene_sequence[:-1], gene_sequence[1:])])
        cpg_observed = gene_dinuc['CG']
        cpg_expected = sum([gene_dinuc[x] for x in cx_dinuc]) *\
                       sum([gene_dinuc[x] for x in gx_dinuc])
        cpg_expected = round(cpg_expected / len(gene_sequence), 2)
        cpg_bias = round(cpg_observed / cpg_expected, 4)
        
        # believe it or not, there are genes WITHOUT a single CpG dinuc...
        if not cpg_observed: continue
        
        # also print out stats for non-methylated genes, there'll just be
        # lots of 0s in its results
        if gene not in cov_data:
            results[gene] = [0, cpg_observed, cpg_expected, cpg_bias, 0, 0]
        else:
            # cpg_observed counts dinucleotides, while cov_data[gene] stores 
            # data on a per-position basis, i.e. lengths can differ up to 2x.
            meth_density = round(len(cov_data[gene]) / 2 / cpg_observed * 100, 2)
            median_meth_pct = round(statistics.median(cov_data[gene]), 2)
            results[gene] = [len(cov_data[gene]), cpg_observed, 
                             cpg_expected, cpg_bias,
                             meth_density, median_meth_pct]

# output
for r in natural_sort.natural_sort(results):
    # replace 'PdaeGene' with 'Pdae'
    print (r.replace('PdaeGene', 'Pdae'), *results[r], sep='\t')
