# find 100% matches to Cirelli genes in the LJI contigs

import glob
import csv
import os.path
from Bio import SeqIO

infiles = 'matches/lji_7g_matches_cirelli.csv'
outfile = 'matches/lij_best_cirelli_matches.csv'

cirelli_ref = 'reference_sets/cirelli_IGH.fasta'
cirelli_germlines = SeqIO.parse(cirelli_ref, 'fasta')
cirelli_rh10 = {}


def gene_number(gene):
    num = gene.split('.')[-1].zfill(3) if gene.split('.')[-1] != 'a' else gene.split('.')[-2].zfill(3) + 'a'
    type = gene.split('_')[1][:4]
    return type + num


for ref in cirelli_germlines:
    cirelli_rh10[ref.description] = str(ref.seq).upper()

best = {}

for k in cirelli_rh10.keys():
    best[k] = {'contig': 'Not found', 'gene_number': gene_number(k), 'cirelli_match': k, 'cirelli_score': '0'}

for k1, v1 in cirelli_rh10.items():
    for k2, v2 in cirelli_rh10.items():
        if k1 != k2 and v1 == v2:
            best[k1]['contig'] = 'Identical to ' + k2
            best[k2]['contig'] = 'Identical to ' + k1

with open(outfile, 'w', newline='') as fo:
    fieldnames = ['contig', 'sense', 'gene_number', 'start', 'end', 'start_rev', 'end_rev', 'evalue', 'matches', 'best', 'imgt_match', 'imgt_score', 'cirelli_match', 'cirelli_score', 'l_part1', 'l_part2', 'v_heptamer', 'v_nonamer', 'j_heptamer', 'j_nonamer', 'j_frame', 'd_3_heptamer', 'd_3_nonamer', 'd_5_heptamer', 'd_5_nonamer', 'functional', 'notes', 'aa', 'v-gene_aligned_aa', 'seq']
    writer = csv.DictWriter(fo, fieldnames=fieldnames, restval='')
    writer.writeheader()

    for fp in glob.glob(infiles):
        with open(fp, 'r') as fi:
            reader = csv.DictReader(fi, dialect='excel')
            contig = os.path.basename(fp).replace('.csv', '').split('_')[-1]
            for row in reader:
                if row['best'] in best:
                    if float(row['cirelli_score']) > float(best[row['best']]['cirelli_score']):
                        row['contig'] = contig
                        row['gene_number'] = best[row['best']]['gene_number']
                        best[row['best']] = row

    for gene in best:
        writer.writerow(best[gene])

