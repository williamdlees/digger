# Functions used by dig commands
import csv
import os
from importlib.resources import files

from Bio import Align, Entrez, SeqIO
from receptor_utils import simple_bio_seq as simple
from receptor_utils.novel_allele_name import name_novel

try:
    from motif import Motif
    from search_motifs import find_compound_motif, MotifResult, SingleMotifResult
    from parse_genes import process_v, process_d, process_j, process_c, VAnnotation
except:
    from digger.motif import Motif
    from digger.search_motifs import find_compound_motif, MotifResult, SingleMotifResult
    from digger.parse_genes import process_v, process_d, process_j, process_c, VAnnotation


def find_target_sequence(assembly, target_seq):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'  # Set alignment mode to global
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    aligner.match_score = 2
    aligner.mismatch_score = 0

    alignments = aligner.align(assembly, target_seq)  # Perform pairwise alignment

    # Extract the best alignment from the alignments list
    best_alignment = alignments[0]

    return best_alignment.aligned[0][0][0], alignments[0].aligned[0][-1][-1], best_alignment.score

# Fetch a sequence from genbank given the accession number
def get_genbank_sequence(acc):
    Entrez.email = 'william@lees.org.uk'
    handle = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')

    if ',' in acc:
        seq_records = SeqIO.parse(handle, "fasta")
        seqs = {}
        for seq_record in seq_records:
            seqs[seq_record.id.split('.')[0]] = (str(seq_record.seq))
        return seqs
    else:
        seq_record = SeqIO.read(handle, "fasta")
        return str(seq_record.seq)


def process_sequence(assembly, genbank_acc, target, germlines, v_gapped_ref, v_ungapped_ref, motifs, conserved_motif_seqs, motif_params):
    gene_type = target[3]
    target_seq = germlines[target]
    sense = '+'
    notes = []

    start, end, score = find_target_sequence(assembly, target_seq)
    rev_start, rev_end, rev_score = find_target_sequence(simple.reverse_complement(assembly), target_seq)

    if rev_score > score:
        start, end, score = rev_start, rev_end, rev_score
        assembly = simple.reverse_complement(assembly)
        sense = '-'

    score = (100*score)/(2*len(target_seq))
    score = round(score, 1)

    if score < 80:
        notes.append('low alignment score')

    if gene_type not in ['V', 'J', 'D', 'C']:
        notes.append('unknown gene type')

    if notes:
        row = {
            'target_allele': target,
            'genbank_acc': genbank_acc,
            'genbank_seq': assembly,
            'alignment_score': score,
            'nt_diff': 'NA',
            'start': 'NA',
            'end': 'NA',
            'end_rev': 'NA',
            'start_rev': 'NA',
            'sense': sense,
            'gene_type': gene_type,
            'notes': ','.join(notes),
            'snps': 'NA',
        }
        return [row]

    if len(target_seq) > (end - start):
        notes.append('Aligned sequence is truncated compared to reference')

    if start == 0:
        notes.append("Sequence to be annotated does not contain 5' regulatory region")

    if end == len(assembly):
        notes.append("Sequence to be annotated does not contain 3' regulatory region")

    # handle 'naked' V-regions separately to avoid spurious RSS hunting
    if gene_type == 'V' and start == 0 and end == len(assembly):
        target_gapped_ref = {target: v_gapped_ref[target]}
        target_ungapped_ref = {target: v_ungapped_ref[target]}
        v_annot = VAnnotation(assembly, 1, len(assembly), refs=[target_gapped_ref, target_ungapped_ref])
        v_annot.annotate()

        nt_diff = simple.nt_diff(assembly, target_seq)

        snps = ''
        if nt_diff > 0:
            snps = name_novel(v_annot.gapped, target_gapped_ref, True)[0]
            snps = snps.replace(target, '')

        row = {
            'target_allele': target,
            'snps': snps,
            'genbank_acc': genbank_acc,
            'genbank_seq': assembly,
            'alignment_score': score,
            'nt_diff': nt_diff,
            'start': 0,
            'end': end,
            'end_rev': rev_start,
            'start_rev': rev_end,
            'sense': sense,
            'gene_type': gene_type,
            'notes': ', '.join(notes),
            'aa': v_annot.gapped_aa.replace('.', ''),
            'v-gene_aligned_aa': v_annot.gapped_aa,
            'seq': v_annot.ungapped,
            'seq_gapped': v_annot.gapped,

        }
        return [row]


    nt_diff = simple.nt_diff(assembly[start:end], target_seq)
    best = {}
    best['subject'] = target
    best['evalue'] = 100
    matches = []
    align = True
    v_parsing_errors = {}
    assembly_rc = simple.reverse_complement(assembly)

    if gene_type == 'V':
        if v_gapped_ref is None:
            print('Error - please specify a gapped reference set for V gene analysis')
            exit(1)
        rows = process_v(assembly, assembly_rc, germlines, v_gapped_ref, v_ungapped_ref, conserved_motif_seqs, motifs, start, end, best, matches, align,
                         motif_params['V_RSS_SPACING'], v_parsing_errors)
    elif gene_type == 'J':
        rows = process_j(assembly, assembly_rc, germlines, conserved_motif_seqs, motifs, start, end, best, matches,
                         motif_params['J_TRP_MOTIF'], motif_params['J_TRP_OFFSET'], motif_params['J_SPLICE'])
    elif gene_type == 'D':
        rows = process_d(assembly, assembly_rc, germlines, conserved_motif_seqs, motifs, start, end, best, matches)
    elif gene_type == 'C':
        rows = process_c(assembly, assembly_rc, germlines, start, end, best, matches)

    for row in rows:
        row['target_allele'] = target
        row['genbank_acc'] = genbank_acc
        row['genbank_seq'] = assembly
        row['alignment_score'] = score
        row['sense'] = sense
        row['nt_diff'] = nt_diff

        snps = ''
        if nt_diff > 0:
            if gene_type == 'V':
                target_gapped_ref = {target: v_gapped_ref[target]}
            else:
                target_gapped_ref = {target: germlines[target]}
            if 'seq_gapped' in row and row['seq_gapped']:
                snps = name_novel(row['seq_gapped'], target_gapped_ref, True)[0]
            else:
                snps = name_novel(row['seq'], target_gapped_ref, True)[0]
            snps = snps.replace(target, '')

        row['snps'] = snps

        if notes:
            row['notes'] = ', '.join(notes) + ',' + row['notes']

    return rows


def print_result(row):
    for field in ['target_allele', 'genbank_acc', 'genbank_seq', 'alignment_score', 'nt_diff', 'snps', 'start', 'end', 'sense']:
        if field == 'genbank_seq' and len(row['genbank_seq']) > 750:
            row['genbank_seq'] = f'length: {len(row["genbank_seq"])}nt'
        print(f'{field}: {row[field]}')

    for field in ['l_part1', 'l_part2', 'v_heptamer', 'v_nonamer', 'j_heptamer', 'j_nonamer', 'j_frame', 'd_3_heptamer', 'd_3_nonamer', 'd_5_heptamer', 'd_5_nonamer', 'notes']:
        if field in row and row[field]:
            print(f'{field}: {row[field]}')

    print('\n\n')


def process_output(args, rows):
    fieldnames = ['target_allele', 'genbank_acc', 'alignment_score', 'nt_diff', 'snps', 'start', 'end', 'start_rev', 'end_rev', 'sense', 'gene_type']
    fieldnames.extend(
        ['functional', 'notes', 'l_part1', 'l_part2', 'v_heptamer', 'v_nonamer', 'j_heptamer', 'j_nonamer', 'j_frame', 'd_3_heptamer', 'd_3_nonamer',
         'd_5_heptamer', 'd_5_nonamer',
         'aa', 'v-gene_aligned_aa', 'seq', 'seq_gapped', '5_rss_start', '5_rss_start_rev', '5_rss_end', '5_rss_end_rev',
         '3_rss_start', '3_rss_start_rev', '3_rss_end', '3_rss_end_rev', 'l_part1_start', 'l_part1_start_rev', 'l_part1_end', 'l_part1_end_rev',
         'l_part2_start', 'l_part2_start_rev', 'l_part2_end', 'l_part2_end_rev'])
    if args.out_file:
        with open(args.out_file, 'w', newline='\n') as fo:
            writer = csv.DictWriter(fo, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
    else:
        for row in rows:
            print_result(row)


def read_vs(args):
    v_gapped_ref = {}
    v_ungapped_ref = {}
    if args.align is not None:
        v_gapped_ref = simple.read_fasta(args.align)
        for name, seq in v_gapped_ref.items():
            v_ungapped_ref[name] = seq.replace('.', '')
    return v_gapped_ref, v_ungapped_ref


def read_motifs(args, locus):
    if not args.motif_dir and not args.species:
        print('Error - please specify either -motif_dir or -species')
        exit(1)
    if args.motif_dir and args.species:
        print('Error - please specify either -motif_dir or -species, not both')
        exit(1)
    motif_dir = args.motif_dir
    if args.species:
        dm = files('digger')
        motif_dir = dm.joinpath(f'motifs/{args.species}/{locus}')
    print(f'Using motif files from {motif_dir}')
    motifs = {}
    for motif_name in ["J-HEPTAMER", "J-NONAMER", 'L-PART1', 'L-PART2', "V-HEPTAMER", "V-NONAMER"]:
        with open(os.path.join(motif_dir, motif_name + '_prob.csv'), 'r') as fi:
            motifs[motif_name] = Motif(motif_name, stream=fi)
    if locus in ['IGH', 'TRB', 'TRG']:
        for motif_name in ["5'D-HEPTAMER", "5'D-NONAMER", "3'D-HEPTAMER", "3'D-NONAMER"]:
            with open(os.path.join(motif_dir, motif_name + '_prob.csv'), 'r') as fi:
                motifs[motif_name] = Motif(motif_name, stream=fi)
    conserved_motif_seqs = {}
    conserved_motif_file = os.path.join(motif_dir, 'conserved_motifs.fasta')
    if os.path.isfile(conserved_motif_file):
        conserved_motif_seqs = simple.read_fasta(conserved_motif_file)
    return conserved_motif_seqs, motifs

def set_motif_params(locus):
    if locus in ['IGK']:
        J_TRP_MOTIF = 'FGXG'
        J_TRP_OFFSET = 10
        J_SPLICE = 'CGT'
    elif locus in ['IGL']:
        J_TRP_MOTIF = 'FGXG'
        J_TRP_OFFSET = 10
        J_SPLICE = 'GGT'
    else:
        J_TRP_MOTIF = 'WGXG'
        J_TRP_OFFSET = 11
        J_SPLICE = 'GGT'

    if locus not in ['IGK']:
        V_RSS_SPACING = 23
    else:
        V_RSS_SPACING = 12

    return {
        'J_TRP_MOTIF': J_TRP_MOTIF,
        'J_TRP_OFFSET': J_TRP_OFFSET,
        'J_SPLICE': J_SPLICE,
        'V_RSS_SPACING': V_RSS_SPACING
    }

