from operator import attrgetter, itemgetter

from Bio.Align import PairwiseAligner
import csv
from collections import defaultdict
from importlib.resources import files

from receptor_utils.number_v import number_ighv, gap_nt_from_aa, nt_diff, gap_sequence
from math import ceil
import glob
import os.path
import argparse
from receptor_utils import simple_bio_seq as simple

try:
    from search_motifs import find_compound_motif, MotifResult, SingleMotifResult
    from parse_genes import process_v, process_d, process_j, process_c
    from dig_utils import reverse_coords, read_motifs
except Exception as e:
    from digger.search_motifs import find_compound_motif, MotifResult, SingleMotifResult
    from digger.parse_genes import process_v, process_d, process_j, process_c
    from digger.dig_utils import reverse_coords, read_motifs

try:
    from motif import Motif
except:
    from digger.motif import Motif

try:
    from slugify import slugify
except:
    from digger.slugify import slugify


def get_parser():
    parser = argparse.ArgumentParser(description='Find valid genes in a contig given blast matches')
    parser.add_argument('germline_file', help='reference set used to produce the blast matches')
    parser.add_argument('assembly_file', help='assembly or contig provided to blast')
    parser.add_argument('blast_file', help='results from blast in the format provided by blastresults_to_csv (can contain wildcards if there are multiple files, will be matched by glob)')
    parser.add_argument('-species', help='use motifs for the specified species provided with the package')
    parser.add_argument('-motif_dir', help='use motif probability files present in the specified directory')
    parser.add_argument('-ref', help='ungapped reference to compare to: name and reference file separated by comma eg mouse,mouse.fasta (may be repeated multiple times)', action='append')
    parser.add_argument('-align', help='gapped reference file to use for V gene alignments (should contain V genes only), otherwise de novo alignment will be attempted')
    parser.add_argument('-locus', help='locus (default is IGH)')
    parser.add_argument('-sense', help='sense in which to read the assembly (forward or reverse) (will select automatically)')
    parser.add_argument('-debug', help='produce parsing_errors file with debug information', action='store_true')
    parser.add_argument('output_file', help='output file (csv)')
    return parser

# global variables initialised in main()

args = None
germlines = None
locus = None
manual_sense = None
J_TRP_MOTIF = None
J_TRP_OFFSET = None
J_SPLICE = None
V_RSS_SPACING = None
J_RSS_SPACING = None
reference_sets = []
v_gapped_ref = {}
v_ungapped_ref = {}
assembly = ''
assembly_rc = ''
assembly_length = 0
motifs = {}
conserved_motif_seqs = {}
motif_params = {}



aligner_global = PairwiseAligner(
    mode = 'global',
    open_gap_score = -1,
    extend_gap_score = -1,
    match_score = 1,
    mismatch_score = 0
)

def calc_matched_refs(row):
    for ref in reference_sets:
        row[ref['name'] + '_match'] = ''
        row[ref['name'] + '_score'] = 0
        row[ref['name'] + '_nt_diffs'] = 999
        best_ref_seq = ''

        for ref_name, ref_seq in ref['seqs'].items():
            if row['gene_type'] in ref_name:
                if row['seq'] and ref_seq:
                    alignment = aligner_global.align(row['seq'], ref_seq)[0]
                    score = alignment.score
                
                    if score > row[ref['name'] + '_score']:
                        row[ref['name'] + '_score'] = score
                        row[ref['name'] + '_match'] = ref_name
                        c = alignment.counts()
                        row[ref['name'] + '_nt_diffs'] = max(len(row['seq']), len(ref_seq)) - c.identities
                        best_ref_seq = ref_seq

        if row[ref['name'] + '_nt_diffs'] == 999:
            row[ref['name'] + '_nt_diffs'] = 0

        if row[ref['name'] + '_score'] > 0:
            row[ref['name'] + '_score'] = round(row[ref['name'] + '_score'] * 100.0 / len(row['seq']), 2)

        # if we get a 100% score, check that there are no diffs and that the two strings are identical
        if row[ref['name'] + '_score'] == 100:
            if row[ref['name'] + '_nt_diffs'] != 0:
                print(f'WARNING: 100% match but nt_diff is nonzero: {row}')
            if row['seq'] != best_ref_seq:
                print(f'WARNING: 100% match reported but sequences are not identical: {row}')
        if row[ref['name'] + '_match'] and row[ref['name'] + '_nt_diffs'] == 0 and row[ref['name'] + '_score'] != 100:
            print(f'WARNING: nt_diff is zero but score is not 100%: {row}')


def strings_to_num(row, fields):
    for field in fields:
        if '.' in row[field]:
            row[field] = float(row[field])
        else:
            row[field] = int(row[field])
    return row


def update_note(base, note):
    if base['notes']:
        if note not in base['notes']:
            base['notes'] = ','.join(base['notes'].split(',') + [note])
    else:
        base['notes'] = note


def process_file(this_blast_file, writer, write_parsing_errors):
    global assembly
    global assembly_rc
    global J_TRP_MOTIF
    global J_TRP_OFFSET
    global J_SPLICE

    assembly_rc = simple.reverse_complement(assembly)

    with open(this_blast_file, 'r') as fi:
        print('Processing %s' % this_blast_file)
        this_contig_name = None
        forward_gene_alignments = defaultdict(list)
        reverse_gene_alignments = defaultdict(list)

        reader = csv.DictReader(fi, dialect='excel')
        for row in reader:
            if not this_contig_name:
                this_contig_name = row['query']
            if row['subject'] not in germlines:
                print('%s is not defined in %s. Please make sure the file has definitions for all sequences used by BLAST.' % (row['subject'], args.germline_file))
                quit()

            row = strings_to_num(row, ['identity', 'alignment length', 'mismatches', 'gap opens', 'q start', 'q end', 's start', 's end', 'evalue', 'bit score'])

            if row['s end'] < row['s start']:
                row['sense'] = '-'
                row['q start'], row['q end'] = assembly_length - row['q end'] + 1, assembly_length - row['q start'] + 1
                row['s start'], row['s end'] = row['s end'], row['s start']
                row['q start'] -= row['s start'] - 1
                row['q end'] = row['q start'] + len(germlines[row['subject']]) -1
                reverse_gene_alignments[row['q end']].append(row)
            else:
                row['sense'] = '+'
                row['q start'] -= row['s start'] - 1
                row['q end'] = row['q start'] + len(germlines[row['subject']]) -1
                forward_gene_alignments[row['q end']].append(row)

        # calculate mean of the best evalue for hits at each position

        def mean_score(gene_alignment):
            total = 0.0
            count = 0

            for end, matches in gene_alignment.items():
                scores = [match['evalue'] if abs(match['s start'] - match['s end']) > 0.9 * len(germlines[match['subject']]) else 100 for match in matches]
                if min(scores) < 100:
                    total += min(scores)
                    count += 1

            return count, (total / count) if count > 0 else 99

        fw_count, fw_score = mean_score(forward_gene_alignments)
        rev_count, rev_score = mean_score(reverse_gene_alignments)

        print('%d forward alignments, mean evalue %0.2g' % (fw_count, fw_score))
        print('%d reverse alignments, mean evalue %0.2g' % (rev_count, rev_score))

        if manual_sense:
            print(f"Using {manual_sense} sense by request")

        co_ordinates_reversed = False

        if manual_sense == 'forward' or (manual_sense is None and fw_count > rev_count):
            gene_alignments = forward_gene_alignments
            print('Using forward alignments')
        else:
            gene_alignments = reverse_gene_alignments
            assembly = simple.reverse_complement(assembly)
            assembly_rc = simple.reverse_complement(assembly)
            co_ordinates_reversed = True
            print('Using reverse alignments')

        results = {}

        # the best match for a location is a complete match against a reference sequence
        # (see overlap reconciliation at the bottom of the loop)
        # otherwise the lowest e-val for a reasonable length
        def match_score(match):
            if match['s end'] - match['s start'] + 1 == len(germlines[match['subject']]) and match['mismatches'] == '0':
                return 0
            elif match['s end'] - match['s start'] + 1 > 0.4 * len(germlines[match['subject']]):
                return match['evalue']
            else:
                return 100

        v_parsing_errors = {}

        for end, matches in gene_alignments.items():
            notes = []
            scores = [match_score(match) for match in matches]
            best = matches[scores.index(min(scores))]

            # make sure we take an alignment length of 90% or more if available, otherwise flag a warning

            #if min(scores) > 99:
            #    notes.append('BLAST alignment was truncated')

            if locus not in best['subject']:
                print(f"Locus {locus} not found in allele name {best['subject']} - was it specified correctly?")
                quit(1)

            if locus + 'C' not in best['subject']:
                gene_type = locus + best['subject'].split(locus)[1][0]
            else:
                gene_type = locus + 'C'

            start = best['q start']

            # print('processing match to %s at %d' % (best['subject'], best['q start']))

            add_rows = True

            if 'V' in gene_type:
                # don't process obviously very truncated records
                if abs(end - start) < 250 or start < 1 or end > assembly_length:
                    continue
                rows = process_v(assembly, assembly_rc, germlines, v_gapped_ref, v_ungapped_ref, conserved_motif_seqs, motifs, motif_params, start, end, best, matches, args.align, V_RSS_SPACING, v_parsing_errors)
            elif 'J' in gene_type:
                rows = process_j(assembly, assembly_rc, germlines, conserved_motif_seqs, motifs, start, end, best, matches, J_TRP_MOTIF, J_TRP_OFFSET, J_SPLICE, J_RSS_SPACING)
            elif 'D' in gene_type:
                rows = process_d(assembly, assembly_rc, germlines, conserved_motif_seqs, motifs, start, end, best, matches, D_5_RSS_SPACING, D_3_RSS_SPACING)
            elif 'C' in gene_type:
                rows = process_c(assembly, assembly_rc, germlines, start, end, best, matches)
            else:
                add_rows = False

            if add_rows:
                for row in rows:
                    row['gene_type'] = gene_type

                    if len(notes) > 0:
                        if len(row['notes']) > 0:
                            row['notes'] = ', '.join(row['notes'].split(', ') + notes)
                        else:
                            row['notes'] = ', '.join(notes)

                    add_rec = True
                    if row['start'] < row['end']:
                        row_length = row['end'] - row['start'] + 1
                    else:
                        row_length = row['start'] - row['end'] + 1
                    for k, result in list(results.items()):
                        if result['start'] < result['end']:
                            result_length = result['end'] - result['start'] + 1
                        else:
                            result_length = result['start'] - result['end'] + 1

                        # if this row overlaps with one already stored in results, keep the one that is functional
                        # if both are functional or both nonfunctional, keep the one that most closely matches a reference, if we have reference sets
                        # Otherwise keep the longer sequence

                        if row['start'] < row['end'] and (result['start'] <= row['start'] <= result['end'] or result['start'] <= row['end'] <= result['end']):
                            replace = False

                            if row['functional'] != 'Functional' and result['functional'] == 'Functional':
                                replace = False
                            elif row['functional'] == 'Functional' and result['functional'] != 'Functional':
                                replace = True
                            else:
                                # if we're comparing against one or more references, use the one with the best score

                                if reference_sets:
                                    if reference_sets[0]['name'] + '_score' not in result:
                                        calc_matched_refs(result)
                                    calc_matched_refs(row)

                                    max_score = 0
                                    for ref in reference_sets:
                                        if result[ref['name'] + '_score'] > max_score:
                                            max_score = result[ref['name'] + '_score']
                                            replace = False
                                        if row[ref['name'] + '_score'] > max_score:
                                            max_score = row[ref['name'] + '_score']
                                            replace = True
                                else:
                                    replace = row_length > result_length

                                base = row if replace else result
                                alt = result if replace else row

                                if row['functional'] == 'Functional' and result['functional'] == 'Functional':
                                    if result['start'] != row['start']:
                                        update_note(base, f'alternate start: {alt["start"]}')
                                    if result['end'] != row['end']:
                                        update_note(base, f'alternate end: {alt["end"]}')

                            if replace:
                                del(results[k])
                            else:
                                add_rec = False
                                break

                    if add_rec:
                        row['sense'] = '+'
                        results[row['start']] = row

        if co_ordinates_reversed:
            reversed_results = {}
            for row in results.values():
                reverse_coords(row)
                row['sense'] = '-'
                reversed_results[row['start']] = row
            results = reversed_results

        for row in results.values():
            calc_matched_refs(row)
            row['contig'] = this_contig_name

        for key in sorted(results.keys()):
            writer.writerow(results[key])

        if write_parsing_errors:
            with open(this_blast_file.replace('.', '_v_parsing_errors.'), 'w', newline='') as fo:
                writer = csv.writer(fo)
                for start in sorted(v_parsing_errors.keys()):
                    writer.writerow(v_parsing_errors[start])


def main():
    global args, assembly, assembly_length, germlines, locus, manual_sense, J_TRP_MOTIF, J_TRP_OFFSET,  J_SPLICE, J_RSS_SPACING, V_RSS_SPACING, D_5_RSS_SPACING, D_3_RSS_SPACING, vreference_sets, conserved_motif_seqs
    global v_gapped_ref
    global motifs
    global motif_params

    args = get_parser().parse_args()
    assembly_file = args.assembly_file
    blast_file = args.blast_file
    match_file = args.output_file
    germlines = simple.read_fasta(args.germline_file)
    locus = args.locus if args.locus is not None else 'IGH'

    if args.sense:
        if args.sense in ['forward', 'reverse']:
            manual_sense = args.sense
        else:
            print("Error - sense must be 'forward' or 'reverse'")
            exit(0)

    conserved_motif_seqs, motifs, motif_params = read_motifs(args, locus)
    J_TRP_MOTIF = motif_params['J_TRP_MOTIF']
    J_TRP_OFFSET = motif_params['J_TRP_OFFSET']
    J_SPLICE = motif_params['J_SPLICE']
    V_RSS_SPACING = motif_params['V_RSS_SPACING']
    J_RSS_SPACING = motif_params['J_RSS_SPACING']
    D_5_RSS_SPACING = motif_params['D_5_RSS_SPACING']
    D_3_RSS_SPACING = motif_params['D_3_RSS_SPACING']


    if args.ref:
        for ref in args.ref:
            if ',' in ref:
                species, filename = ref.split(',')
                reference_sets.append({'name': species, 'file': filename})
            else:
                print('ref argument should consist of a species name and filename separated by a comma')

    for ref in reference_sets:
        ref['seqs'] = simple.read_fasta(ref['file'])

    if args.align is not None:
        v_gapped_ref = simple.read_fasta(args.align)
        for name, seq in v_gapped_ref.items():
            v_ungapped_ref[name] = seq.replace('.', '')

    with open(match_file, 'w', newline='') as fo:
        fieldnames = ['contig', 'start', 'end', 'start_rev', 'end_rev', 'sense', 'gene_type', 'gene_start', 'gene_end', 'gene_start_rev', 'gene_end_rev']
        for ref in reference_sets:
            fieldnames.extend([ref['name'] + '_match', ref['name'] + '_score', ref['name'] + '_nt_diffs'])

        fieldnames.extend(
            ['functional', 'notes', 'likelihood', 'tata_box', 'octamer', 'l_part1', 'l_part2', 'v_heptamer', 'v_nonamer', 'j_heptamer', 'j_nonamer', 'j_frame', 'd_3_heptamer', 'd_3_nonamer', 'd_5_heptamer', 'd_5_nonamer',
            'aa', 'v-gene_aligned_aa', 'gene_seq', 'seq', 'seq_gapped', '5_rss_start', '5_rss_start_rev', '5_rss_end', '5_rss_end_rev',
            '3_rss_start', '3_rss_start_rev', '3_rss_end', '3_rss_end_rev', 'l_part1_start', 'l_part1_start_rev', 'l_part1_end', 'l_part1_end_rev',
            'l_part2_start', 'l_part2_start_rev', 'l_part2_end', 'l_part2_end_rev', 'octamer_start', 'octamer_start_rev', 'octamer_end', 'octamer_end_rev',
         'tata_box_start', 'tata_box_start_rev', 'tata_box_end', 'tata_box_end_rev',])

        fieldnames.extend(['matches', 'blast_match', 'blast_score', 'blast_nt_diffs', 'evalue'])

        writer = csv.DictWriter(fo, fieldnames=fieldnames, restval='')
        writer.writeheader()

        if '*' not in blast_file:
            assembly = simple.read_single_fasta(assembly_file)
            assembly_length = len(assembly)
            process_file(blast_file, writer, args.debug)
        else:
            assemblies = {}
            for name, seq in simple.read_fasta(assembly_file).items():
                assemblies[name] = seq

            for blast_file_name in glob.glob(blast_file):
                suffix = os.path.splitext(os.path.basename(blast_file_name))[0].replace(blast_file.split('*')[0], '')
                assembly = None
                for assembly_name in assemblies.keys():
                    s_name = slugify(assembly_name.split(' ')[0].replace('\n', '').replace('|', '')) + '.csv'
                    if s_name in blast_file_name or suffix in assembly_name:
                        assembly = assemblies[assembly_name]
                        assembly_length = len(assembly)
                        break

                if assembly is not None:
                    process_file(blast_file_name, writer, args.debug)
                else:
                    print('no corresponding assembly found for output file %s' % blast_file_name)

if __name__ == "__main__":
    main()
