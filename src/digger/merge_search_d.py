from Bio import SeqIO, pairwise2
import csv
import argparse
import sys
from receptor_utils import simple_bio_seq as simple


# The sequence presented to search-d must be in forward orientation. If the assembly is in reverse sense, the sequence passed
# to search-d will have to be reverse-complemented. In this case, d_locus start should be > d_locus_end, to indicate the orientation

def get_parser():
    parser = argparse.ArgumentParser(description='Merge SEARCH-D results into search file produced by find_alignments.py')
    parser.add_argument('find_alignments_file', help='file produced by find_alignments.py')
    parser.add_argument('search_d_match_file', help='file produced by search-d')
    parser.add_argument('sense', help='sense to be reported for the search-d output (+ or -)')
    parser.add_argument('d_locus_start', help='1-based start co-ord of d locus')
    parser.add_argument('d_locus_end', help='1-based end co-ord of d locus')
    parser.add_argument('entire_locus_end', help='1-based end co-ord of entire locus')
    parser.add_argument('entire_assembly', help='entire assembly sequence passed to find_alignments.py')
    parser.add_argument('-ref', help='reference to compare to: name and reference file separated by comma eg mouse,mouse,fasta (may be repeated multiple times)', action='append')
    parser.add_argument('output_file', help='output file')
    return parser


reference_sets = None

def nt_diff(s1, s2):
    diffs = 0
    row = list(zip(s1.upper(), s2.upper()))

    for p in range(len(row)):
        if row[p][0] != row[p][1]:
            diffs += 1

    diffs += abs(len(s1) - len(s2))
    return diffs


def calc_matched_refs(row):
    for ref in reference_sets:
        row[ref['name'] + '_match'] = ''
        row[ref['name'] + '_score'] = 0
        row[ref['name'] + '_nt_match'] = ''
        row[ref['name'] + '_nt_diffs'] = 999

        for ref_name, ref_seq in ref['seqs'].items():
            if row['gene_type'] in ref_name:
                score = pairwise2.align.globalms(row['seq'], ref_seq, 1, 0, -1, -1, score_only=True)
                if isinstance(score, float) and score > row[ref['name'] + '_score']:
                    row[ref['name'] + '_score'] = score
                    row[ref['name'] + '_match'] = ref_name

                score = nt_diff(row['seq'], ref_seq)
                if score < row[ref['name'] + '_nt_diffs']:
                    row[ref['name'] + '_nt_diffs'] = score
                    row[ref['name'] + '_nt_match'] = ref_name

        #if row[ref['name'] + '_nt_match'] != '':
        #    print('%s / %s\nseq: %s\nref: %s\ndiff: %d\n\n' %
        #      (row['best'], row[ref['name'] + '_nt_match'], row['seq'], ref['seqs'][row[ref['name'] + '_nt_match']], row[ref['name'] + '_nt_diffs']))

        if row[ref['name'] + '_nt_diffs'] == 999:
            row[ref['name'] + '_nt_diffs'] = 0

        if row[ref['name'] + '_score'] > 0:
            row[ref['name'] + '_score'] = round(row[ref['name'] + '_score'] * 100.0 / len(row['seq']), 2)


# co-ordinates of the D locus passsed to SEARCH-D relative to the whole assembly processed by blast
# (sequence that is processed by SEARCH-D is rc compared to the whole assembly)

def main():
    global reference_sets

    args = get_parser().parse_args()

    find_alignments_file = args.find_alignments_file
    search_d_match_file = args.search_d_match_file
    sense = args.sense
    output_file = args.output_file
    d_locus_start = int(args.d_locus_start)
    d_locus_end = int(args.d_locus_end)
    entire_locus_end = int(args.entire_locus_end)
    entire_assembly = simple.read_single_fasta(args.entire_assembly)

    # Reference sets to which we wish to make comparisons
    reference_sets = []
    for ref in args.ref:
        if ',' in ref:
            species, filename = ref.split(',')
            reference_sets.append({'name': species, 'file': filename})
        else:
            print('ref argument should consist of a species name and filename separated by a comma')

    for ref in reference_sets:
        ref['seqs'] = {}
        seqs = SeqIO.parse(ref['file'], 'fasta')
        for seq in seqs:
            ref['seqs'][seq.description] = str(seq.seq).upper()

    bl_recs = []
    with open(find_alignments_file, 'r') as fi_bl:
        reader = csv.DictReader(fi_bl)
        for row in reader:
            bl_recs.append(row)

    sd_recs = SeqIO.parse(search_d_match_file, 'fasta')
    sd_rows = []

    for rec in sd_recs:
        row = {}

        desc = {}
        for item in rec.description.split('|'):
            if ':' not in item:
                print('The FASTA headers do not include RSS co-ordinates. Please use the fork of SEARCH-D at https://github.com/williamdlees/SEARCH-D')
                quit()
            k, v = item.split(':')
            if ',' in v:
                v = v.split(',')
                desc[k] = v[0]
                desc[k + '_start'] = v[1]
            else:
                desc[k] = v

        # turn co-ordinates into integers and reverse them if the sequence passed to d-search is in reverse sense compared to the
        # assembly as a whole

        for k in desc.keys():
            if 'start' in k:
                desc[k] = int(desc[k]) + 1   # adjust for 1-based numbering

        desc['end'] = int(desc['end'])

        # check the sequence of the L9 against the entire assembly to make sure we're ok

        if d_locus_end > d_locus_start:
            ent_L9 = entire_assembly[desc['L9_start'] - 1:desc['L9_start'] + 8]
        else:
            ent_end = entire_locus_end - desc['L9_start']
            ent_L9 = simple.reverse_complement(entire_assembly[ent_end - 8:ent_end + 1])


        if ent_L9 != desc['L9']:
            print(f"Error: in id {desc['id']}, L9 sequence {desc['L9']} does not agree with entire assembly sequence {ent_L9}")

        row = {
            'functional': 'functional',
            'notes': '',
            'seq': str(rec.seq),
            'gene_type': 'IGHD',
            'd_3_heptamer': desc['R7'],
            'd_3_nonamer': desc['R9'],
            'd_5_heptamer': desc['L7'],
            'd_5_nonamer': desc['L9'],
            'sense': sense,
            'contig': bl_recs[0]['contig'],
        }

        def swap(row, a, b):
            row[a], row[b] = row[b], row[a]

        row['start'] = desc['start']
        row['end'] = desc['end']
        row['start_rev'] = entire_locus_end - row['start'] + 1
        row['end_rev'] = entire_locus_end - row['end'] + 1

        row['3_rss_start'] = desc['R7_start']
        row['3_rss_end'] = desc['R9_start'] + 9
        row['3_rss_start_rev'] = entire_locus_end - row['3_rss_start'] + 1
        row['3_rss_end_rev'] = entire_locus_end - row['3_rss_end'] + 1

        row['5_rss_start'] = desc['L9_start']
        row['5_rss_end'] = desc['L7_start'] + 7
        row['5_rss_start_rev'] = entire_locus_end - row['5_rss_start'] + 1
        row['5_rss_end_rev'] = entire_locus_end - row['5_rss_end'] + 1

        if d_locus_end < d_locus_start:
            swap(row, 'start', 'start_rev')
            swap(row, 'end', 'end_rev')

            swap(row, '3_rss_start', '3_rss_start_rev')
            swap(row, '3_rss_end', '3_rss_end_rev')

            swap(row, '5_rss_start', '5_rss_start_rev')
            swap(row, '5_rss_end', '5_rss_end_rev')


        # check the sequence is recorded at the correct location in the full assembly

        if d_locus_end > d_locus_start:
            seq = entire_assembly[row['start'] - 1: row['end']]
        else:
            seq = simple.reverse_complement(entire_assembly[row['end'] - 1:row['start']])

        if rec.seq != seq:
            print(f"Error: in id {desc['id']}, sequence {rec.seq} does not agree with entire assembly sequence {seq}")

        # check L9 is recorded at the correct location in the full assembly

        if d_locus_end > d_locus_start:
            seq_R7 = entire_assembly[row['3_rss_start'] - 1: row['3_rss_start'] + 6]
        else:
            seq_R7 = simple.reverse_complement(entire_assembly[row['3_rss_start'] - 7:row['3_rss_start']])

        if desc['R7'] != seq_R7:
            print(f"Error: in id {desc['id']}, 3' heptamer {desc['R7']} does not agree with entire assembly sequence {seq_R7}")

        # check for overlaps

        add_rec = True
        for sd_row in sd_rows:
            if sd_row['start'] <= row['start'] <= sd_row['end'] or sd_row['start'] <= row['end'] <= sd_row['end']:
                if sd_row['end'] - sd_row['start'] < row['end'] - row['start']:
                    sd_rows.remove(sd_row)
                else:
                    add_rec = False
                    break

        if add_rec:
            calc_matched_refs(row)
            sd_rows.append(row)


    bl_recs.extend(sd_rows)

    fieldnames = bl_recs[0].keys()
    bl_recs.sort(key=lambda x : int(x['start']))

    with open(output_file, 'w', newline='') as fo:
        writer = csv.DictWriter(fo, fieldnames=fieldnames)
        writer.writeheader()
        for rec in bl_recs:
            writer.writerow(rec)


if __name__ == "__main__":
    main()

