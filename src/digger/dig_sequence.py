# Process a single sequence, given a guide gene to look for
# The sequence may be nominated as a fasta file, or by genbank accession number, or by a list of genbank accession numbers in a CSV file
import argparse
import csv

from receptor_utils import simple_bio_seq as simple

try:
    from dig_utils import read_motifs, read_vs, process_sequence, process_output, get_genbank_sequence
except:
    from digger.dig_utils import read_motifs, read_vs, process_sequence, process_output, get_genbank_sequence


def get_parser():
    parser = argparse.ArgumentParser(description='Annotate a genomic sequence representing the nominated receptor gene')
    subparsers = parser.add_subparsers(required=True)

    parser_fasta = subparsers.add_parser('fasta', help='Annotate a single genomic sequence in a FASTA file')
    parser_fasta.add_argument('target', help='Name of nominated sequence in reference set')
    parser_fasta.add_argument('germline_file', help='ungapped reference set containing the nominated sequence (FASTA)')
    parser_fasta.add_argument('query_file', help='file containing the sequence to annotate (FASTA)')
    parser_fasta.add_argument('-align', help='gapped reference set to use for V gene alignments (required for V gene analysis')
    parser_fasta.add_argument('-species', help='use motifs for the specified species provided with the package')
    parser_fasta.add_argument('-motif_dir', help='use motif probability files present in the specified directory')
    parser_fasta.add_argument('-out_file', help='output file (CSV)')
    parser_fasta.add_argument('-debug', help='produce parsing_errors file with debug information', action='store_true')
    parser_fasta.set_defaults(func=main_fasta)

    parser_single = subparsers.add_parser('single', help='Annotate a single sequence given its genbank accession number')
    parser_single.add_argument('target', help='Name of nominated sequence')
    parser_single.add_argument('germline_file', help='ungapped reference set containing the nominated sequence (FASTA)')
    parser_single.add_argument('genbank_acc', help='genbank accession number of the sequence to annotate')
    parser_single.add_argument('email_addr', help='email address to provide to genbank')
    parser_single.add_argument('-align', help='gapped reference set to use for V gene alignments (required for V gene analysis')
    parser_single.add_argument('-species', help='use motifs for the specified species provided with the package')
    parser_single.add_argument('-motif_dir', help='use motif probability files present in the specified directory')
    parser_single.add_argument('-out_file', help='output file (CSV)')
    parser_single.set_defaults(func=main_single)

    parser_multi = subparsers.add_parser('multi', help='Read allele names and corresponding genbank accession numbers from a CSV file')
    parser_multi.add_argument('locus', help='Locus of nominated sequences')
    parser_multi.add_argument('germline_file', help='ungapped reference set containing the nominated sequence (FASTA)')
    parser_multi.add_argument('query_file', help='File containing list of targets and associated genbank accession numbers (CSV)')
    parser_multi.add_argument('email_addr', help='email address to provide to genbank')
    parser_multi.add_argument('-align', help='gapped reference set to use for V gene alignments (required for V gene analysis')
    parser_multi.add_argument('-species', help='use motifs for the specified species provided with the package')
    parser_multi.add_argument('-motif_dir', help='use motif probability files present in the specified directory')
    parser_multi.add_argument('-out_file', help='output file (CSV)')
    parser_multi.set_defaults(func=main_multi)

    parser_multi = subparsers.add_parser('multi_seq', help='Read allele names and genomic sequences from a CSV file')
    parser_multi.add_argument('locus', help='Locus of nominated sequences')
    parser_multi.add_argument('germline_file', help='ungapped reference set containing the nominated sequence (FASTA)')
    parser_multi.add_argument('query_file', help='File containing list of targets and genomic sequences (CSV)')
    parser_multi.add_argument('-align', help='gapped reference set to use for V gene alignments (required for V gene analysis')
    parser_multi.add_argument('-species', help='use motifs for the specified species provided with the package')
    parser_multi.add_argument('-motif_dir', help='use motif probability files present in the specified directory')
    parser_multi.add_argument('-out_file', help='output file (CSV)')
    parser_multi.set_defaults(func=main_multi_seq)

    return parser


def main():
    args = get_parser().parse_args()
    args.func(args)

def main_fasta(args):
    germlines = simple.read_fasta(args.germline_file)
    target = args.target
    locus = target[:3]
    conserved_motif_seqs, motifs, motif_params = read_motifs(args, locus)
    v_gapped_ref, v_ungapped_ref = read_vs(args)

    assembly = simple.read_single_fasta(args.query_file)
    patch = ''

    genbank_acc = 'NA'
    row = process_sequence(assembly, genbank_acc, patch, target, germlines, v_gapped_ref, v_ungapped_ref, motifs, conserved_motif_seqs, motif_params)

    process_output(args, [row])


def main_single(args):
    germlines = simple.read_fasta(args.germline_file)
    target = args.target
    locus = target[:3]
    conserved_motif_seqs, motifs, motif_params = read_motifs(args, locus)
    v_gapped_ref, v_ungapped_ref = read_vs(args)

    if target not in germlines:
        print('Target sequence not found in reference set')
        return

    genbank_res = get_genbank_sequence(args.genbank_acc, args.email_addr)

    if args.genbank_acc in genbank_res:
        genbank_res = genbank_res[args.genbank_acc]
        row = process_sequence(genbank_res['seq'], genbank_res['acc'], genbank_res['patch'], target, germlines, v_gapped_ref, v_ungapped_ref, motifs, conserved_motif_seqs, motif_params)
        process_output(args, [row])
    else:
        print('Accession number not found in Genbank')



def main_multi(args):
    multi_driver(args, True)

def main_multi_seq(args):
    multi_driver(args, False)

def multi_driver(args, genomic):
    germlines = simple.read_fasta(args.germline_file)
    locus = args.locus
    conserved_motif_seqs, motifs, motif_params = read_motifs(args, locus)
    v_gapped_ref, v_ungapped_ref = read_vs(args)

    query_list = {}
    if genomic:
        wanted_accs = []
        with open(args.query_file, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0] not in query_list:
                    query_list[row[0]] = []
                query_list[row[0]].append({'accession': row[1]})
                if row[1] not in wanted_accs:
                    wanted_accs.append(row[1])

        # Fetch the sequences from genbank and process them
        genbank_res = get_genbank_sequence(','.join(wanted_accs), args.email_addr)
        accession_results = {acc: row for acc, row in genbank_res.items()}

        # add sequences to query_list
        for allele in query_list:
            for row in query_list[allele]:
                if row['accession'] not in accession_results:
                    print(f'No sequence found for accession number {row["accession"]}')
                else:
                    row['seq'] = accession_results[row['accession']]['seq']
                    row['patch'] = accession_results[row['accession']]['patch']
    else:
        for row in simple.read_csv(args.query_file):
            if row['allele'] not in query_list:
                query_list[row['allele']] = []
            rec = {'seq': row['sequence'], 'accession': row['accession'], 'patch': ''}
            if 'assembly_offset' in row and 'sense' and 'assembly_length' in row:
                rec['assembly_offset'] = row['assembly_offset']
                rec['sense'] = row['sense']
                rec['assembly_length'] = row['assembly_length']
            query_list[row['allele']].append(rec)

    # Process the sequences

    rows = []
    for target_allele, reqs in query_list.items():
        if target_allele in germlines:
            for req in reqs:
                if 'seq' in req and req['seq']:
                    print(f'Processing {target_allele} {req["accession"]}')
                    row = process_sequence(req['seq'], req['accession'], req['patch'], target_allele, germlines, v_gapped_ref, v_ungapped_ref, motifs, conserved_motif_seqs, motif_params)
                    rows.append(row)

                    if row['start'] != 'NA' and 'assembly_offset' in req:
                        row['sense'] = req['sense']
                        if req['sense'] == '+':
                            for el in ['start', 'end', 'gene_start', 'gene_end', 'l_part1_start', 'l_part1_end', 'l_part2_start', 'l_part2_end', '3_rss_start', '3_rss_end', '5_rss_start', '5_rss_end']:
                                if el in row:
                                    row[el] += int(req['assembly_offset'])
                        else:
                            for el in ['start', 'end', 'gene_start', 'gene_end', 'l_part1_start', 'l_part1_end', 'l_part2_start', 'l_part2_end', '3_rss_start', '3_rss_end', '5_rss_start', '5_rss_end']:
                                if el in row:
                                    row[el] = int(req['assembly_length']) - row[el] + 1

                        for el in list(row.keys()):
                            if 'start' in el or 'end' in el and 'rev' not in el:
                                row[el + '_rev'] = int(req['assembly_length']) - row[el] + 1

        else:
            print(f'No germline sequence found for {target_allele}')

    process_output(args, rows)


if __name__ == "__main__":
    main()
