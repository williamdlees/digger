# Process a single sequence, given a guide gene to look for
# The sequence may be nominated as a fasta file, or by genbank accession number, or by a list of genbank accession numbers in a CSV file
import argparse
from receptor_utils import simple_bio_seq as simple
import csv

try:
    from dig_utils import set_motif_params, read_motifs, read_vs, process_sequence, process_output, get_genbank_sequence
except:
    from digger.dig_utils import set_motif_params, read_motifs, read_vs, process_sequence, process_output



def get_parser():
    parser = argparse.ArgumentParser(description='Annotate a genomic sequence representing the nominated receptor gene')
    subparsers = parser.add_subparsers(required=True)

    parser_fasta = subparsers.add_parser('fasta', help='Annotate a single sequence in a FASTA file')
    parser_fasta.add_argument('target', help='Name of nominated sequence')
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
    parser_single.add_argument('-align', help='gapped reference set to use for V gene alignments (required for V gene analysis')
    parser_single.add_argument('-species', help='use motifs for the specified species provided with the package')
    parser_single.add_argument('-motif_dir', help='use motif probability files present in the specified directory')
    parser_single.add_argument('-out_file', help='output file (CSV)')
    parser_single.set_defaults(func=main_single)

    parser_multi = subparsers.add_parser('multi', help='Read allele names and corresponding genbank accession numbers from a CSV file')
    parser_multi.add_argument('locus', help='Locus of nominated sequences')
    parser_multi.add_argument('germline_file', help='ungapped reference set containing the nominated sequence (FASTA)')
    parser_multi.add_argument('genbank_acc_file', help='File containing list of targets and associated genbank accession numbers (CSV)')
    parser_multi.add_argument('-align', help='gapped reference set to use for V gene alignments (required for V gene analysis')
    parser_multi.add_argument('-species', help='use motifs for the specified species provided with the package')
    parser_multi.add_argument('-motif_dir', help='use motif probability files present in the specified directory')
    parser_multi.add_argument('-out_file', help='output file (CSV)')
    parser_multi.set_defaults(func=main_multi)


    return parser


def main():
    args = get_parser().parse_args()
    args.func(args)

def main_fasta(args):
    germlines = simple.read_fasta(args.germline_file)
    target = args.target
    locus = target[:3]
    motif_params = set_motif_params(locus)
    conserved_motif_seqs, motifs = read_motifs(args, locus)
    v_gapped_ref, v_ungapped_ref = read_vs(args)

    assembly = list(simple.read_fasta(args.query_file).values())[0]

    genbank_acc = 'NA'
    rows = process_sequence(assembly, genbank_acc, target, germlines, v_gapped_ref, v_ungapped_ref, motifs, conserved_motif_seqs, motif_params)

    process_output(args, rows)


def main_single(args):
    germlines = simple.read_fasta(args.germline_file)
    target = args.target
    locus = target[:3]
    motif_params = set_motif_params(locus)
    conserved_motif_seqs, motifs = read_motifs(args, locus)
    v_gapped_ref, v_ungapped_ref = read_vs(args)

    assembly = get_genbank_sequence(args.genbank_acc)

    genbank_acc = args.genbank_acc
    rows = process_sequence(assembly, genbank_acc, target, germlines, v_gapped_ref, v_ungapped_ref, motifs, conserved_motif_seqs, motif_params)

    process_output(args, rows)


def main_multi(args):
    germlines = simple.read_fasta(args.germline_file)
    locus = args.locus
    motif_params = set_motif_params(locus)
    conserved_motif_seqs, motifs = read_motifs(args, locus)
    v_gapped_ref, v_ungapped_ref = read_vs(args)

    # make a dict of target sequences and associated accession numbers
    genbank_acc_dict = {}
    wanted_accs = []
    with open(args.genbank_acc_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] not in genbank_acc_dict:
                genbank_acc_dict[row[0]] = []
            genbank_acc_dict[row[0]].append(row[1])
            if row[1] not in wanted_accs:
                wanted_accs.append(row[1])

    # Fetch the sequences from genbank and process them
    fetched_accessions = get_genbank_sequence(','.join(wanted_accs))

    # Process the sequences

    rows = []
    for target, accs in genbank_acc_dict.items():
        if target in germlines:
            for genbank_acc in accs:
                if genbank_acc in fetched_accessions:
                    assembly = fetched_accessions[genbank_acc]
                    # print(f'Processing {target} {genbank_acc}')
                    rows.extend(process_sequence(assembly, genbank_acc, target, germlines, v_gapped_ref, v_ungapped_ref, motifs, conserved_motif_seqs, motif_params))
                else:
                    print(f'No sequence found for {genbank_acc}')
        else:
            print(f'No germline sequence found for {target}')

    process_output(args, rows)


if __name__ == "__main__":
    main()
