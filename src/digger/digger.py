#!python
# Annotate a reference assembly given a starting germline set

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE

import argparse
import os.path
import subprocess
from receptor_utils import simple_bio_seq as simple
import pathlib
import glob
from digger.slugify import slugify


def get_parser():
    parser = argparse.ArgumentParser(description='Find functional and nonfunctional genes in a single assembly sequence')
    parser.add_argument('assembly_file', help='assembly sequence to search')
    parser.add_argument('motif_dir', help='pathname to directory containing motif probability files')
    parser.add_argument('-locus', help='locus (default is IGH)')
    parser.add_argument('-v_ref', help='set of V reference genes to use as starting point for search')
    parser.add_argument('-d_ref', help='set of D reference genes to use as starting point for search')
    parser.add_argument('-j_ref', help='set of J reference genes to use as starting point for search')
    parser.add_argument('-v_ref_gapped', help='IMGT-gapped v-reference set used to determine alignment of novel sequences')
    parser.add_argument('-ref_comp', help='ungapped reference set(s) to compare to: name and reference file separated by comma eg mouse,mouse.fasta (may be repeated multiple times)', action='append')
    parser.add_argument('-sense', help='sense in which to read the assembly (forward or reverse) (if omitted will select automatically)')
    parser.add_argument('-searchd', help='directory in which SEARCHD can be found (if specified, this will be used for finding D genes in place of blastn)')
    parser.add_argument('-keepwf', help='keep working files after processing has completed', action='store_true')
    parser.add_argument('output_file', help='output file (csv)')
    return parser


# A function to remove intermediate files

def remove_working_files():
    for fn in ['full_germline_set.fasta', 'assembly.fasta', 'assembly_rc.fasta', 'blast_results_assembly.csv', 'results_ex_searchd.csv', 'searchd_results.fasta']:
        if os.path.isfile(fn):
            os.remove(fn)

    for fn in glob.glob('assembly*.out'):
        if os.path.isfile(fn):
            os.remove(fn)

    for ext in ['.ndb', '.nhr', '.nin', '.njs', '.not', '.nsq', '.ntf', '.nto']:
        for fn in glob.glob('*' + ext):
            if os.path.isfile(fn):
                os.remove(fn)


def main():
    args = get_parser().parse_args()

    cwd = pathlib.Path().resolve()
    searchd_path = None
    remove_working_files()

    if os.path.isfile(args.output_file):
        os.remove(args.output_file)

    # check arguments

    if args.searchd:
        searchd_path = os.path.join(args.searchd, 'search_d.py')
        if not os.path.isfile(searchd_path):
            print(f"{searchd_path} - file not found")
            exit(1)

    for fn in [args.assembly_file, args.v_ref]:
        if not os.path.isfile(fn):
            print(f"{fn} - file not found")
            exit(1)

    if not os.path.isdir(args.motif_dir):
        print(f"{args.motif_dir} - file not found")
        exit(1)

    full_germline_set = {}
    for fn in [args.v_ref, args.d_ref, args.j_ref]:
        if fn and not os.path.isfile(fn):
            print(f"{fn} - file not found")
            exit(1)

        if fn:
            full_germline_set |= simple.read_fasta(fn)

    simple.write_fasta(full_germline_set, 'full_germline_set.fasta')

    if args.ref_comp:
        for ref_arg in args.ref_comp:
            ref_name, fn = ref_arg.split(',')
            if not os.path.isfile(fn):
                print(f"{fn} - file not found")
                exit(1)

    if not(args.v_ref or args.d_ref or args.j_ref):
        print('Please specify at least one starting point reference set')
        exit(1)

    assembly_contents = simple.read_fasta(args.assembly_file)

    if len(assembly_contents) != 1:
        print('The assembly file must contain exactly one sequence.')
        exit(1)

    simple.write_fasta({list(assembly_contents.keys())[0]: list(assembly_contents.values())[0]}, 'assembly.fasta')
    assembly_name = slugify(list(assembly_contents.keys())[0])
    assembly_length = len(list(assembly_contents.values())[0])

    locus = 'IGH'

    if args.locus:
        if args.locus not in ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG']:
            print(f"locus must be one of {','.join(['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG'])}")
            exit(1)
        locus = args.locus

    # Make blast databases

    for fn in [args.v_ref, args.d_ref, args.j_ref]:
        if fn:
            cmd = [
                    'makeblastdb',
                    '-in', f'{fn}',
                    '-dbtype', 'nucl',
            ]
            print(f"\n-- executing {' '.join(cmd)}\n")
            process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            #if not os.path.isfile(fn+'.ndb'):
            #    print('Blast database not created - quitting')
            #    exit(1)

    # Run blast and process output

    first_time = True

    for fn, gene_type in [(args.v_ref, 'V'), (args.d_ref, 'D'), (args.j_ref, 'J')]:
        if fn:
            if gene_type == 'D' and searchd_path:
                print('skipping BLAST of D genes as SEARCHD has been specified')
                continue

            if gene_type == 'J':
                word_size = '7'
            else:
                word_size = '11'

            cmd = [
                'blastn',
                '-db', f'{fn}',
                '-query', 'assembly.fasta',
                '-out', f"assembly_{gene_type}.out",
                '-outfmt', '7',
                '-gapopen', '5',
                '-gapextend', '5',
                '-penalty', '-1',
                '-word_size', word_size,
            ]
            print(f"\n-- executing {' '.join(cmd)}\n")
            process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            if not os.path.isfile(f"assembly_{gene_type}.out"):
                print('Blast output not created - quitting')
                exit(1)

            cmd = [
                'blastresults_to_csv',
                f"assembly_{gene_type}.out",
                'blast_results_',
            ]

            if not first_time:
                cmd.append('-a')

            first_time = False

            print(f"\n-- executing {' '.join(cmd)}\n")
            process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    if not os.path.isfile(f'blast_results_{assembly_name}.csv'):
        print('Blast results not converted to csv - quitting')
        exit(1)

    # Run find_alignments

    cmd = [
        'find_alignments',
        'full_germline_set.fasta',
        'assembly.fasta',
        f'blast_results_{assembly_name}.csv',
        args.motif_dir,
        '-locus', locus,
        ]

    if args.ref_comp:
        for ref_arg in args.ref_comp:
            cmd.extend(['-ref', ref_arg])

    if args.sense:
        cmd.extend(['-sense', args.sense])

    if args.v_ref_gapped:
        cmd.extend(['-align', args.v_ref_gapped])

    cmd.append(args.output_file)

    print(f"\n-- executing {' '.join(cmd)}\n")
    process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(process.stdout.decode("utf-8"))
    reverse_sense = 'Using reverse' in process.stdout.decode("utf-8")

    if not os.path.isfile(args.output_file):
        print(f'Result file {args.output_file} was not produced by find_alignments.py')
        exit(1)

    if not args.searchd:
        if not args.keepwf:
            remove_working_files()
        exit(0)

    # Run SEARCH-D if requested

    os.rename(args.output_file, 'results_ex_searchd.csv')

    if reverse_sense:
        simple.write_fasta({'assembly_rc': simple.reverse_complement(list(assembly_contents.values())[0])}, 'assembly_rc.fasta')
        abs_assembly = os.path.abspath('assembly_rc.fasta')
    else:
        abs_assembly = os.path.abspath('assembly.fasta')

    abs_res = os.path.abspath('searchd_results.fasta')
    os.chdir(args.searchd)

    cmd = [
        'python',
        'search_d.py',
        abs_assembly,
        abs_res,
    ]

    print(f"\n-- executing {' '.join(cmd)}\n")
    process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    os.chdir(cwd)

    if not os.path.isfile('searchd_results.fasta'):
        print('SEARCH-D results not found - quitting')
        exit(1)

    if reverse_sense:
        cmd = [
            'merge_search_d',
            'results_ex_searchd.csv',
            'searchd_results.fasta',
            '-',
            f"{assembly_length}", '1', f"{assembly_length}",
            'assembly.fasta',
            args.output_file,
        ]
    else:
        cmd = [
            'merge_search_d',
            'results_ex_searchd.csv',
            'searchd_results.fasta',
            '+',
            '1', f"{assembly_length}", f"{assembly_length}",
            'assembly.fasta',
            args.output_file,
        ]

    if args.ref_comp:
        for ref_arg in args.ref_comp:
            cmd.extend(['-ref', ref_arg])

    print(f"\n-- executing {' '.join(cmd)}\n")
    process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    if not args.keepwf:
        remove_working_files()

if __name__ == "__main__":
    main()


