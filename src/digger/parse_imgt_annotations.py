# this uses the representation at http://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000064

# pull out the V-regions, functional/pseudo and locs

import csv
import os

from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SeqRecord import SeqRecord
from collections import Counter
import argparse
import urllib.request
import html
from receptor_utils import simple_bio_seq as simple

assembly = None

def get_parser():
    parser = argparse.ArgumentParser(description='Given a set of IMGT annotations, build a CSV file containing gene names and co-ordinates')
    parser.add_argument('imgt_url', help='URL of IMGT annotation, e.g. http://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000064, or name of text file containing its contents')
    parser.add_argument('outfile', help='Output file (CSV)')
    parser.add_argument('locus', help='one of IGH, IGK, IGL, TRA, TRB, TRD, TRG')
    parser.add_argument('--save_download', help='Save contents of annotation to specified file')
    parser.add_argument('--save_sequence', help='Save sequence to specified file')
    parser.add_argument('--save_imgt_annots', help='Save IMGT annotations to specified file')
    return parser


def find_range(line):
    if '..' not in line:
        return (None, None)


    (start, end) = line[25:].split('..')
    start = int(start.replace('<', '').replace('>', ''))
    end = int(end.replace('\n', '').replace('<', '').replace('>', ''))
    return (start, end)

def find_seq(line, extra=0):
    sense = 'forward'
    complement = 'complement' in line
    if complement:
        line = line.replace('complement(','').replace(')', '')
        sense = 'reverse'

    (start, end) = find_range(line)

    if start is None:
        return '', sense

    if complement:
        seq = assembly[start-1-extra:end]
    else:
        seq = assembly[start-1:end+extra]

    if complement:
        seq = simple.reverse_complement(seq)

    return seq, sense


def process_V(imgt_annots, parsed_genes, accession):
    in_v_region = False
    parsing_gene = False
    l_part1s = []
    l_part2s = []
    v_heptamers = []
    v_nonamers = []


    row = {'accessions': accession, 'utr': None, 'l-part1': None, 'l-part2': None, 'v-heptamer': None, 'v-nonamer': None, 'v-region': None, 'sense': '', 'v-donor-splice': None, 'v-acceptor-splice': None}

    for line in imgt_annots:
        if 'FT   V-GENE ' in line:
            if parsing_gene:
                # We have a partially parsed truncated gene
                for el in ['utr', 'l-part1', 'l-part2', 'v-heptamer', 'v-nonamer', 'v-region']:
                    row[el] = str(row[el]) if row[el] is not None else ''
                parsed_genes.append(row)
                row = {'accessions': accession, 'utr': None, 'l-part1': None, 'l-part2': None, 'v-heptamer': None, 'v-nonamer': None, 'v-region': None}

        elif 'FT   V-REGION' in line:
            in_v_region = True

            parsing_gene = True
            row['v-region'], row['sense'] = find_seq(line.replace('<', '').replace('>', ''))
        elif len(line) > 5 and line[5] != ' ' and in_v_region:
            in_v_region = False

        if 'FT   L-PART1' in line:
            row['l-part1'], _ = find_seq(line, 2) # include donor splice

        if 'FT   L-PART2' in line:
            row['l-part2'], _ = find_seq(line)

        if "FT   5'UTR" in line:
            row['utr'], _ = find_seq(line)

        if 'FT   V-HEPTAMER' in line:
            row['v-heptamer'], _ = find_seq(line)

        if 'FT   DONOR-SPLICE' in line:
            row['v-donor-splice'], _ = find_seq(line)

        if 'FT   ACCEPTOR-SPLICE' in line:
            row['v-acceptor-splice'], _ = find_seq(line)

        if 'FT   V-NONAMER' in line:
            row['v-nonamer'], _ = find_seq(line)
            l_part1s.append(row['l-part1'])
            l_part2s.append(row['l-part2'])
            v_heptamers.append(row['v-heptamer'])
            v_nonamers.append(row['v-nonamer'])
            for el in ['utr', 'l-part1', 'l-part2', 'v-heptamer', 'v-nonamer', 'v-region']:
                row[el] = str(row[el]) if row[el] is not None else ''
            parsed_genes.append(row)
            row = {'accessions': accession, 'utr': None, 'l-part1': None, 'l-part2': None, 'v-heptamer': None, 'v-nonamer': None, 'v-region': None}
            parsing_gene = False

        if in_v_region:
            if 'V-REGION' in line and '..' in line:
                reg = line.replace('FT   V-REGION', '')
                if 'complement' in reg:
                    reg = reg.replace('complement(', '').replace(')', '')
                reg = reg.split('..')
                row['start'] = reg[0].replace(' ', '')
                row['end'] = reg[1].replace('\n', '')
            elif 'IMGT_allele=' in line:
                row['allele'] = line.split('=')[1].replace('"', '').replace('\n', '')
            elif '/functional' in line:
                row['functional'] = 'functional'
            elif '/pseudo' in line:
                row['functional'] = 'pseudo'
            elif '/ORF' in line:
                row['functional'] = 'ORF'

    if parsing_gene:
        # We have a partially parsed truncated gene
        for el in ['utr', 'l-part1', 'l-part2', 'v-heptamer', 'v-nonamer', 'v-region']:
            row[el] = str(row[el]) if row[el] is not None else ''
        parsed_genes.append(row)

def process_D(imgt_annots, parsed_genes, accession):
    in_d_region = False
    parsing_gene = False

    d_3_heptamers = []
    d_3_nonamers = []
    d_3_spacers = []
    d_5_heptamers = []
    d_5_nonamers = []
    d_5_spacers = []

    row = {'accessions': accession, 'd-3-nonamer': None, 'd-3-spacer': None, 'd-3-heptamer': None, 'd-5-heptamer': None, 'd-5-spacer': None, 'd-5-nonamer': None}

    for line in imgt_annots:
        if 'FT   D-REGION' in line:
            if parsing_gene:
                for el in ['d-3-nonamer', 'd-3-spacer', 'd-3-heptamer', 'd-5-heptamer', 'd-5-spacer', 'd-5-nonamer']:
                    row[el] = str(row[el]) if row[el] is not None else ''
                parsed_genes.append(row)
                row = {'accessions': accession, 'd-3-nonamer': None, 'd-3-spacer': None, 'd-3-heptamer': None, 'd-5-heptamer': None, 'd-5-spacer': None,
                       'd-5-nonamer': None}
            parsing_gene = True
            in_d_region = True
        elif len(line) > 5 and line[5] != ' ' and in_d_region:
            in_d_region = False

        if "FT   5'D-HEPTAMER" in line:
            row['d-5-heptamer'], _ = find_seq(line)
        if "FT   5'D-SPACER" in line:
            row['d-5-spacer'], _ = find_seq(line)
        if "FT   5'D-NONAMER" in line:
            row['d-5-nonamer'], _ = find_seq(line)

        if "FT   3'D-HEPTAMER" in line:
            row['d-3-heptamer'], _ = find_seq(line)
        if "FT   3'D-SPACER" in line:
            row['d-3-spacer'], _ = find_seq(line)

        if "FT   3'D-NONAMER" in line:
            row['d-3-nonamer'], _ = find_seq(line)
            d_3_heptamers.append(row['d-3-heptamer'])
            d_3_nonamers.append(row['d-3-nonamer'])
            d_3_spacers.append(row['d-3-spacer'])
            d_5_heptamers.append(row['d-5-heptamer'])
            d_5_nonamers.append(row['d-5-nonamer'])
            d_5_spacers.append(row['d-5-spacer'])
            for el in ['d-3-nonamer', 'd-3-spacer', 'd-3-heptamer', 'd-5-heptamer', 'd-5-spacer', 'd-5-nonamer']:
                row[el] = str(row[el]) if row[el] is not None else ''
            parsed_genes.append(row)
            row = {'accessions': accession, 'd-3-nonamer': None, 'd-3-spacer': None, 'd-3-heptamer': None, 'd-5-heptamer': None, 'd-5-spacer': None, 'd-5-nonamer': None}
            parsing_gene = False

        if in_d_region:
            if 'FT   D-REGION' in line and '..' in line:
                reg = line.replace('FT   D-REGION', '')
                if 'complement' in reg:
                    reg = reg.replace('complement(', '').replace(')', '')
                reg = reg.split('..')
                row['start'] = reg[0].replace(' ', '')
                row['end'] = reg[1].replace('\n', '')
                row['d-region'], row['sense'] = find_seq(line)
            elif 'IMGT_allele=' in line:
                row['allele'] = line.split('=')[1].replace('"', '').replace('\n', '')
            elif '/functional' in line:
                row['functional'] = 'functional'
            elif '/pseudo' in line:
                row['functional'] = 'pseudo'
            elif '/ORF' in line:
                row['functional'] = 'ORF'

    if parsing_gene:
        for el in ['d-3-nonamer', 'd-3-spacer', 'd-3-heptamer', 'd-5-heptamer', 'd-5-spacer', 'd-5-nonamer']:
            row[el] = str(row[el]) if row[el] is not None else ''
        parsed_genes.append(row)


def process_J(imgt_annots, parsed_genes, accession):
    in_j_region = False
    in_j_gene = False
    j_heptamers = []
    j_nonamers = []
    j_spacers = []

    row = {'accessions': accession, 'start': None, 'end': None, 'allele': None, 'functional': None, 'j-nonamer': None, 'j-spacer': None, 'j-heptamer': None, 'j-region': None}

    for line in imgt_annots:
        if 'FT   J-GENE ' in line:
            in_j_gene = True
        if 'FT   J-REGION ' in line:
            in_j_region = True
        elif len(line) > 5 and line[5] != ' ' and in_j_region:
            in_j_region = False
            in_j_gene = False
            parsed_genes.append(row)
            row = {'accessions': accession, 'start': None, 'end': None, 'allele': None, 'sense': None, 'functional': None, 'j-nonamer': None, 'j-spacer': None, 'j-heptamer': None, 'j-region': None}

        if 'FT   J-HEPTAMER' in line:
            row['j-heptamer'], _ = find_seq(line)
            j_heptamers.append(row['j-heptamer'])
            j_nonamers.append(row['j-nonamer'])
            j_spacers.append(row['j-spacer'])
            for el in ['j-heptamer', 'j-nonamer', 'j-spacer']:
                row[el] = str(row[el]) if row[el] is not None else ''

        if 'FT   J-NONAMER' in line:
            row['j-nonamer'], _ = find_seq(line)

        if 'FT   J-SPACER' in line:
            row['j-spacer'], _ = find_seq(line)


        if in_j_region or in_j_gene:
            if 'FT   J-REGION ' in line and '..' in line:
                reg = line.replace('FT   J-REGION ', '')
                if 'complement' in reg:
                    reg = reg.replace('complement(', '').replace(')', '')
                reg = reg.split('..')
                row['start'] = reg[0].replace(' ', '')
                row['end'] = reg[1].replace('\n', '')
                row['j-region'], row['sense'] = find_seq(line)
            elif 'IMGT_allele=' in line:
                row['allele'] = line.split('=')[1].replace('"', '').replace('\n', '')
            elif '/functional' in line:
                row['functional'] = 'functional'
            elif '/pseudo' in line:
                row['functional'] = 'pseudo'
            elif '/ORF' in line:
                row['functional'] = 'ORF'


def main():
    global assembly
    sequences = {}
    args = get_parser().parse_args()

    if args.locus not in ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG']:
        print('error: locus must be one of IGH, IGK, IGL, TRA, TRB, TRD, TRG')
        exit(0)

    if 'http' in args.imgt_url:
        try:
            with urllib.request.urlopen(args.imgt_url) as fi:
                imgt_text = html.unescape(fi.read().decode('utf-8'))

                if args.save_download:
                    with open(args.save_download, 'w') as fo:
                        fo.write(imgt_text)
        except urllib.error.URLError as e:
            print(f'Error downloading data: {e.reason}')
            print('Please check the URL in a browser. If it is reachable, please try this command again. Otherwise it is likely that the URL is incorrect or the server is down.')
            exit(0)
    else:
        with open(args.imgt_url, 'r') as fi:
            imgt_text = fi.read()

    all_annotations = []
    parsed_genes = []
    annotations = []
    sequence = ""
    reading = 'preamble'
    accessions = ''
    fieldnames = None

    for line in imgt_text.split('\n'):
        if reading == 'preamble':
            annotations.append(line)
            if len(line) > 1 and line[0:2] == 'FT':
                reading = 'annotations'
            if len(line) > 5 and line[0:2] == 'AC':
                accessions = line[5:]
                print(f'Processing {accessions}')

        elif reading == 'annotations':
            annotations.append(line)
            if len(line) < 2 or line[0:2] != 'FT':
                reading = 'postamble'

        elif reading == 'postamble':
            if len(line) > 1 and line[0:2] == 'SQ':
                reading = 'sequence'
            else:
                annotations.append(line)

        elif reading == 'sequence':
            if line[0:5] == '     ':
                sequence += line[5:70].replace(' ', '')
            else:
                assembly = sequence

                if args.save_sequence and accessions:
                    sequences[accessions] = assembly

                process_V(annotations, parsed_genes, accessions)

                if args.locus in ['IGH', 'TRB', 'TRD']:
                    process_D(annotations, parsed_genes, accessions)

                process_J(annotations, parsed_genes, accessions)

                fieldnames = [
                   'accessions', 'allele', 'functional',  'start', 'end', 'sense',
                   'utr', 'l-part1', 'l-part2', 'v-region',  'v-nonamer', 'v-heptamer', 'v-donor-splice', 'v-acceptor-splice',
                   'd-5-heptamer', 'd-5-spacer', 'd-5-nonamer', 'd-region','d-3-heptamer', 'd-3-spacer', 'd-3-nonamer',
                   'j-region', 'j-heptamer', 'j-spacer', 'j-nonamer',
                ]

                sequence = ""
                reading = 'preamble'
                accessions = ''
                all_annotations.extend(annotations)
                annotations = []

    if args.save_imgt_annots:
        with open(args.save_imgt_annots, 'w') as fo:
            fo.write('\n'.join(all_annotations))

    if args.save_sequence:
        simple.write_fasta(args.save_sequence, sequences)

    if fieldnames:
        with open(args.outfile, 'w', newline='') as fo:
            writer = csv.DictWriter(fo, fieldnames=fieldnames, restval='')
            writer.writeheader()
            for parsed_gene in parsed_genes:
                writer.writerow(parsed_gene)
    else:
        print("No annotation records found.")
        print('Please check the URL in a browser and verify that the page contains an IMGT-annotated assembly.')



if __name__ == "__main__":
    main()
