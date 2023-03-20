# this uses the representation at http://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000064

# pull out the V-regions, functional/pseudo and locs

import csv
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
    parser.add_argument('imgt_url', help='URL of IMGT annotation, e.g. http://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000064, of name of text file containing its contents')
    parser.add_argument('outfile', help='Output file (CSV)')
    parser.add_argument('locus', help='IGH, IGK or IGL')
    parser.add_argument('--save_download', help='Save download to specified file')
    parser.add_argument('--save_sequence', help='Save sequence to specified file')
    parser.add_argument('--save_imgt_annots', help='Save IMGT annotations to specified file')
    parser.add_argument('--summarise_motifs', help='Print motifs to stdout', action='store_true')
    return parser


def find_range(line):
    if '..' not in line:
        return (None, None)


    (start, end) = line[25:].split('..')
    start = int(start)
    end = int(end.replace('\n', ''))
    return (start, end)

def find_seq(line, extra=0):
    sense = 'forward'
    complement = 'complement' in line
    if complement:
        line = line.replace('complement(','').replace(')', '')
        sense = 'reverse'

    (start, end) = find_range(line)

    if start is None:
        return None

    seq = assembly[start-1:end+extra]

    if complement:
        seq = simple.reverse_complement(seq)

    return seq, sense


def process_V(imgt_annots, parsed_genes, summarise_motifs):
    in_v_region = False
    l_part1s = []
    l_part2s = []
    v_heptamers = []
    v_nonamers = []


    row = {'utr': None, 'l-part1': None, 'l-part2': None, 'v-heptamer': None, 'v-nonamer': None, 'v-region': None, 'sense': ''}

    for line in imgt_annots:
        if 'FT   V-REGION' in line:
            in_v_region = True
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

        if 'FT   V-NONAMER' in line:
            row['v-nonamer'], _ = find_seq(line)
            l_part1s.append(row['l-part1'])
            l_part2s.append(row['l-part2'])
            v_heptamers.append(row['v-heptamer'])
            v_nonamers.append(row['v-nonamer'])
            for el in ['utr', 'l-part1', 'l-part2', 'v-heptamer', 'v-nonamer', 'v-region']:
                row[el] = str(row[el]) if row[el] is not None else ''
            parsed_genes.append(row)
            row = {'utr': None, 'l-part1': None, 'l-part2': None, 'v-heptamer': None, 'v-nonamer': None, 'v-region': None}

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

    if summarise_motifs:
        for (name, val) in zip(['l_part1', 'l_part2', 'v_heptamer', 'v_nonamer'], [l_part1s, l_part2s, v_heptamers, v_nonamers]):
            # chuck out any weird lengths
            lengths = Counter([len(s) for s in val if s is not None])
            l = lengths.most_common()[0][0]
            val = [SeqRecord(v) for v in val if v is not None and len(v) == l]
            alignment = MultipleSeqAlignment(val)
            summary = AlignInfo.SummaryInfo(alignment)
            consensus = summary.dumb_consensus(threshold=0.5)
            for v in val:
                print(str(v.seq))
            print('%s consensus: %s' % (name, str(consensus)))
        print('l-part1 consensus includes ending donor splice')


def process_D(imgt_annots, parsed_genes, summarise_motifs):
    in_d_region = False

    d_3_heptamers = []
    d_3_nonamers = []
    d_3_spacers = []
    d_5_heptamers = []
    d_5_nonamers = []
    d_5_spacers = []

    row = {'d-3-nonamer': None, 'd-3-spacer': None, 'd-3-heptamer': None, 'd-5-heptamer': None, 'd-5-spacer': None, 'd-5-nonamer': None}

    for line in imgt_annots:
        if 'FT   D-REGION' in line:
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
            row = {'d-3-nonamer': None, 'd-3-spacer': None, 'd-3-heptamer': None, 'd-5-heptamer': None, 'd-5-spacer': None, 'd-5-nonamer': None}

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

    if summarise_motifs:
        for (name, val) in zip(['d-3-nonamer', 'd-3-spacer', 'd-3-heptamer', 'd-5-heptamer', 'd-5-spacer', 'd-5-nonamer'],
                               [d_3_nonamers, d_3_spacers, d_3_heptamers, d_5_heptamers, d_5_spacers, d_5_nonamers]):
            # chuck out any weird lengths
            lengths = Counter([len(s) for s in val if s is not None])
            l = lengths.most_common()[0][0]
            val = [SeqRecord(v) for v in val if v is not None and len(v) == l]
            alignment = MultipleSeqAlignment(val)
            summary = AlignInfo.SummaryInfo(alignment)
            consensus = summary.dumb_consensus(threshold=0.5)
            for v in val:
                print(str(v.seq))
            print('%s consensus: %s' % (name, str(consensus)))


def process_J(imgt_annots, parsed_genes, summarise_motifs):
    in_j_region = False
    in_j_heptamer = False
    in_j_nonamer = False

    j_heptamers = []
    j_nonamers = []
    j_spacers = []

    row = {'start': None, 'end': None, 'allele': None, 'functional': None, 'j-nonamer': None, 'j-spacer': None, 'j-heptamer': None, 'j-region': None}

    for line in imgt_annots:
        if 'FT   J-REGION ' in line:
            in_j_region = True
        elif len(line) > 5 and line[5] != ' ' and in_j_region:
            in_j_region = False
            parsed_genes.append(row)
            row = {'start': None, 'end': None, 'allele': None, 'sense': None, 'functional': None, 'j-nonamer': None, 'j-spacer': None, 'j-heptamer': None, 'j-region': None}

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


        if in_j_region:
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

    if summarise_motifs:
        for (name, val) in zip(['j_heptamer', 'j-spacer', 'j_nonamer'], [j_heptamers, j_spacers, j_nonamers]):
            # chuck out any weird lengths
            lengths = Counter([len(s) for s in val if s is not None])
            l = lengths.most_common()[0][0]
            val = [SeqRecord(v) for v in val if v is not None and len(v) == l]
            alignment = MultipleSeqAlignment(val)
            summary = AlignInfo.SummaryInfo(alignment)
            consensus = summary.dumb_consensus(threshold=0.5)
            for v in val:
                print(str(v.seq))
            print('%s consensus: %s' % (name, str(consensus)))


def main():
    global assembly
    args = get_parser().parse_args()

    if args.locus not in ['IGH', 'IGK', 'IGL']:
        print('error: locus must me IGH, IGK or IGL')
        exit(0)

    if 'http' in args.imgt_url:
        with urllib.request.urlopen(args.imgt_url) as fi:
            imgt_text = html.unescape(fi.read().decode('utf-8'))

            if args.save_download:
                with open(args.save_download, 'w') as fo:
                    fo.write(imgt_text)
    else:
        with open(args.imgt_url, 'r') as fi:
            imgt_text = fi.read()

    annotations = []
    sequence = ""
    reading = 'preamble'

    for line in imgt_text.split('\n'):
        if reading == 'preamble':
            if len(line) > 1 and line[0:2] == 'FT':
                annotations.append(line)
                reading = 'annotations'

        elif reading == 'annotations':
            if len(line) < 2 or line[0:2] != 'FT':
                reading = 'postamble'
            else:
                annotations.append(line)

        elif reading == 'postamble':
            if len(line) > 1 and line[0:2] == 'SQ':
                reading = 'sequence'

        elif reading == 'sequence':
            if line[0:5] == '     ':
                sequence += line[5:70].replace(' ', '')
            else:
                break

    assembly = sequence

    if args.save_sequence:
        simple.write_fasta({args.locus: sequence}, args.save_sequence)

    if args.save_imgt_annots:
        with open(args.save_imgt_annots, 'w') as fo:
            fo.write('\n'.join(annotations))

    parsed_genes = []
    process_V(annotations, parsed_genes, args.summarise_motifs)

    if args.locus == 'IGH':
        process_D(annotations, parsed_genes, args.summarise_motifs)

    process_J(annotations, parsed_genes, args.summarise_motifs)

    fieldnames = [
       'allele', 'functional',  'start', 'end', 'sense',
       'utr', 'l-part1', 'l-part2', 'v-region',  'v-nonamer', 'v-heptamer',
       'd-5-heptamer', 'd-5-spacer', 'd-5-nonamer', 'd-region','d-3-heptamer', 'd-3-spacer', 'd-3-nonamer',
       'j-region', 'j-heptamer', 'j-spacer', 'j-nonamer'
    ]

    with open(args.outfile, 'w', newline='') as fo:
        writer = csv.DictWriter(fo, fieldnames=fieldnames, restval='')
        writer.writeheader()
        for parsed_gene in parsed_genes:
            writer.writerow(parsed_gene)


if __name__ == "__main__":
    main()
