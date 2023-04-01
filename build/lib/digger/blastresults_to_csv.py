# Convert blast file format 7 to one or more CSVs

import argparse
import os
from digger.slugify import slugify


def get_parser():
    parser = argparse.ArgumentParser(description='Convert blast file format 7 to one or more CSVs')
    parser.add_argument('infile', help='the blast file')
    parser.add_argument('out_prefix', help='prefix for csv files')
    parser.add_argument('-a', '--append', help='append to existing output files', action='store_true')
    return parser


def main():
    args = get_parser().parse_args()

    cols = ['query', 'subject', 'identity', 'alignment length', 'mismatches', 'gap opens', 'q start', 'q end', 's start', 's end', 'evalue', 'bit score']

    fo = None

    with open(args.infile, 'r') as fi:
        for row in fi:
            if '# Query: ' in row:
                if fo is not None:
                    fo.close()
                fn = slugify(args.out_prefix + row.replace('# Query: ', '').split(' ')[0].replace('\n', '').replace('|', '')) + '.csv'
                if args.append and os.path.isfile(fn):
                    fo = open(fn, 'a')
                else:
                    fo = open(fn, 'w')
                    fo.write(','.join(cols) + '\n')
            elif row[0] != '#':
                fo.write(row.replace('\t', ','))


if __name__ == "__main__":
    main()
