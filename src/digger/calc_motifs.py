# Create motif matrices from gene features

import argparse
from collections import defaultdict
import random

from digger.motif import Motif
import csv


def get_parser():
    parser = argparse.ArgumentParser(description='Given a set of gene features, create motif matrices')
    parser.add_argument('feat_file', help='feature file, created, for example, by parse_imgt_annotations')
    return parser


feats_in_feat_file = {
    "J-HEPTAMER": 'j-heptamer',
    "J-NONAMER": 'j-nonamer',
    "5'D-HEPTAMER": 'd-5-heptamer',
    "5'D-NONAMER": 'd-5-nonamer',
    "3'D-HEPTAMER": 'd-3-heptamer',
    "3'D-NONAMER": 'd-3-nonamer',
    'L-PART1': 'l-part1',
    'L-PART2': 'l-part2',
    "V-HEPTAMER": 'v-heptamer',
    "V-NONAMER": 'v-nonamer',
}

'''
feats_in_feat_file = {
    "J-HEPTAMER": 'j_heptamer',
    "J-NONAMER": 'j_nonamer',
    "5'D-HEPTAMER": 'd_5_heptamer',
    "5'D-NONAMER": 'd_5_nonamer',
    "3'D-HEPTAMER": 'd_3_heptamer',
    "3'D-NONAMER": 'd_3_nonamer',
    'L-PART1': 'l_part1',
    'L-PART2': 'l_part2',
    "V-HEPTAMER": 'v_heptamer',
    "V-NONAMER": 'v_nonamer',
}
'''

def random_seq(length):
    bases = ['A', 'C', 'G', 'T']
    ret = ''
    for i in range(length):
        ret += bases[random.randrange(4)]
    return ret


def main():
    args = get_parser().parse_args()
    csv.field_size_limit(10000000)      # Reported UTRs can be very long sometimes

    with open(args.feat_file, 'r') as fi:
        reader = csv.DictReader(fi)
        feature_rows = list(reader)

        gene_num = 1

        for feature_name in ["J-HEPTAMER", "J-NONAMER", "5'D-HEPTAMER", "5'D-NONAMER", "3'D-HEPTAMER", "3'D-NONAMER", 'L-PART1', 'L-PART2', "V-HEPTAMER", "V-NONAMER"]:
            seqs = {}
            fname = feats_in_feat_file[feature_name]

            for feature_row in feature_rows:
                if fname in feature_row and len(feature_row[fname]) > 0 and feature_row['functional'].lower() == 'functional':
                    gene_name = feature_row['allele'].split('*')[0] if '*' in feature_row['allele'] else str(gene_num)
                    gene_num += 1
                    if gene_name not in seqs:
                        seqs[gene_name] = []
                    if feature_row[fname].upper() not in seqs[gene_name]:
                        seqs[gene_name].append(feature_row[fname].upper())

            all_seqs = []

            for seq_list in seqs.values():
                for s in seq_list:
                    all_seqs.append(s)

            seqs = all_seqs

            if len(seqs) == 0:
                continue

            lengths = defaultdict(int)
            for seq in seqs:
                lengths[len(seq)] += 1
            max_count = max(lengths.values())
            maj_length = 0
            for k, v in lengths.items():
                if v == max_count:
                    maj_length = k
                    break

            seqs = [list(s) for s in seqs if len(s) == maj_length]
            matrix = Motif(seqs)
            print('%s    %s    %0.4f    %s    %0.2f' % (feature_name, matrix.consensus, matrix.consensus_prob, matrix.conserved_consensus, matrix.calc_prob(matrix.conserved_consensus)))

            probs = []
            length = len(matrix.consensus)
            for i in range(10000):
                probs.append(matrix.calc_prob(random_seq(length))/matrix.consensus_prob)

            probs.sort()
            matrix.likelihood_threshold = probs[int(len(probs)*0.95)]
            print('5pc cutoff: %.2E' % matrix.likelihood_threshold)

            with open(feature_name + '_prob.csv', 'w', newline='') as fo:
                matrix.write_prob_matrix(fo)


if __name__ == "__main__":
    main()
