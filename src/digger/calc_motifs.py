# Create motif matrices from gene features

import argparse
import os
from collections import defaultdict
import random
from importlib.resources import files
from shutil import copyfile
import csv
import statistics
from Bio.Cluster import distancematrix, kcluster, treecluster
import numpy as np
from receptor_utils import simple_bio_seq as simple

try:
    from motif import Motif
except:
    from digger.motif import Motif



def get_parser():
    parser = argparse.ArgumentParser(description='Given a set of gene features, create motif matrices')
    parser.add_argument('locus', help='locus (e.g. IGH, TRA)')
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
    "OCTAMER": 'octamer',   
}


def random_seq(length):
    bases = ['A', 'C', 'G', 'T']
    ret = ''
    for i in range(length):
        ret += bases[random.randrange(4)]
    return ret


def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def distance_matrix(seqs):
    n = len(seqs)
    d = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d[i, j] = hamming_distance(seqs[i], seqs[j])
            d[j, i] = d[i, j]
    return d


def aa_consensus(aas):
    consensus = list(aas[0])
    for i in range(1, len(aas)):
        row = aas[i]
        for j in range(len(consensus)):
            if consensus[j] != row[j]:
                consensus[j] = '-'
    return ''.join(consensus)


def aa_logo(aas):
    consensus = [list(a) for a in list(aas[0])]
    for i in range(1, len(aas)):
        row = aas[i]
        for j in range(len(consensus)):
            if row[j] not in consensus[j]:
                consensus[j].append(row[j])
    return consensus


def main():
    args = get_parser().parse_args()

    accepted_loci = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG']
    if args.locus not in accepted_loci:
        print(f"locus must be one of {', '.join(accepted_loci)}")
        return

    try:
        dm = files('digger')
        motif_dir = dm.joinpath(f'motifs/human/{args.locus}')
    except TypeError:
        path = os.path.abspath(__file__)
        dm = os.path.dirname(path)
        motif_dir = os.path.join(dm, f'motifs/human/{args.locus}')

    copyfile(f'{motif_dir}/motif_params.json', './motif_params.json')

    csv.field_size_limit(10000000)      # Reported UTRs can be very long sometimes

    with open(args.feat_file, 'r') as fi:
        reader = csv.DictReader(fi)
        feature_rows = list(reader)

        gene_num = 1

        for feature_name in ["J-HEPTAMER", "J-NONAMER", "5'D-HEPTAMER", "5'D-NONAMER", "3'D-HEPTAMER", "3'D-NONAMER", 'L-PART1', 'L-PART2', "V-HEPTAMER", "V-NONAMER", "OCTAMER"]:
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

            if feature_name == 'L-PART1':
                seqs = [s for s in seqs if s.startswith('ATG')]

            print(f'Processing {len(seqs)} sequences for motif {feature_name}')

            if len(seqs) == 0:
                if feature_name != 'OCTAMER':
                    print(f'WARNING: no sequences found for {feature_name}. Digger will not function correctly without representative sequences for this feature.')
                continue

            lengths = [len(s) for s in seqs]    

            # for L-PART1, only consider sequences that start with ATG

            # determine the lengths of all sequences for this motif

            maj_length = int(statistics.mode(lengths))

            # for all except L-PART1 and L-PART2, only consider sequences of the most common length

            res = []
            if feature_name in ['L-PART1', 'L-PART2']:
                lengths = list(set(lengths))
            else:
                lengths = [maj_length]

            # see whether we should split any length into multiple groups

            for length in lengths:
                lseqs = list(set([s for s in seqs if len(s) == length]))
                if feature_name == 'L-PART1':
                    logo = aa_logo([simple.translate(l) for l in lseqs])
                else:
                    logo = None

                if len(lseqs) > 2 and feature_name in ['L-PART1', 'L-PART2']:
                    dm = distance_matrix(lseqs)
                    for nk in range(1, 5):
                        if nk > 1:
                            clusterid, error, nfound = kcluster(dm, nclusters=nk)
                            #tree = treecluster(None, distancematrix=dm)
                            #clusterid = tree.cut(nk)

                            # split sequences into clusters by clusterid
                            clusters = defaultdict(list)
                            for i, c in enumerate(clusterid):
                                clusters[c].append(lseqs[i])
                        else:
                            clusters = {0: lseqs}
                        
                        # calculate the maximum distance within each cluster
                        dists = []
                        for c in clusters:
                            mdm = distance_matrix(clusters[c])
                            #print(f'cluster {c} has {len(clusters[c])} sequences max dist {np.amax(mdm)}')
                            dists.append(np.amax(mdm))

                        if max(dists) <= 12:
                            break

                    for k, v in clusters.items():
                        res.append((v, logo))
                else:
                    res.append((lseqs, logo))

            cl_ind = 1
            for seqs, logo in res:
                if not len(seqs):
                    continue 

                matrix = Motif(feature_name, seqs)
                print('%s    %s    %0.4f    %s    %0.2f' % (feature_name, matrix.consensus, matrix.consensus_prob, matrix.conserved_consensus, matrix.calc_prob(matrix.conserved_consensus)))

                probs = []
                length = len(matrix.consensus)
                for i in range(10000):
                    probs.append(matrix.calc_prob(random_seq(length))/matrix.consensus_prob)

                probs.sort()
                matrix.likelihood_threshold = probs[int(len(probs)*0.95)]
                print('5pc cutoff: %.2E' % matrix.likelihood_threshold)

                if feature_name in ['L-PART1', 'L-PART2']:
                    fn = feature_name + '_' + str(cl_ind) + '_' + str(len(seqs[0]))
                    cl_ind += 1
                else:
                    fn = feature_name

                with open(fn + '_prob.csv', 'w', newline='') as fo:
                    matrix.logo = logo
                    matrix.write_prob_matrix(fo)


if __name__ == "__main__":
    main()
