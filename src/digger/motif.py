# A class to manage and use motif matrices

import csv
import json
from receptor_utils import simple_bio_seq as simple


class Motif:
    def __init__(self, name, seqs=None, stream=None):
        self.name = name
        self.matrix = []
        self.consensus = None
        self.conserved_consensus = None
        self.consensus_prob = None
        self.likelihood_threshold = None
        self.logo = None

        if seqs is not None:
            self.calc_prob_matrix(seqs)
        elif stream is not None:
            self.read_prob_matrix(stream)

    def calc_update(self):
        self.calc_consensus()
        self.consensus_prob = self.calc_prob(self.consensus)
        self.calc_conserved_consensus()

    def calc_prob_matrix(self, seqs):
        al = zip(*seqs)

        self.matrix = []
        for z in al:
            pseudo = min(len(seqs)/10, 1)
            dist = {'A': pseudo, 'C': pseudo, 'G': pseudo, 'T': pseudo}  # inc. pseudocount
            total = 4*pseudo
            for base in z:
                if base in ['A', 'C', 'G', 'T']:
                    dist[base] += 1.0
                    total += 1
            for base in ['A', 'C', 'G', 'T']:
                dist[base] = dist[base]/total
            dist_sum = 0
            for base in ['A', 'C', 'G', 'T']:
                dist_sum += dist[base]
            if dist_sum < 0.99 or dist_sum > 1.01:
                print('error in pcm calc')
            self.matrix.append(dist)
        self.calc_update()

    def write_prob_matrix(self, fo):
        fo.write('%.6E\n' % self.likelihood_threshold)
        if self.logo:
            fo.write(f'Logo: {json.dumps(self.logo)}\n')

        writer = csv.writer(fo)
        writer.writerow(['A', 'C', 'G', 'T'])
        for row in self.matrix:
            writer.writerow([round(row[base], 4) for base in ['A', 'C', 'G', 'T']])

    def read_prob_matrix(self, fi):
        self.likelihood_threshold = float(fi.readline().replace('\n', ''))
        pos = fi.tell()
        logo = fi.readline().replace('\n', '')

        if logo and 'Logo: ' in logo:
            self.logo = json.loads(logo.split('Logo: ')[1])
        else:
            fi.seek(pos)

        reader = csv.DictReader(fi)
        self.matrix = list(reader)
        for row in self.matrix:
            for k,v in row.items():
                row[k] = float(v)
        self.calc_update()

    def calc_consensus(self):
        self.consensus = ''
        for row in self.matrix:
            c_pos = None
            c_prob = 0.0
            for base in ['A', 'C', 'G', 'T']:
                if row[base] > c_prob:
                    c_prob = row[base]
                    c_pos = base
            self.consensus += c_pos if c_pos is not None else '-'

    def calc_conserved_consensus(self):
        self.conserved_consensus = ''
        for row in self.matrix:
            c_pos = None
            for base in ['A', 'C', 'G', 'T']:
                if row[base] >= 0.75:
                    c_pos = base
                    break
            self.conserved_consensus += c_pos if c_pos is not None else '-'

    def calc_prob(self, seq):
        prob = 1.0
        for b, m in zip(list(seq),self.matrix):
            if b in m:
                prob = prob * m[b]
        return prob

    def calc_likelihood(self, seq):
        return self.calc_prob(seq) / self.consensus_prob


    def find_non_conserved(self, seq):
        non_conserved = 0
        res = ''
        for s, ref in zip(list(seq), list(self.conserved_consensus)):
            if ref != '-' and s != ref:
                non_conserved += 1
                res += s
            else:
                res += '-'

        return non_conserved, res

    def check_logo(self, seq):
        if self.logo is None:
            return True
        
        seq = simple.translate(seq)
        
        if len(seq) != len(self.logo):
            return None

        ex_count = 0
        for i, base in enumerate(seq):
            if base not in self.logo[i]:
                ex_count += 1

        #print(ex_count)

        #if ex_count == 5:
        #    breakpoint()

        return ex_count < 3
