from operator import attrgetter, itemgetter

from Bio import pairwise2
import csv
from collections import defaultdict
from importlib.resources import files

from receptor_utils.number_ighv import number_ighv, gap_nt_from_aa, nt_diff, gap_sequence
from math import ceil
import glob
import os.path
import argparse
from receptor_utils import simple_bio_seq as simple
try:
    from motif import Motif
except:
    from digger.motif import Motif

try:
    from slugify import slugify
except:
    from digger.slugify import slugify


def get_parser():
    parser = argparse.ArgumentParser(description='Find valid genes in a contig given blast matches')
    parser.add_argument('germline_file', help='reference set used to produce the blast matches')
    parser.add_argument('assembly_file', help='assembly or contig provided to blast')
    parser.add_argument('blast_file', help='results from blast in the format provided by blastresults_to_csv (can contain wildcards if there are multiple files, will be matched by glob)')
    parser.add_argument('-species', help='use motifs for the specified species provided with the package')
    parser.add_argument('-motif_dir', help='use motif probability files present in the specified directory')
    parser.add_argument('-ref', help='ungapped reference to compare to: name and reference file separated by comma eg mouse,mouse.fasta (may be repeated multiple times)', action='append')
    parser.add_argument('-align', help='gapped reference file to use for V gene alignments (should contain V genes only), otherwise de novo alignment will be attempted')
    parser.add_argument('-locus', help='locus (default is IGH)')
    parser.add_argument('-sense', help='sense in which to read the assembly (forward or reverse) (will select automatically)')
    parser.add_argument('output_file', help='output file (csv)')
    return parser

# global variables initialised in main()

args = None
germlines = None
locus = None
manual_sense = None
J_TRP_MOTIF = None
J_TRP_OFFSET = None
J_SPLICE = None
V_RSS_SPACING = None
reference_sets = []
v_gapped_ref = {}
v_ungapped_ref = {}
assembly = ''
assembly_rc = ''
assembly_length = 0
motifs = {}
conserved_motif_seqs = {}


def find_all_matches(seq, pattern, thresh=0.7):
    dists = []
    pattern_length = len(pattern)

    significant_bps = len([x for x in pattern if x != 'X'])
    min_result = ceil((1.0-thresh) * significant_bps)

    for i in range(len(seq) - len(pattern)):
        d = nt_diff(seq[i:i+pattern_length], pattern)
        if d <= min_result:
            dists.append((i, d))

    return dists


class DAnnotation:
    def __init__(self, start, end, left_motif, right_motif):
        self.start = start
        self.end = end
        self.seq = assembly[start-1:end]
        self.notes = []
        left_motif.notes = [f"5' {n}" for n in left_motif.notes]
        right_motif.notes = [f"3' {n}" for n in right_motif.notes]
        self.notes.extend(left_motif.notes)
        self.notes.extend(right_motif.notes)
        self.functionality = None
        self.d_5_nonamer = left_motif.left
        self.d_5_heptamer = left_motif.right
        self.d_3_nonamer = right_motif.right
        self.d_3_heptamer = right_motif.left
        self.d_5_rss_start = left_motif.start
        self.d_5_rss_end = left_motif.end
        self.d_3_rss_start = right_motif.start
        self.d_3_rss_end = right_motif.end
        self.likelihood = left_motif.likelihood * right_motif.likelihood

    def annotate(self):
        self.functionality = 'Functional'
        for note in self.notes:
            note = note.lower()
            if 'not found' in note or 'conserved' in note:
                self.functionality = 'ORF'
                break


def process_d(start, end, best, matches):
    left_rss = find_compound_motif("5'D-NONAMER", "5'D-HEPTAMER", 12, 12, 5, end=start-1)
    right_rss = find_compound_motif("3'D-HEPTAMER", "3'D-NONAMER", 12, 12, 5, start=end)

    results = []

    for left in left_rss:
        for right in right_rss:
            results.append(DAnnotation(left.end + 1, right.start - 1, left, right))

    rows = []
    for result in results:
        result.annotate()
        best_match, best_score, best_nt_diffs = calc_best_match_score(best, result.seq)

        row = {
            'start': result.start,
            'end': result.end,
            'end_rev': assembly_length - result.start + 1,
            'start_rev': assembly_length - result.end + 1,
            'evalue': best['evalue'],
            'matches': len(matches),
            'blast_match': best_match,
            'blast_score': best_score,
            'blast_nt_diffs': best_nt_diffs,
            'functional': result.functionality,
            'notes': ', '.join(result.notes),
            'seq': result.seq,
            'd_3_heptamer': result.d_3_heptamer,
            'd_3_nonamer': result.d_3_nonamer,
            'd_5_heptamer': result.d_5_heptamer,
            'd_5_nonamer': result.d_5_nonamer,
            'likelihood': result.likelihood,
        }

        row['3_rss_start'] = result.d_3_rss_start
        row['3_rss_end'] = result.d_3_rss_end
        row['3_rss_start_rev'], row['3_rss_end_rev'] = assembly_length - row['3_rss_end'] + 1, assembly_length - row['3_rss_start'] + 1

        row['5_rss_start'] = result.d_5_rss_start
        row['5_rss_end'] = result.d_5_rss_end
        row['5_rss_start_rev'], row['5_rss_end_rev'] = assembly_length - row['5_rss_end'] + 1, assembly_length - row['5_rss_start'] + 1

        rows.append(row)

    for row in rows:
        check_seq(row['seq'], row['start'], row['end'], False)
        check_seq(row['seq'], row['start_rev'], row['end_rev'], True)
        check_seq(row['d_5_nonamer'], row['5_rss_start'], row['5_rss_start'] + 8, False)
        check_seq(row['d_5_heptamer'], row['5_rss_end'] - 6, row['5_rss_end'], False)
        check_seq(assembly[row['5_rss_start'] - 1:row['5_rss_end']], row['5_rss_start_rev'], row['5_rss_end_rev'], True)
        check_seq(row['d_3_heptamer'], row['3_rss_start'], row['3_rss_start'] + 6, False)
        check_seq(row['d_3_nonamer'], row['3_rss_end'], row['3_rss_end'] - 8, False)
        check_seq(assembly[row['3_rss_start'] - 1:row['3_rss_end']], row['3_rss_start_rev'], row['3_rss_end_rev'], True)

    return rows

class CAnnotation:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.seq = assembly[start-1:end]
        self.notes = []
        self.functionality = None

    def annotate(self):
        self.functionality = 'Functional'
        self.notes = []


def process_c(start, end, best, matches):
    # find optimal starting point for ungapped alignment

    best_match = -1
    best_match_score = 999
    for i in range(start - (int(best['s start'])) - 10, start - (int(best['s start'])) + 10):
        score = simple.nt_diff(assembly[i - 1:i + len(germlines[best['subject']]) - 1], germlines[best['subject']])
        if score < best_match_score:
            best_match_score = score
            best_match = i

    result_start = best_match
    result_end = best_match + len(germlines[best['subject']]) - 1
    result_seq = assembly[best_match - 1:best_match + len(germlines[best['subject']]) - 1]

    best_match, best_score, best_nt_diffs = calc_best_match_score(best, result_seq)

    if best_score < 0.1: # just too many Ns
        return []

    row = {
        'start': result_start,
        'end': result_end,
        'end_rev': assembly_length - result_start + 1,
        'start_rev': assembly_length - result_end + 1,
        'matches': len(matches),
        'blast_match': best_match,
        'blast_score': best_score,
        'blast_nt_diffs': best_nt_diffs,
        'functional': 'Functional',
        'notes': '',
        'seq': result_seq,
    }


    check_seq(row['seq'], row['start'], row['end'], False)
    check_seq(row['seq'], row['start_rev'], row['end_rev'], True)

    return[row]

class JAnnotation:
    def __init__(self, start, end, motif):
        self.start = start
        self.end = end
        self.aa = ''
        self.seq = assembly[start-1:end]
        self.notes = []
        self.functionality = None
        self.j_frame = None
        self.nonamer = motif.left
        self.heptamer = motif.right
        self.rss_start = motif.start
        self.rss_end = motif.end
        self.likelihood = motif.likelihood

    def annotate(self):
        # consider each frame
        hit = 0
        for i in range(0, 3):
            j_codons = self.seq[i:]
            self.aa = simple.translate(j_codons)
            hit = find_best_match(self.aa, J_TRP_MOTIF, thresh=1.0)
            if hit >= 0:
                break

        if hit >= 0:
            self.end = self.start + i + (hit + J_TRP_OFFSET)*3      # include the donor splice
            self.seq = assembly[self.start - 1:self.end]
            if assembly[self.end - 1:self.end + 2] != J_SPLICE:
                self.notes.append('Donor splice not found')
                self.j_frame = 0
                self.functionality = 'pseudo'
            else:
                self.j_frame = i
                self.functionality = 'Functional'
        else:
            self.notes.append('J-TRP not found')
            self.j_frame = 0
            self.functionality = 'pseudo'

def process_j(start, end, best, matches):
    j_rss = find_compound_motif('J-NONAMER', 'J-HEPTAMER', 22, 23, 5, end=start-1)

    if len(j_rss) == 0:
        return []

    # we need to make a call between the rss on the basis of likelihood
    # otherwise we'll get overlapping rss and call js that only differ in a bp or two at the 5' end

    j_rss.sort(key=attrgetter('likelihood'), reverse=True)
    j_rs = j_rss[0]

    result = (JAnnotation(j_rs.end + 1, end, j_rs))
    result.annotate()

    best_match, best_score, best_nt_diffs = calc_best_match_score(best, result.seq)

    row = {
        'start': result.start,
        'end': result.end,
        'end_rev': assembly_length - result.start + 1,
        'start_rev': assembly_length - result.end + 1,
        'evalue': best['evalue'],
        'matches': len(matches),
        'blast_match': best_match,
        'blast_score': best_score,
        'blast_nt_diffs': best_nt_diffs,
        'functional': result.functionality,
        'notes': ', '.join(result.notes),
        'seq': result.seq,
        'j_heptamer': result.heptamer,
        'j_nonamer': result.nonamer,
        'j_frame': result.j_frame,
        'aa': result.aa,
        'likelihood': result.likelihood,
    }

    row['5_rss_start'] = result.rss_start
    row['5_rss_end'] = result.rss_end
    row['5_rss_start_rev'], row['5_rss_end_rev'] = assembly_length - row['5_rss_end'] + 1, assembly_length - row['5_rss_start'] + 1

    check_seq(row['seq'], row['start'], row['end'], False)
    check_seq(row['seq'], row['start_rev'], row['end_rev'], True)
    check_seq(row['j_nonamer'], row['5_rss_start'], row['5_rss_start'] + 8, False)
    check_seq(row['j_heptamer'], row['5_rss_end'] - 6, row['5_rss_end'], False)
    check_seq(assembly[row['5_rss_start'] - 1:row['5_rss_end']], row['5_rss_start_rev'], row['5_rss_end_rev'], True)

    return [row]


class VAnnotation:
    def __init__(self, start, end, refs=None):
        self.start = start
        self.end = end
        self.ungapped = assembly[start-1:end]
        self.gapped = None
        self.gapped_aa = None
        self.notes = []
        self.functionality = None
        self.likelihood = None
        if refs is not None:
            self.align_refs = True
            self.v_gapped_ref = refs[0]
            self.v_ungapped_ref = refs[1]
        else:
            self.align_refs = False

    def annotate(self):
        v_p = simple.translate(self.ungapped)

        if self.align_refs:
            self.gapped, self.gapped_aa, errors = gap_sequence(self.ungapped, self.v_gapped_ref, self.v_ungapped_ref)
        else:
            self.gapped_aa, errors = number_ighv(str(v_p))
            if self.gapped_aa is not None and len(self.gapped_aa) > 0:
                self.gapped = gap_nt_from_aa(self.ungapped, self.gapped_aa)

        if len(errors) > 0:
            if 'stop' in errors.lower():
                self.functionality = 'pseudo'
            else:
                self.functionality = 'ORF'

            self.notes.append(errors)
        else:
            self.functionality = 'Functional'


# assembly is now forward-sense
# want to find all possible leaders within a range, with their probs
# then find all possible rss within range
# evaluate every combination - joint prob and functionality
def process_v(start, end, best, matches, v_parsing_errors):
    #if start == 686835:
    #    breakpoint()

    leaders = find_compound_motif('L-PART1', 'L-PART2', 10, 400, 8, end=start-1, right_force=start-len(motifs['L-PART2'].consensus))

    # restrain length to between 270 and 320 nt, allow a window anywhere within that range
    rights = find_compound_motif('V-HEPTAMER', 'V-NONAMER', V_RSS_SPACING-1, V_RSS_SPACING, 25, start=start+295)

    # put in a dummy leader if we found an RSS but no leader

    if rights and not leaders:
        leaders.append(MotifResult(motifs['L-PART1'],
                                   0,
                                   start-len(motifs['L-PART2'].consensus)-10-len(motifs['L-PART2'].consensus),
                                   motifs['L-PART2'],
                                   0,
                                   start-len(motifs['L-PART2'].consensus),
                                   ['Leader not found']))

    best_rights = find_best_rss(rights)
    best_leaders = find_best_leaders(leaders)

    max_likelihood_rec = None
    results = []

    # report one or more functional V-GENEs found between the identified leaders and rss.
    # if none, report the non-functional sequence with highest combined (leader, rss) likelihood

    if len(best_leaders) == 0:
        v_parsing_errors[start] = ((start, end, best['subject'], 'leader not found', assembly[start-500:start]))

    if len(best_rights) == 0:
        v_parsing_errors[start] = ((start, end, best['subject'], 'rss not found', assembly[end+1:end+24]))

    # if we have multiple possible rights, discard any that are substantially weaker than the best

    if len(best_rights) > 1:
        best_rl = max([x.likelihood for x in best_rights.values()])
        for k in list(best_rights.keys()):
            if best_rights[k].likelihood < (best_rl / 100):
                del best_rights[k]

    for left in best_leaders.values():
        #if left.end == 1114541:
        #    breakpoint()
        for right in best_rights.values():
            if args.align:
                v_annot = VAnnotation(left.end + 1, right.start - 1, refs=[v_gapped_ref, v_ungapped_ref])
            else:
                v_annot = VAnnotation(left.end + 1, right.start - 1)
            v_annot.annotate()
            v_annot.likelihood = left.likelihood * right.likelihood

            if v_annot.functionality == 'Functional':
                results.append((left, v_annot, right))
            if max_likelihood_rec is None or v_annot.likelihood > max_likelihood_rec[1].likelihood:
                max_likelihood_rec = (left, v_annot, right)


    if len(results) == 0 and max_likelihood_rec is not None:
        results = [max_likelihood_rec]

    #if len(results) == 0:
    #    breakpoint()

    rows = []

    for leader, v_gene, rss in results:
        # Mark the gene as non-functional if there is a problem in the leader or RSS.

        if v_gene.functionality != 'pseudo':
            for note in leader.notes:
                note = note.lower()
                if 'not found' in note or 'conserved' in note:
                    v_gene.functionality = 'ORF'
            for note in rss.notes:
                note = note.lower()
                if 'conserved' in note or 'not found' in note:
                    v_gene.functionality = 'ORF'
            for note in leader.notes:
                note = note.lower()
                if 'stop codon' in note or 'atg' in note:
                    v_gene.functionality = 'pseudo'

        best_match, best_score, best_nt_diffs = calc_best_match_score(best, v_gene.ungapped)

        #if '1-56' in best_match:
        #    breakpoint()

        row = {
            'start': v_gene.start,
            'end': v_gene.end,
            'end_rev': assembly_length - v_gene.start + 1,
            'start_rev': assembly_length - v_gene.end + 1,
            'evalue': best['evalue'],
            'matches': len(matches),
            'blast_match': best_match,
            'blast_score': best_score,
            'blast_nt_diffs': best_nt_diffs,
            'l_part1': leader.left,
            'l_part2': leader.right,
            'v_heptamer': rss.left,
            'v_nonamer': rss.right,
            'functional': v_gene.functionality,
            'notes': ', '.join(leader.notes + v_gene.notes + rss.notes),
            'v-gene_aligned_aa': v_gene.gapped_aa,
            'seq': v_gene.ungapped,
            'seq_gapped': v_gene.gapped,
            'likelihood': v_gene.likelihood,
        }

        row['3_rss_start'] = rss.start
        row['3_rss_start_rev'] = assembly_length - rss.end + 1
        row['3_rss_end'] = rss.end
        row['3_rss_end_rev'] = assembly_length - rss.start + 1

        row['l_part1_start'] = leader.start
        row['l_part1_end'] = leader.start + len(leader.left) - 1
        row['l_part1_start_rev'], row['l_part1_end_rev'] = assembly_length - row['l_part1_end'] + 1, assembly_length - row['l_part1_start'] + 1

        row['l_part2_start'] = leader.end - len(leader.right) + 1
        row['l_part2_end'] = leader.end
        row['l_part2_start_rev'], row['l_part2_end_rev'] = assembly_length - row['l_part2_end'] + 1, assembly_length - row['l_part2_start'] + 1

        row['aa'] = simple.translate(v_gene.ungapped)

        rows.append(row)

    for row in rows:
        check_seq(row['seq'], row['start'], row['end'], False)
        check_seq(row['seq'], row['start_rev'], row['end_rev'], True)
        check_seq(row['v_heptamer'], row['3_rss_start'], row['3_rss_start'] + 6, False)
        check_seq(row['v_nonamer'], row['3_rss_end'], row['3_rss_end'] - 8, False)
        check_seq(assembly[row['3_rss_start'] - 1:row['3_rss_end']], row['3_rss_start_rev'], row['3_rss_end_rev'], True)
        check_seq(row['l_part1'], row['l_part1_start'], row['l_part1_end'], False)
        check_seq(assembly[row['l_part1_start'] - 1:row['l_part1_end']], row['l_part1_start_rev'], row['l_part1_end_rev'], True)
        check_seq(row['l_part2'], row['l_part2_start'], row['l_part2_end'], False)
        check_seq(assembly[row['l_part2_start'] - 1:row['l_part2_end']], row['l_part2_start_rev'], row['l_part2_end_rev'], True)

    return rows


# check co-ordinates of each element for sanity
# all co-ordinates should be 1-based!

def check_seq(seq, start, end, rev):
    if len(seq) == 0:
        return  #   allow empty sequences which may occur if the assembly is truncated

    if end < start:
        end, start = start, end

    if end >= len(assembly) or start <= 0:
        return  # allow co-ords to go off the end - can happen if we're extrapolating the leader, say.
                # TODO - fix this, but it's not that big an issue

    if not rev:
        if seq != assembly[start-1:end]:
            print('Error: sequence at (%d, %d) failed co-ordinate check.' % (start, end))
    else:
        if simple.reverse_complement(seq) != assembly_rc[start-1:end]:
            print('Error: reverse sequence at (%d, %d) failed co-ordinate check.' % (start, end))


# Look for a single motif within the specified range of start positions
# return a list of all discovered motifs, and their likelihoods
# min_start, max_start - 1-based co-ords of allowed start positions
def find_single_motif(motif, min_start, max_start):
    res = []
    max_found = 0
    max_hit = ''
    # print('looking for single motif %s in %s' % (motif.consensus, assembly[min_start-1:min_start + len(motif.consensus) + (max_start-min_start)]))
    consensus_len = len(motif.consensus)
    for p in range(min_start, max_start+1):
        seq = assembly[p-1:p-1 + consensus_len]
        likelihood = motif.calc_likelihood(seq)

        if likelihood > max_found:
            max_found = likelihood

        if likelihood >= motif.likelihood_threshold:
            res.append(SingleMotifResult(motif, p, likelihood))

    return res


class SingleMotifResult:
    def __init__(self, motif, position, likelihood):
        self.start = position
        self.name = motif.name
        self.end = position + len(motif.consensus) - 1
        self.seq = assembly[position - 1:position - 1 + len(motif.consensus)]
        self.likelihood = likelihood
        self.notes = []

        if likelihood:      # don't check dummy motifs
            self.check_motif_consensus()


    # check residues that are strongly conserved across all loci
    def check_motif_consensus(self):
        cons = None

        if self.name in conserved_motif_seqs:
            cons = conserved_motif_seqs[self.name]

            non_conserved = 0
            res = ''
            if cons:
                for s, ref in zip(list(self.seq), list(cons)):
                    if ref != '-' and s != ref:
                        res += s
                        non_conserved += 1
                    else:
                        res += '-'

            if non_conserved:
                self.notes.append(f'{self.name} has variation at strongly conserved residue(s): {res}')


# Find compound motifs - ie heptamer plus nonamer, or l-part1 and l-part2.
# left, right - the two motifs
# min_gap, max_gap - the gap between them
# start co-ord, end co-ord: 1-based co-ordinates of either the desired ending position of left, or starting position of right (specify ome or the other)
# window - the function will look either for a left that ends at end +/- window, or a right that starts at start +/- window
# returns [{start, end, left_seq, gap_seq, right_seq, likelihood}] for each left, right combination that meets the threshold criterion
# If force has a value, always report a right starting at that position, with the best available left. This is used to force an L-PART2 to be discovered
# at the point that the BLAST alignment indicates the v-region should start
def find_compound_motif(left_motif_name, right_motif_name, min_gap, max_gap, window, start=None, end=None, right_force=None):
    left_motif = motifs[left_motif_name]
    right_motif = motifs[right_motif_name]
    max_length = len(left_motif.consensus) + max_gap + len(right_motif.consensus)
    min_length = len(left_motif.consensus) + min_gap + len(right_motif.consensus)

    if start is None:
        min_start = end - max_length - window
    else:
        min_start = start - window
    max_start = min_start + (max_length - min_length) + 2 * window

    min_end = min_start + min_length
    max_end = max_start + max_length

    if min_end > len(assembly):
        return []

    if max_end > len(assembly):
        max_end = len(assembly)

    # start by finding the motif with the highest threshold

    res = []

    if left_motif.likelihood_threshold > right_motif.likelihood_threshold:
        # use the specified start position if we have one, otherwise extrapolate from the end and min/max gaps - less accurate
        if start is not None:
            left_results = find_single_motif(left_motif, start - window, start + window)
        else:
            left_results = find_single_motif(left_motif, min_start, min_start + 2*window)

        right_results = []
        for left_result in left_results:
            right_results = find_single_motif(right_motif, left_result.start + len(left_motif.consensus) + min_gap, left_result.start + len(left_result.seq) + max_gap)
            for right_result in right_results:
                res.append(MotifResult(left_result, right_result))

        # put in a dummy if we missed the lower likelihood motif
        if left_results and not right_results:
            best_left = sorted(left_results, key=attrgetter('likelihood'), reverse=True)[0]
            res.append(MotifResult(best_left, SingleMotifResult(right_motif, best_left.start + len(best_left.seq) + min_gap, 0), notes=[f'{right_motif.name} not found']))
    else:
        if end is not None:
            right_results = find_single_motif(right_motif, end - len(right_motif.consensus) - window, end + window)
        else:
            right_results = find_single_motif(right_motif, min_end - len(right_motif.consensus), max_end - len(right_motif.consensus))

        # if right_force is set, if necessary, force discovery of a right motif at the forced position

        if right_force and right_force not in [r.start for r in right_results]:
            right_results.append(SingleMotifResult(right_motif, right_force, 0))

        left_results = []
        for right_result in right_results:
            left_results = find_single_motif(left_motif, right_result.start - max_gap - len(left_motif.consensus), right_result.start - min_gap - len(left_motif.consensus))
            for left_result in left_results:
                res.append(MotifResult(left_result, right_result))

        # put in a dummy if we missed the lower likelihood motif
        if right_results and not left_results:
            best_right = sorted(right_results, key=attrgetter('likelihood'), reverse=True)[0]
            res.append(MotifResult(SingleMotifResult(left_motif, best_right.start - min_gap - len(left_motif.consensus), 0), best_right, notes=[f'{left_motif.name} not found']))

    return res


class MotifResult:
    def __init__(self, left_motif, right_motif, notes=None):
        self.start = left_motif.start
        self.end = right_motif.start + len(right_motif.seq) - 1
        self.left = left_motif.seq
        self.gap = assembly[self.start + len(self.left):right_motif.start - 1]
        self.right = right_motif.seq
        self.likelihood = left_motif.likelihood * right_motif.likelihood
        if notes:
            self.notes = [n for n in notes]
        else:
            self.notes = []
        self.notes.extend(left_motif.notes)
        self.notes.extend(right_motif.notes)



# find the best leader PART1 for each PART2 starting position:
# if PART2 is in-frame: the in-frame PART1 giving best likelihood
# otherwise, the PART1 giving best likelihood regardless of frame
def find_best_leaders(leaders):
    leader_choices = defaultdict(list)
    for leader in leaders:
        leader_choices[leader.end].append(leader)
    best_leaders = {}
    for position, choices in leader_choices.items():
        l2p = simple.translate(choices[0].right[2:])
        choices.sort(key=attrgetter('likelihood'), reverse=True)
        if 'X' not in l2p and '*' not in l2p:
            for choice in choices:
                l12p = simple.translate(choice.left + choice.right)
                if 'X' not in l12p and '*' not in l12p and l12p[0] == 'M':
                    best_leaders[position] = choice
                    break
                else:
                    if 'X' in l12p or '*' in l12p and 'Stop codon in leader' not in choice.notes:
                        choice.notes.append('Stop codon in leader')
                    if l12p[0] != 'M':
                        choice.notes.append('Leader missing initial ATG')
        if position not in best_leaders:
            best_leaders[position] = choices[0]
            if 'Stop codon in leader' not in best_leaders[position].notes:
                best_leaders[position].notes.append('Stop codon in leader')


    return best_leaders


# find the highest-likelihood heptamer/nonamer combination for each starting position
def find_best_rss(rights):
    right_choices = defaultdict(list)
    for right in rights:
        right_choices[right.start].append(right)
    best_rights = {}
    for position, choice in right_choices.items():
        choice.sort(key=attrgetter('likelihood'), reverse=True)
        best_rights[position] = choice[0]

    return best_rights


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

                if 'V' in row['gene_type'] and len(row['seq']) > 290 and len(ref_seq) > 290:
                    score = nt_diff(row['seq'][:291], ref_seq[:291])
                else:
                    score = nt_diff(row['seq'], ref_seq)
                if score < row[ref['name'] + '_nt_diffs']:
                    row[ref['name'] + '_nt_diffs'] = score

        #if row[ref['name'] + '_nt_match'] != '':
        #    print('%s / %s\nseq: %s\nref: %s\ndiff: %d\n\n' %
        #      (row['best'], row[ref['name'] + '_nt_match'], row['seq'], ref['seqs'][row[ref['name'] + '_nt_match']], row[ref['name'] + '_nt_diffs']))

        if row[ref['name'] + '_nt_diffs'] == 999:
            row[ref['name'] + '_nt_diffs'] = 0

        if row[ref['name'] + '_score'] > 0:
            row[ref['name'] + '_score'] = round(row[ref['name'] + '_score'] * 100.0 / len(row['seq']), 2)


def calc_best_match_score(best, seq):
    score = pairwise2.align.globalms(seq, germlines[best['subject']], 1, 0, -1, -1, score_only=True)
    if isinstance(score, float) and score > 0:
        best_score = round(100 * score / len(seq), 2)
        best_match = best['subject']
    else:
        best_score = 0.0
        best_match = ''

    best_nt_diffs = nt_diff(seq, germlines[best['subject']])

    return best_match, best_score, best_nt_diffs


def strings_to_num(row, fields):
    for field in fields:
        if '.' in row[field]:
            row[field] = float(row[field])
        else:
            row[field] = int(row[field])
    return row


def find_best_match(seq, pattern, thresh=0.7):
    dists = find_all_matches(seq, pattern, thresh)

    if len(dists) == 0:
        return -1

    scores = [d[1] for d in dists]

    return dists[scores.index(min(scores))][0]


def reverse_coords(row):
    start_revs = [x for x in row.keys() if 'start_rev' in x]
    for start_rev in start_revs:
        end_rev = start_rev.replace('start', 'end')
        start = start_rev.replace('_rev', '')
        end = end_rev.replace('_rev', '')
        row[start], row[end_rev] = row[end_rev], row[start]
        row[end], row[start_rev] = row[start_rev], row[end]


def process_file(this_blast_file, writer):
    global assembly
    global assembly_rc
    global J_TRP_MOTIF
    global J_TRP_OFFSET
    global J_SPLICE

    assembly_rc = simple.reverse_complement(assembly)

    with open(this_blast_file, 'r') as fi:
        print('Processing %s' % this_blast_file)
        this_contig_name = None
        forward_gene_alignments = defaultdict(list)
        reverse_gene_alignments = defaultdict(list)

        reader = csv.DictReader(fi, dialect='excel')
        for row in reader:
            if not this_contig_name:
                this_contig_name = row['query']
            if row['subject'] not in germlines:
                print('%s is not defined in %s. Please make sure the file has definitions for all sequences used by BLAST.' % (row['subject'], args.germline_file))
                quit()

            row = strings_to_num(row, ['identity', 'alignment length', 'mismatches', 'gap opens', 'q start', 'q end', 's start', 's end', 'evalue', 'bit score'])
            if row['s end'] < row['s start']:
                row['sense'] = '-'
                row['q start'], row['q end'] = assembly_length - row['q end'] + 1, assembly_length - row['q start'] + 1
                row['s start'], row['s end'] = row['s end'], row['s start']
                row['q start'] -= row['s start'] - 1
                row['q end'] = row['q start'] + len(germlines[row['subject']]) -1
                reverse_gene_alignments[row['q end']].append(row)
            else:
                row['sense'] = '+'
                row['q start'] -= row['s start'] - 1
                row['q end'] = row['q start'] + len(germlines[row['subject']]) -1
                forward_gene_alignments[row['q end']].append(row)

        # calculate mean of the best evalue for hits at each position

        def mean_score(gene_alignment):
            total = 0.0
            count = 0

            for end, matches in gene_alignment.items():
                scores = [match['evalue'] if abs(match['s start'] - match['s end']) > 0.9 * len(germlines[match['subject']]) else 100 for match in matches]
                if min(scores) < 100:
                    total += min(scores)
                    count += 1

            return count, (total / count) if count > 0 else 99

        print('%d forward alignments, mean evalue %0.2g' % mean_score(forward_gene_alignments))
        print('%d reverse alignments, mean evalue %0.2g' % mean_score(reverse_gene_alignments))

        if manual_sense:
            print(f"Using {manual_sense} sense by request")

        co_ordinates_reversed = False

        if manual_sense == 'forward' or (manual_sense is None and len(forward_gene_alignments) > len(reverse_gene_alignments)):
            gene_alignments = forward_gene_alignments
            print('Using forward alignments')
        else:
            gene_alignments = reverse_gene_alignments
            assembly = simple.reverse_complement(assembly)
            assembly_rc = simple.reverse_complement(assembly)
            co_ordinates_reversed = True
            print('Using reverse alignments')

        results = {}

        # the best match for a location is a complete match against a reference sequence
        # (see overlap reconciliation at the bottom of the loop)
        # otherwise the lowest e-val for a reasonable length
        def match_score(match):
            if match['s end'] - match['s start'] + 1 == len(germlines[match['subject']]) and match['mismatches'] == '0':
                return 0
            elif match['s end'] - match['s start'] + 1 > 0.4 * len(germlines[match['subject']]):
                return match['evalue']
            else:
                return 100

        v_parsing_errors = {}

        for end, matches in gene_alignments.items():
            notes = []
            scores = [match_score(match) for match in matches]
            best = matches[scores.index(min(scores))]

            # make sure we take an alignment length of 90% or more if available, otherwise flag a warning

            #if min(scores) > 99:
            #    notes.append('BLAST alignment was truncated')

            if locus + 'C' not in best['subject']:
                gene_type = locus + best['subject'].split(locus)[1][0]
            else:
                gene_type = locus + 'C'

            start = best['q start']

            #print('processing match to %s at %d' % (best['subject'], best['q start']))

            add_rows = True

            if 'V' in gene_type:
                # don't process obviously very truncated records
                if abs(end - start) < 250:
                    continue
                rows = process_v(start, end, best, matches, v_parsing_errors)
            elif 'J' in gene_type:
                rows = process_j(start, end, best, matches)
            elif 'D' in gene_type:
                rows = process_d(start, end, best, matches)
            elif 'C' in gene_type:
                rows = process_c(start, end, best, matches)
            else:
                add_rows = False

            if add_rows:
                for row in rows:
                    row['gene_type'] = gene_type

                    if len(notes) > 0:
                        if len(row['notes']) > 0:
                            row['notes'] = ', '.join(row['notes'].split(', ') + notes)
                        else:
                            row['notes'] = ', '.join(notes)

                    add_rec = True
                    for k, result in list(results.items()):
                        # if this row overlaps with one already stored in results, keep the one that is longer provided it is functional and there are no 'not found' elements
                        if result['start'] <= row['start'] <= result['end'] or result['start'] <= row['end'] <= result['end']:
                            # if row['likelihood'] > result['likelihood']:

                            #if 'not found' in results[k]['notes'] or 'not found' in row['notes']:
                            #    breakpoint()

                            if row['functional'] == 'Functional' \
                                    and (row['end'] - row['start'] > result['end'] - result['start'] or result['functional'] != 'Functional')\
                                    and ('not found' in results[k]['notes'] or 'not found' not in row['notes']):
                                del(results[k])
                            else:
                                add_rec = False
                                break

                    if add_rec:
                        row['sense'] = '-' if co_ordinates_reversed else '+'
                        if co_ordinates_reversed:
                            reverse_coords(row)
                        results[row['start']] = row

        for row in results.values():
            calc_matched_refs(row)
            row['contig'] = this_contig_name

        for key in sorted(results.keys()):
            writer.writerow(results[key])

        with open(this_blast_file.replace('.', '_v_parsing_errors.'), 'w', newline='') as fo:
            writer = csv.writer(fo)
            for start in sorted(v_parsing_errors.keys()):
                writer.writerow(v_parsing_errors[start])


def main():
    global args, assembly, assembly_length, germlines, locus, manual_sense, J_TRP_MOTIF, J_TRP_OFFSET,  J_SPLICE, V_RSS_SPACING, reference_sets, conserved_motif_seqs
    global v_gapped_ref

    args = get_parser().parse_args()
    assembly_file = args.assembly_file
    blast_file = args.blast_file
    match_file = args.output_file
    germlines = simple.read_fasta(args.germline_file)
    locus = args.locus if args.locus is not None else 'IGH'

    if args.sense:
        if args.sense in ['forward', 'reverse']:
            manual_sense = args.sense
        else:
            print("Error - sense must be 'forward' or 'reverse'")
            exit(0)

    if locus in ['IGK']:
        J_TRP_MOTIF = 'FGXG'
        J_TRP_OFFSET = 10
        J_SPLICE = 'CGT'
    elif locus in ['IGL']:
        J_TRP_MOTIF = 'FGXG'
        J_TRP_OFFSET = 10
        J_SPLICE = 'GGT'
    else:
        J_TRP_MOTIF = 'WGXG'
        J_TRP_OFFSET = 11
        J_SPLICE = 'GGT'

    if locus not in ['IGK']:
        V_RSS_SPACING = 23
    else:
        V_RSS_SPACING = 12

    if not args.motif_dir and not args.species:
        print('Error - please specify either -motif_dir or -species')
        exit(1)

    if args.motif_dir and args.species:
        print('Error - please specify either -motif_dir or -species, not both')
        exit(1)

    motif_dir = args.motif_dir

    if args.species:
        motif_dir = files('digger.motifs').joinpath(f'{args.species}/{locus}')

    print(f'Using motif files from {motif_dir}')

    for motif_name in ["J-HEPTAMER", "J-NONAMER", 'L-PART1', 'L-PART2', "V-HEPTAMER", "V-NONAMER"]:
        with open(os.path.join(motif_dir, motif_name + '_prob.csv'), 'r') as fi:
            motifs[motif_name] = Motif(motif_name, stream=fi)

    if locus in ['IGH', 'TRB', 'TRG']:
        for motif_name in ["5'D-HEPTAMER", "5'D-NONAMER", "3'D-HEPTAMER", "3'D-NONAMER"]:
            with open(os.path.join(motif_dir, motif_name + '_prob.csv'), 'r') as fi:
                motifs[motif_name] = Motif(motif_name, stream=fi)

    conserved_motif_file = os.path.join(motif_dir, 'conserved_motifs.fasta')
    if os.path.isfile(conserved_motif_file):
        conserved_motif_seqs = simple.read_fasta(conserved_motif_file)

    for ref in args.ref:
        if ',' in ref:
            species, filename = ref.split(',')
            reference_sets.append({'name': species, 'file': filename})
        else:
            print('ref argument should consist of a species name and filename separated by a comma')

    for ref in reference_sets:
        ref['seqs'] = simple.read_fasta(ref['file'])

    if args.align is not None:
        v_gapped_ref = simple.read_fasta(args.align)
        for name, seq in v_gapped_ref.items():
            v_ungapped_ref[name] = seq.replace('.', '')

    with open(match_file, 'w', newline='') as fo:
        fieldnames = ['contig', 'start', 'end', 'start_rev', 'end_rev', 'sense']
        for ref in reference_sets:
            fieldnames.extend([ref['name'] + '_match', ref['name'] + '_score', ref['name'] + '_nt_match', ref['name'] + '_nt_diffs'])

        fieldnames.extend(
            ['likelihood', 'l_part1', 'l_part2', 'v_heptamer', 'v_nonamer', 'j_heptamer', 'j_nonamer', 'j_frame', 'd_3_heptamer', 'd_3_nonamer', 'd_5_heptamer', 'd_5_nonamer',
            'functional', 'notes', 'aa', 'v-gene_aligned_aa', 'seq', 'seq_gapped', '5_rss_start', '5_rss_start_rev', '5_rss_end', '5_rss_end_rev',
            '3_rss_start', '3_rss_start_rev', '3_rss_end', '3_rss_end_rev', 'l_part1_start', 'l_part1_start_rev', 'l_part1_end', 'l_part1_end_rev',
            'l_part2_start', 'l_part2_start_rev', 'l_part2_end', 'l_part2_end_rev'])

        fieldnames.extend(['matches', 'gene_type', 'blast_match', 'blast_score', 'blast_nt_diffs', 'evalue'])

        writer = csv.DictWriter(fo, fieldnames=fieldnames, restval='')
        writer.writeheader()

        if '*' not in blast_file:
            assembly = simple.read_single_fasta(assembly_file)
            assembly_length = len(assembly)
            process_file(blast_file, writer)
        else:
            assemblies = {}
            for name, seq in simple.read_fasta(assembly_file).items():
                assemblies[name] = seq

            for blast_file_name in glob.glob(blast_file):
                suffix = os.path.splitext(os.path.basename(blast_file_name))[0].replace(blast_file.split('*')[0], '')
                assembly = None
                for assembly_name in assemblies.keys():
                    s_name = slugify(assembly_name.split(' ')[0].replace('\n', '').replace('|', '')) + '.csv'
                    if s_name in blast_file_name or suffix in assembly_name:
                        assembly = assemblies[assembly_name]
                        assembly_length = len(assembly)
                        break

                if assembly is not None:
                    process_file(blast_file_name, writer)
                else:
                    print('no corresponding assembly found for output file %s' % blast_file_name)

if __name__ == "__main__":
    main()
