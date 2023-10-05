from collections import defaultdict
from math import ceil
from operator import attrgetter

from Bio.Align import PairwiseAligner
from receptor_utils import simple_bio_seq as simple
from receptor_utils.number_v import gap_sequence, number_ighv, gap_nt_from_aa

try:
    from motif import Motif
    from search_motifs import find_compound_motif, MotifResult, SingleMotifResult
except:
    from digger.motif import Motif
    from digger.search_motifs import find_compound_motif, MotifResult, SingleMotifResult


class DAnnotation:
    """
    Represents an annotation for a D-segment.

    :param assembly: The assembly sequence.
    :type assembly: str
    :param start: The starting position of the D-segment.
    :type start: int
    :param end: The ending position of the D-segment.
    :type end: int
    :param left_motif: The left motif information.
    :type left_motif: MotifResult
    :param right_motif: The right motif information.
    :type right_motif: MotifResult

    :ivar start: The starting position of the D-segment.
    :vartype start: int
    :ivar end: The ending position of the D-segment.
    :vartype end: int
    :ivar seq: The sequence of the D-segment.
    :vartype seq: str
    :ivar notes: List of notes related to the D-segment.
    :vartype notes: list
    :ivar functionality: The functionality of the D-segment.
    :vartype functionality: str
    :ivar d_5_nonamer: The 5' nonamer sequence of the D-segment.
    :vartype d_5_nonamer: str
    :ivar d_5_heptamer: The 5' heptamer sequence of the D-segment.
    :vartype d_5_heptamer: str
    :ivar d_3_nonamer: The 3' nonamer sequence of the D-segment.
    :vartype d_3_nonamer: str
    :ivar d_3_heptamer: The 3' heptamer sequence of the D-segment.
    :vartype d_3_heptamer: str
    :ivar d_5_rss_start: The starting position of the 5' RSS of the D-segment.
    :vartype d_5_rss_start: int
    :ivar d_5_rss_end: The ending position of the 5' RSS of the D-segment.
    :vartype d_5_rss_end: int
    :ivar d_3_rss_start: The starting position of the 3' RSS of the D-segment.
    :vartype d_3_rss_start: int
    :ivar d_3_rss_end: The ending position of the 3' RSS of the D-segment.
    :vartype d_3_rss_end: int
    :ivar likelihood: The likelihood of the D-segment.
    :vartype likelihood: float

    :method annotate: Annotate the D-segment by determining its functionality.

    .. note::
        The class represents an annotation for a D-segment in the assembly sequence.
        It provides information about the D-segment's start and end positions, sequence,
        notes, functionality, motif sequences, RSS positions, and likelihood.
    """
    def __init__(self, assembly, start, end, left_motif, right_motif):
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
        """
        Annotate the D-segment by determining its functionality.

        :return: None
        """
        self.functionality = 'Functional'
        for note in self.notes:
            note = note.lower()
            if 'not found' in note or 'conserved' in note:
                self.functionality = 'ORF'
                break


def process_d(assembly, assembly_rc, germlines, conserved_motif_seqs, motifs, start, end, best, matches, D_5_RSS_SPACING, D_3_RSS_SPACING):
    """
    Process D-segment annotations.

    :param assembly: The assembly sequence.
    :type assembly: str
    :param assembly_rc: The reverse complement of the assembly sequence.
    :type assembly_rc: str
    :param germlines: Dictionary containing the germline sequences.
    :type germlines: dict
    :param conserved_motif_seqs: Dictionary containing conserved motif sequences.
    :type conserved_motif_seqs: dict
    :param motifs: Dictionary containing motifs.
    :type motifs: dict
    :param start: The starting position of the D-segment search.
    :type start: int
    :param end: The ending position of the D-segment search.
    :type end: int
    :param best: The best match information.
    :type best: dict
    :param matches: List of matches.
    :type matches: list

    :return: List of processed D-segment annotations.
    :rtype: list of dict

    .. note::
        The function processes D-segment annotations based on the provided parameters.
        It returns a list of dictionaries representing the processed D-segment annotations.

    .. seealso::
        - `find_compound_motif`: Function for finding compound motifs.
        - `calc_best_match_score`: Function for calculating the best match score.
        - `check_seq`: Function for checking a sequence.
   """
    #if start == 1678577:
    #   breakpoint()

    left_rss = find_compound_motif(assembly, conserved_motif_seqs, motifs["5'D-NONAMER"], motifs["5'D-HEPTAMER"], D_5_RSS_SPACING, D_5_RSS_SPACING, 15, end=start-1)
    right_rss = find_compound_motif(assembly, conserved_motif_seqs, motifs["3'D-HEPTAMER"], motifs["3'D-NONAMER"], D_3_RSS_SPACING, D_3_RSS_SPACING, 15, start=end)

    results = []

    for left in left_rss:
        for right in right_rss:
            annot = DAnnotation(assembly, left.end + 1, right.start - 1, left, right)
            if annot.likelihood:
                results.append(annot)    # explicitly don't add Ds with missing RSS components, there are too many

    rows = []
    for result in results:
        result.annotate()
        best_match, best_score, best_nt_diffs = calc_best_match_score(germlines, best, result.seq)

        row = {
            'start': result.start,
            'end': result.end,
            'end_rev': len(assembly) - result.start + 1,
            'start_rev': len(assembly) - result.end + 1,
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
        row['3_rss_start_rev'], row['3_rss_end_rev'] = len(assembly) - row['3_rss_end'] + 1, len(assembly) - row['3_rss_start'] + 1

        row['5_rss_start'] = result.d_5_rss_start
        row['5_rss_end'] = result.d_5_rss_end
        row['5_rss_start_rev'], row['5_rss_end_rev'] = len(assembly) - row['5_rss_end'] + 1, len(assembly) - row['5_rss_start'] + 1

        add_gene_coords(row, assembly)
        rows.append(row)

    for row in rows:
        check_seq(assembly, assembly_rc, row['seq'], row['start'], row['end'], 'seq', False)
        check_seq(assembly, assembly_rc, row['seq'], row['start_rev'], row['end_rev'], 'seq', True)
        check_seq(assembly, assembly_rc, row['d_5_nonamer'], row['5_rss_start'], row['5_rss_start'] + 8, 'd_5_nonamer', False)
        check_seq(assembly, assembly_rc, row['d_5_heptamer'], row['5_rss_end'] - 6, row['5_rss_end'], 'd_5_heptamer', False)
        check_seq(assembly, assembly_rc, assembly[row['5_rss_start'] - 1:row['5_rss_end']], row['5_rss_start_rev'], row['5_rss_end_rev'], '5_rss', True)
        check_seq(assembly, assembly_rc, row['d_3_heptamer'], row['3_rss_start'], row['3_rss_start'] + 6, 'd_3_heptamer', False)
        check_seq(assembly, assembly_rc, row['d_3_nonamer'], row['3_rss_end'], row['3_rss_end'] - 8, 'd_3_nonamer', False)
        check_seq(assembly, assembly_rc, assembly[row['3_rss_start'] - 1:row['3_rss_end']], row['3_rss_start_rev'], row['3_rss_end_rev'], '3_rss', True)

    return rows

class CAnnotation:
    """
    Represents an annotation for a C-segment.

    :param assembly: The assembly sequence.
    :type assembly: str
    :param start: The starting position of the C-segment.
    :type start: int
    :param end: The ending position of the C-segment.
    :type end: int

    :ivar start: The starting position of the C-segment.
    :vartype start: int
    :ivar end: The ending position of the C-segment.
    :vartype end: int
    :ivar seq: The sequence of the C-segment.
    :vartype seq: str
    :ivar notes: List of notes related to the C-segment.
    :vartype notes: list
    :ivar functionality: The functionality of the C-segment.
    :vartype functionality: str

    :method annotate: Annotate the C-segment by determining its functionality.

    .. note::
        The class represents an annotation for a C-segment in the assembly sequence.
        It provides information about the C-segment's start and end positions, sequence,
        notes, and functionality.
    """
    def __init__(self, assembly, start, end):
        self.start = start
        self.end = end
        self.seq = assembly[start-1:end]
        self.notes = []
        self.functionality = None

    """
    Annotate the C-segment by determining its functionality.

    :return: None
    """
    def annotate(self):
        self.functionality = 'Functional'
        self.notes = []


def process_c(assembly, assembly_rc, germlines, start, end, best, matches):
    """
    Process C-segment annotations.

    :param assembly: The assembly sequence.
    :type assembly: str
    :param assembly_rc: The reverse complement of the assembly sequence.
    :type assembly_rc: str
    :param germlines: The germline sequences.
    :type germlines: dict
    :param start: The starting position of the C-segment.
    :type start: int
    :param end: The ending position of the C-segment.
    :type end: int
    :param best: The best match information.
    :type best: dict
    :param matches: List of matches.
    :type matches: list

    :return: List of processed C-segment annotations.
    :rtype: list of dict

    .. note::
        The function processes C-segment annotations based on the provided parameters.
        It returns a list of dictionaries representing the processed C-segment annotations.

    .. seealso::
        - `calc_best_match_score`: Function for calculating the best match score.
        - `check_seq`: Function for checking a sequence.
    """
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

    best_match, best_score, best_nt_diffs = calc_best_match_score(germlines, best, result_seq)

    if best_score < 0.1: # just too many Ns
        return []

    row = {
        'start': result_start,
        'end': result_end,
        'end_rev': len(assembly) - result_start + 1,
        'start_rev': len(assembly) - result_end + 1,
        'matches': len(matches),
        'blast_match': best_match,
        'blast_score': best_score,
        'blast_nt_diffs': best_nt_diffs,
        'functional': 'Functional',
        'notes': '',
        'seq': result_seq,
    }

    add_gene_coords(row, assembly)

    check_seq(assembly, assembly_rc, row['seq'], row['start'], row['end'], 'seq', False)
    check_seq(assembly, assembly_rc, row['seq'], row['start_rev'], row['end_rev'], 'seq', True)

    return[row]

class JAnnotation:
    """
    Represents an annotation for a J-segment.

    :param assembly: The assembly sequence.
    :type assembly: str
    :param start: The starting position of the J-segment.
    :type start: int
    :param end: The ending position of the J-segment.
    :type end: int
    :param motif: The motif representing the J-segment.
    :type motif: MotifResult

    :ivar start: The starting position of the J-segment.
    :vartype start: int
    :ivar end: The ending position of the J-segment.
    :vartype end: int
    :ivar aa: The translated amino acid sequence of the J-segment.
    :vartype aa: str
    :ivar seq: The sequence of the J-segment.
    :vartype seq: str
    :ivar assembly: The assembly sequence.
    :vartype assembly: str
    :ivar notes: List of notes related to the J-segment.
    :vartype notes: list
    :ivar functionality: The functionality of the J-segment.
    :vartype functionality: str
    :ivar j_frame: The reading frame of the J-segment.
    :vartype j_frame: int
    :ivar nonamer: The nonamer sequence of the J-segment.
    :vartype nonamer: str
    :ivar heptamer: The heptamer sequence of the J-segment.
    :vartype heptamer: str
    :ivar rss_start: The start position of the RSS associated with the J-segment.
    :vartype rss_start: int
    :ivar rss_end: The end position of the RSS associated with the J-segment.
    :vartype rss_end: int
    :ivar likelihood: The likelihood of the J-segment.
    :vartype likelihood: float

    :method annotate: Annotate the J-segment by determining its functionality.

    .. note::
        The class represents an annotation for a J-segment in the assembly sequence.
        It provides information about the J-segment's start and end positions, amino acid
        sequence, sequence, assembly sequence, notes, functionality, reading frame, nonamer,
        heptamer, RSS positions, and likelihood.
    """
    def __init__(self, assembly, start, end, motif):
        self.start = start
        self.end = end
        self.aa = ''
        self.seq = assembly[start-1:end]
        self.assembly = assembly
        self.notes = motif.notes
        self.functionality = None
        self.j_frame = None
        self.nonamer = motif.left
        self.heptamer = motif.right
        self.rss_start = motif.start
        self.rss_end = motif.end
        self.likelihood = motif.likelihood

    def annotate(self, assembly, J_TRP_MOTIF, J_TRP_OFFSET, J_SPLICE):
        """
        Annotate the J-segment by determining its functionality.

        :param J_TRP_MOTIF: The J-TRP motif to search for.
        :type J_TRP_MOTIF: str
        :param J_TRP_OFFSET: The number of codons between the start of the J-TRP motif and the end of the J-REGION.
        :type J_TRP_OFFSET: int
        :param J_SPLICE: The J-splice sequence to search for.
        :type J_SPLICE: str

        :return: None
        """
        # consider each frame
        hit = 0
        for i in range(0, 3):
            j_codons = self.seq[i:]
            self.aa = simple.translate(j_codons)

            if not isinstance(J_TRP_MOTIF, list):
                J_TRP_MOTIF = [J_TRP_MOTIF]

            for motif in J_TRP_MOTIF:
                hit = find_best_match(self.aa, motif, thresh=1.0)
                if hit >= 0:
                    break
            
            if hit >= 0:
                break

        def check_splice(assembly, end, splice):
            return (splice[0] != '*' and assembly[end - 1:end + 2] == splice) or (splice[0] == '*' and assembly[end:end + 2] == splice[1:])
        
        if hit >= 0:
            self.j_frame = i
            self.end = min(self.start + i + (hit + J_TRP_OFFSET)*3, len(assembly))
            self.seq = self.assembly[self.start - 1:self.end]
            self.aa = simple.translate(self.assembly[self.start + i - 1:self.end])
            splice_found = False

            # This provides, if necesary, scope to search for a donor splice at a non-caonical position
            # Not used at present because it is not clear whether to use it, and, if so, how far to search
            for p in range(-3, 30, 3):
                if self.end + p < len(assembly) and check_splice(assembly, self.end + p, J_SPLICE):
                    self.end += p
                    self.seq = self.assembly[self.start - 1:self.end]
                    splice_found = True
                    break
            
            if not splice_found:
                self.notes.append('Donor splice not found')
                self.functionality = 'pseudo'
            else:
                self.functionality = 'Functional'
                if p != 0:
                    self.notes.append(f'Donor splice found at non-canonical position')
                    self.functionality = 'ORF'
                for note in self.notes:
                    if 'conserved' in note:
                        self.functionality = 'ORF'
                        break

            if self.functionality != 'pseudo':
                if '*' in self.aa[hit:]:
                    self.notes.append('Stop codon after J-TRP')
                    self.functionality = 'pseudo'
        else:
            self.notes.append('J-TRP not found')
            self.j_frame = 0
            self.end = min(self.start + (hit + J_TRP_OFFSET)*3, len(assembly))
            self.seq = self.assembly[self.start - 1:self.end]
            self.functionality = 'pseudo'

def process_j(assembly, assembly_rc, germlines, conserved_motif_seqs, motifs, start, end, best, matches, J_TRP_MOTIF, J_TRP_OFFSET, J_SPLICE, J_RSS_SPACING):
    """
    Process J-segment annotations.

    :param assembly: The assembly sequence.
    :type assembly: str
    :param assembly_rc: The reverse complement of the assembly sequence.
    :type assembly_rc: str
    :param germlines: The germline sequences.
    :type germlines: dict
    :param conserved_motif_seqs: Dictionary containing conserved motif sequences.
    :type conserved_motif_seqs: dict
    :param motifs: Dictionary containing motifs.
    :type motifs: dict
    :param start: The starting position of the J-segment search.
    :type start: int
    :param end: The ending position of the J-segment search.
    :type end: int
    :param best: The best match information.
    :type best: dict
    :param matches: List of matches.
    :type matches: list

    :return: List of processed J-segment annotations.
    :rtype: list of dict

    .. note::
        The function processes J-segment annotations based on the provided parameters.
        It returns a list of dictionaries representing the processed J-segment annotations.

    .. seealso::
        - `find_compound_motif`: Function for finding compound motifs.
        - `calc_best_match_score`: Function for calculating the best match score.
        - `check_seq`: Function for checking a sequence.

    """
    j_rss = find_compound_motif(assembly, conserved_motif_seqs, motifs['J-NONAMER'], motifs['J-HEPTAMER'], J_RSS_SPACING-1, J_RSS_SPACING, 15, end=start-1)

    if len(j_rss) == 0:
        return []


    # we need to make a call between the rss on the basis of likelihood
    # otherwise we'll get overlapping rss and call js that only differ in a bp or two at the 5' end

    j_rss.sort(key=attrgetter('likelihood'), reverse=True)
    j_rs = j_rss[0]

    result = (JAnnotation(assembly, j_rs.end + 1, end, j_rs))
    result.annotate(assembly, J_TRP_MOTIF, J_TRP_OFFSET, J_SPLICE)

    if result.functionality == 'Functional' and 'not found' in ', '.join(j_rs.notes):
        result.functionality = 'ORF'

    best_match, best_score, best_nt_diffs = calc_best_match_score(germlines, best, result.seq)

    row = {
        'start': result.start,
        'end': result.end,
        'end_rev': len(assembly) - result.start + 1,
        'start_rev': len(assembly) - result.end + 1,
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
    row['5_rss_start_rev'], row['5_rss_end_rev'] = len(assembly) - row['5_rss_end'] + 1, len(assembly) - row['5_rss_start'] + 1

    add_gene_coords(row, assembly)

    check_seq(assembly, assembly_rc, row['seq'], row['start'], row['end'], 'seq', False)
    check_seq(assembly, assembly_rc, row['seq'], row['start_rev'], row['end_rev'], 'seq', True)
    check_seq(assembly, assembly_rc, row['j_nonamer'], row['5_rss_start'], row['5_rss_start'] + 8, 'j_nonamer', False)
    check_seq(assembly, assembly_rc, row['j_heptamer'], row['5_rss_end'] - 6, row['5_rss_end'], 'j_heptamer', False)
    check_seq(assembly, assembly_rc, assembly[row['5_rss_start'] - 1:row['5_rss_end']], row['5_rss_start_rev'], row['5_rss_end_rev'], '5_rss', True)

    return [row]


class VAnnotation:
    """
    Represents a V-segment annotation.

    :param assembly: The assembly sequence.
    :type assembly: str
    :param start: The starting position of the V-segment annotation (1-based coordinate).
    :type start: int
    :param end: The ending position of the V-segment annotation (1-based coordinate).
    :type end: int
    :param refs: Optional references for sequence alignment.
    :type refs: list, optional

    :ivar start: The starting position of the V-segment annotation.
    :vartype start: int
    :ivar end: The ending position of the V-segment annotation.
    :vartype end: int
    :ivar ungapped: The ungapped sequence of the V-segment annotation.
    :vartype ungapped: str
    :ivar gapped: The gapped sequence of the V-segment annotation.
    :vartype gapped: str
    :ivar gapped_aa: The gapped amino acid sequence of the V-segment annotation.
    :vartype gapped_aa: str
    :ivar notes: List of additional notes related to the V-segment annotation.
    :vartype notes: list
    :ivar functionality: The functionality of the V-segment annotation ('Functional', 'ORF', 'pseudo').
    :vartype functionality: str
    :ivar likelihood: The likelihood value of the V-segment annotation.
    :vartype likelihood: None or float
    :ivar align_refs: Indicates if references are used for sequence alignment.
    :vartype align_refs: bool
    :ivar v_gapped_ref: The gapped reference sequence for alignment.
    :vartype v_gapped_ref: str or None
    :ivar v_ungapped_ref: The ungapped reference sequence for alignment.
    :vartype v_ungapped_ref: str or None

    .. note::
        If `refs` parameter is provided, the references for sequence alignment will be used.
        Otherwise, the alignment process will be skipped.

    """
    def __init__(self, assembly, start, end, refs=None):
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
        """
        Annotate the V-segment.

        This method performs sequence alignment and annotation of the V-segment.
        It updates the `gapped`, `gapped_aa`, `functionality`, and `notes` attributes.
        """
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
def process_v(assembly, assembly_rc, germlines, v_gapped_ref, v_ungapped_ref, conserved_motif_seqs, motifs, start, end, best, matches, align, V_RSS_SPACING, v_parsing_errors):
    """
     Process V-segment annotations.

     :param assembly: The assembly sequence.
     :type assembly: str
     :param assembly_rc: The reverse complement of the assembly sequence.
     :type assembly_rc: str
     :param germlines: The germline sequences.
     :type germlines: dict
     :param v_gapped_ref: The gapped reference sequence for V-segment alignment.
     :type v_gapped_ref: str
     :param v_ungapped_ref: The ungapped reference sequence for V-segment alignment.
     :type v_ungapped_ref: str
     :param conserved_motif_seqs: Dictionary containing conserved motif sequences.
     :type conserved_motif_seqs: dict
     :param motifs: Dictionary containing motifs.
     :type motifs: dict
     :param start: The starting position of the V-segment search (1-based).
     :type start: int
     :param end: The ending position of the V-segment search (1-based).
     :type end: int
     :param best: The best match information.
     :type best: dict
     :param matches: List of matches.
     :type matches: list
     :param align: Flag indicating whether sequence alignment should be performed.
     :type align: bool
     :param V_RSS_SPACING: The spacing between V-segment RSS motifs.
     :type V_RSS_SPACING: int
     :param v_parsing_errors: Dictionary for storing V-segment parsing errors.
     :type v_parsing_errors: dict

     :return: List of processed V-segment annotations.
     :rtype: list of dict

     .. note::
         The function processes V-segment annotations based on the provided parameters.
         It returns a list of dictionaries representing the processed V-segment annotations.

     .. seealso::
         - `find_compound_motif`: Function for finding compound motifs.
         - `VAnnotation`: Class representing V-segment annotations.
         - `calc_best_match_score`: Function for calculating the best match score.
         - `check_seq`: Function for checking a sequence.
   """
    leaders = find_compound_motif(assembly, conserved_motif_seqs, motifs['L-PART1'], motifs['L-PART2'], 10, 800, 8, end=start-1, right_force=start-len(motifs['L-PART2'].consensus))

    # restrain length to between 270 and 320 nt, allow a window anywhere within that range
    rights = find_compound_motif(assembly, conserved_motif_seqs, motifs['V-HEPTAMER'], motifs['V-NONAMER'], V_RSS_SPACING-1, V_RSS_SPACING, 25, start=start+295)

    # put in a dummy leader if we found an RSS but no leader

    if rights and not leaders:
        leaders.append(
            MotifResult(assembly,
                   SingleMotifResult(assembly, conserved_motif_seqs, motifs['L-PART1'], start-len(motifs['L-PART1'].consensus)-10-len(motifs['L-PART2'].consensus), 0.1),
                   SingleMotifResult(assembly, conserved_motif_seqs, motifs['L-PART2'], start-len(motifs['L-PART2'].consensus), 0.1),
                   ['Leader not found'])
        )

    # put in a dummy RSS if we found leader but no RSS

    if leaders and not rights:
        rights.append(
            MotifResult(assembly,
                   SingleMotifResult(assembly, conserved_motif_seqs, motifs['V-HEPTAMER'], min(end+1, len(assembly)), 0),
                   SingleMotifResult(assembly, conserved_motif_seqs, motifs['V-NONAMER'], min(end+V_RSS_SPACING+1, len(assembly)), 0),
                   ['RSS not found'])
        )

    best_rights = find_best_rss(rights)
    best_leaders = find_best_leaders(assembly, leaders)

    max_likelihood_rec = None
    max_likelihood_rec_at_start = None
    results = []

    # report one or more functional V-GENEs found between the identified leaders and rss.
    # if none, report the non-functional sequence with highest combined (leader, rss) likelihood

    if len(best_leaders) == 0:
        err_seq = assembly[max(start-500, 0):start]
        v_parsing_errors[start] = ((start, end, best['subject'], 'leader not found', err_seq))

    if len(best_rights) == 0:
        if len(assembly) >= end+24:
            err_seq = assembly[end+1:end+24]
        else:
            err_seq = assembly[end+1:]
        v_parsing_errors[start] = ((start, end, best['subject'], 'rss not found', err_seq))

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
            if align:
                v_annot = VAnnotation(assembly, left.end + 1, right.start - 1, refs=[v_gapped_ref, v_ungapped_ref])
            else:
                v_annot = VAnnotation(assembly, left.end + 1, right.start - 1)
            v_annot.annotate()
            v_annot.likelihood = left.likelihood * right.likelihood

            if v_annot.functionality == 'Functional':
                results.append((left, v_annot, right))
            if max_likelihood_rec is None or v_annot.likelihood > max_likelihood_rec[1].likelihood:
                max_likelihood_rec = (left, v_annot, right)
            if v_annot.start == start and (max_likelihood_rec_at_start is None or v_annot.likelihood > max_likelihood_rec[1].likelihood):
                max_likelihood_rec_at_start = (left, v_annot, right)

    # If we don't have any functional results, favour the one with the suggested start position.
    if len(results) == 0 and max_likelihood_rec_at_start is not None:
        results = [max_likelihood_rec_at_start]

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
                if 'stop codon' in note or 'atg' in note or 'donor-splice' in note or 'acceptor-splice' in note:
                    v_gene.functionality = 'pseudo'

        best_match, best_score, best_nt_diffs = calc_best_match_score(germlines, best, v_gene.ungapped)

        #if v_gene.start == 556809:
        #    breakpoint()

        row = {
            'start': v_gene.start,
            'end': v_gene.end,
            'end_rev': len(assembly) - v_gene.start + 1,
            'start_rev': len(assembly) - v_gene.end + 1,
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
        row['3_rss_start_rev'] = len(assembly) - rss.end + 1
        row['3_rss_end'] = rss.end
        row['3_rss_end_rev'] = len(assembly) - rss.start + 1

        row['l_part1_start'] = leader.start
        row['l_part1_end'] = leader.start + len(leader.left) - 1
        row['l_part1_start_rev'], row['l_part1_end_rev'] = len(assembly) - row['l_part1_end'] + 1, len(assembly) - row['l_part1_start'] + 1

        row['l_part2_start'] = leader.end - len(leader.right) + 1
        row['l_part2_end'] = leader.end
        row['l_part2_start_rev'], row['l_part2_end_rev'] = len(assembly) - row['l_part2_end'] + 1, len(assembly) - row['l_part2_start'] + 1

        row['aa'] = simple.translate(v_gene.ungapped)

        add_gene_coords(row, assembly)

        rows.append(row)

    for row in rows:
        check_seq(assembly, assembly_rc, row['seq'], row['start'], row['end'], 'seq', False)
        check_seq(assembly, assembly_rc, row['seq'], row['start_rev'], row['end_rev'], 'seq', True)
        check_seq(assembly, assembly_rc, row['v_heptamer'], row['3_rss_start'], row['3_rss_start'] + 6, 'v_heptamer', False)
        check_seq(assembly, assembly_rc, row['v_nonamer'], row['3_rss_end'], row['3_rss_end'] - 8, 'v_nonamer', False)
        check_seq(assembly, assembly_rc, assembly[row['3_rss_start'] - 1:row['3_rss_end']], row['3_rss_start_rev'], row['3_rss_end_rev'], '3_rss', True)
        check_seq(assembly, assembly_rc, row['l_part1'], row['l_part1_start'], row['l_part1_end'], 'l_part1', False)
        check_seq(assembly, assembly_rc, assembly[row['l_part1_start'] - 1:row['l_part1_end']], row['l_part1_start_rev'], row['l_part1_end_rev'], 'l_part1', True)
        check_seq(assembly, assembly_rc, row['l_part2'], row['l_part2_start'], row['l_part2_end'], 'l_part2', False)
        check_seq(assembly, assembly_rc, assembly[row['l_part2_start'] - 1:row['l_part2_end']], row['l_part2_start_rev'], row['l_part2_end_rev'], 'l_part1', True)
        remove_notfound_coords(row)

    return rows


# check co-ordinates of each element for sanity
# all co-ordinates should be 1-based!

def check_seq(assembly, assembly_rc, seq, start, end, name, rev):
    if len(seq) == 0:
        return  #   allow empty sequences which may occur if the assembly is truncated

    if end < start:
        end, start = start, end

    if not rev:
        if seq != assembly[start-1:end]:
            print(f'Error: sequence {name} ({start}, {end}) failed co-ordinate check.')
    else:
        if simple.reverse_complement(seq) != assembly_rc[start-1:end]:
            print(f'Error: reverse sequence {name} ({start}, {end}) failed co-ordinate check.')



# find the best leader PART1 for each PART2 starting position:
# if PART2 is in-frame: the in-frame PART1 giving best likelihood
# otherwise, the PART1 giving best likelihood regardless of frame
def find_best_leaders(assembly, leaders):
    def spread(n, inc=1):
        ret = []
        for i in range(0, n, inc):
            ret.append(i)
            if i != 0:
                ret.append(0 - i)
        return ret

    # return the needed change in the size of the right, to keep in frame with a change in the size of the left
    # i.e. change_in_right + change_in_left = [some multiple of 3 up to a limit]
    def change_in_right(change_in_left):
        ret = []
        for i in spread(12, 3):
            c = i - (change_in_left % 3)
            if c > -7:
                ret.append(c)
        return list(ret)

    leader_choices = defaultdict(list)
    for leader in leaders:
        leader_choices[leader.end].append(leader)

    best_leaders = {}
    for position, choices in leader_choices.items():
        good_leaders = []
        bad_leaders = []

        for choice in choices:
            if choice.left:
                bad_start = choice.left[:3] != 'ATG'

                #if position == 742228 and choice.start == 742085:
                #    breakpoint()

                for i in spread(16):
                    donor = assembly[choice.start - 1 + len(choice.left) + i:choice.start - 1 + len(choice.left) + 2 + i]
                    bad_donor = donor != 'GT'
                    bad_acceptor = None

                    if not bad_donor:
                        # a negative j should make the right smaller, i.e. the start co-ordinate bigger
                        for j in change_in_right(len(choice.left) + len(choice.right) + i):
                            acceptor = assembly[choice.end - len(choice.right) - j - 2:choice.end - len(choice.right) - j]
                            bad_acceptor = acceptor != 'AG'
                            if not bad_acceptor:
                                cl = assembly[choice.start - 1:choice.start - 1 + len(choice.left) + i]
                                cr = assembly[choice.end - len(choice.right) - j: choice.end]
                                l12p = simple.translate(cl + cr)

                                if not ('X' in l12p or '*' in l12p):
                                    choice.right = cr
                                    choice.left = cl
                                    break

                    if not bad_acceptor and not bad_donor:     # we found a solution
                        break
            else:
                bad_donor = True
                bad_start = True
                donor = ''

            if bad_donor:    # we didn't find a solution: just annotate the right as it stands, knowing the left is bad
                acceptor = assembly[choice.end - len(choice.right):choice.end - len(choice.right) + 2]
                bad_acceptor = acceptor != 'AG'

            l12p = simple.translate(choice.left + choice.right)
            stop_codon = 'X' in l12p or '*' in l12p

            if bad_start and choice.left:
                choice.notes.append('Leader missing initial ATG')
            if bad_donor and choice.left:
                choice.notes.append(f'Bad DONOR-SPLICE: {donor}')
            if bad_acceptor:
                choice.notes.append(f'Bad ACCEPTOR-SPLICE: {acceptor}')
            if stop_codon:
                choice.notes.append('Stop codon in leader')

            if bad_start or bad_donor or bad_acceptor or stop_codon:
                bad_leaders.append(choice)
            else:
                good_leaders.append(choice)

        if good_leaders:
            best_leaders[position] = sorted(good_leaders, key=attrgetter('likelihood'), reverse=True)[0]
        else:
            bad_leaders.sort(key=attrgetter('likelihood'), reverse=True)
            for leader in bad_leaders:
                l2p = simple.translate(leader.right[2:])
                if not ('X' in l12p or '*' in l12p):
                    best_leaders[position] = leader
                    break
            if position not in best_leaders:
                best_leaders[position] = bad_leaders[0]

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


aligner_global = PairwiseAligner(
    mode = 'global',
    open_gap_score = -1,
    extend_gap_score = -1,
    match_score = 1,
    mismatch_score = 0
)


def calc_best_match_score(germlines, best, seq):
    if not seq or not best['subject']:
        return '', 0.0, 0

    score = aligner_global.align(seq, germlines[best['subject']]).score
    
    if isinstance(score, float) and score > 0:
        best_score = round(100 * score / len(seq), 2)
        best_match = best['subject']
    else:
        best_score = 0.0
        best_match = ''

    best_nt_diffs = simple.nt_diff(seq, germlines[best['subject']])

    return best_match, best_score, best_nt_diffs


def find_best_match(seq, pattern, thresh=0.7):
    dists = find_all_matches(seq, pattern, thresh)

    if len(dists) == 0:
        return -1

    scores = [d[1] for d in dists]

    return dists[scores.index(min(scores))][0]

def matches(s1, s2):
    return sum([1 for i in range(len(s1)) if s1[i] == s2[i] and s2[i] != 'X'])

def find_all_matches(seq, pattern, thresh=0.7):
    dists = []
    pattern_length = len(pattern)

    significant_bps = len([x for x in pattern if x != 'X'])
    min_result = ceil(thresh * significant_bps)

    for i in range(len(seq) - len(pattern)):
        d = matches(seq[i:i+pattern_length], pattern)
        if d >= min_result:
            dists.append((i, d))

    return dists


# Find all 'start' or 'end' coords in a record
def find_coords(rec, keyword):
    coords = []

    for k in rec.keys():
        if keyword in k and 'rev' not in k:
            coords.append(int(rec[k]))

    return coords


# Add overall gene coordinates to a row
def add_gene_coords(row, assembly):
    row['gene_start'] = max(min(find_coords(row, 'start')), 1)
    row['gene_end'] = min(max(find_coords(row, 'end')), len(assembly))
    row['gene_start_rev'], row['gene_end_rev'] = len(assembly) - row['gene_end'] + 1, len(assembly) - row['gene_start'] + 1
    row['gene_seq'] = assembly[row['gene_start'] - 1:row['gene_end']]


def reverse_coords(row):
    start_revs = [x for x in row.keys() if 'start_rev' in x]
    for start_rev in start_revs:
        end_rev = start_rev.replace('start', 'end')
        start = start_rev.replace('_rev', '')
        end = end_rev.replace('_rev', '')
        row[start], row[start_rev] = row[start_rev], row[start]
        row[end], row[end_rev] = row[end_rev], row[end]


# Remove any coordinates for sequences that were not found
def remove_notfound_coords(row):
    start_revs = [x for x in row.keys() if 'start_rev' in x]
    for start_rev in start_revs:
        seq_name = start_rev.replace('_start_rev', '')
        if seq_name in row and not row[seq_name]:
            del row[start_rev]
            del row[start_rev.replace('start_rev', 'end_rev')]
            del row[start_rev.replace('start_rev', 'start')]
            del row[start_rev.replace('start_rev', 'end')]

