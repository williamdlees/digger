# Search for single and double motifs
from operator import attrgetter


class SingleMotifResult:
    """
    Represents a single motif result.

    :param assembly: The assembly sequence.
    :type assembly: str
    :param conserved_motif_seqs: Dictionary containing conserved motif sequences.
    :type conserved_motif_seqs: dict
    :param motif: The motif object.
    :type motif: Motif
    :param position: The starting position of the motif (1-based).
    :type position: int
    :param likelihood: The likelihood value.
    :type likelihood: float

    :ivar start: The starting position of the motif  (1-based).
    :vartype start: int
    :ivar name: The name of the motif.
    :vartype name: str
    :ivar end: The ending position of the motif  (1-based).
    :vartype end: int
    :ivar seq: The sequence of the motif.
    :vartype seq: str
    :ivar likelihood: The likelihood value.
    :vartype likelihood: float
    :ivar notes: List of notes related to the motif.
    :vartype notes: list

    .. note::
        If `likelihood` is truthy, the `check_motif_consensus` method is called to check for motif consensus.

    .. method:: check_motif_consensus(conserved_motif_seqs)

        Check residues that are strongly conserved across all loci.

        :param conserved_motif_seqs: Dictionary containing conserved motif sequences.
        :type conserved_motif_seqs: dict
    """

    def __init__(self, assembly, conserved_motif_seqs, motif, position, likelihood):
        self.start = position
        self.name = motif.name
        self.end = position + len(motif.consensus) - 1
        self.seq = assembly[position - 1:position - 1 + len(motif.consensus)]
        self.likelihood = likelihood
        self.notes = []

        if likelihood:      # don't check dummy motifs
            self.check_motif_consensus(conserved_motif_seqs)


    # check residues that are strongly conserved across all loci
    def check_motif_consensus(self, conserved_motif_seqs):
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


class MotifResult:
    """
    Represents a motif result.

    :param assembly: The assembly sequence.
    :type assembly: str
    :param left_motif: The left motif object.
    :type left_motif: Motif
    :param right_motif: The right motif object.
    :type right_motif: Motif
    :param notes: Optional list of notes.
    :type notes: list, optional

    :ivar start: The starting position of the motif result.
    :vartype start: int
    :ivar end: The ending position of the motif result.
    :vartype end: int
    :ivar left: The sequence of the left motif.
    :vartype left: str
    :ivar gap: The sequence between the left and right motifs.
    :vartype gap: str
    :ivar right: The sequence of the right motif.
    :vartype right: str
    :ivar likelihood: The likelihood value of the motif result.
    :vartype likelihood: float
    :ivar notes: List of notes related to the motif result.
    :vartype notes: list

    .. note::
        If `notes` is provided, it will be used as the initial value for the `notes` attribute.
        Otherwise, an empty list will be used as the initial value.
    """
    def __init__(self, assembly, left_motif, right_motif, notes=None):
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


def find_single_motif(assembly, conserved_motif_seqs, motif, min_start, max_start):
    """
    Find a single motif within the specified range of start positions.

    :param assembly: The assembly sequence.
    :type assembly: str
    :param conserved_motif_seqs: Dictionary containing conserved motif sequences.
    :type conserved_motif_seqs: dict
    :param motif: The motif object.
    :type motif: Motif
    :param min_start: The minimum start position (1-based coordinate) for searching the motif.
    :type min_start: int
    :param max_start: The maximum start position (1-based coordinate) for searching the motif (inclusive).
    :type max_start: int

    :return: List of discovered motifs and their likelihoods.
    :rtype: list of SingleMotifResult

    .. note::
        The function searches for the motif within the specified range of start positions in the assembly sequence.
        It returns a list of `SingleMotifResult` objects representing the discovered motifs and their likelihoods.

    .. warning::
        The function assumes that the `Motif` class has a `calc_likelihood` method to calculate the likelihood of a given sequence.
    """
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
            res.append(SingleMotifResult(assembly, conserved_motif_seqs, motif, p, likelihood))

    return res




# Find compound motifs - ie heptamer plus nonamer, or l-part1 and l-part2.
# left, right - the two motifs
# min_gap, max_gap - the gap between them
# start co-ord, end co-ord: 1-based co-ordinates of either the desired ending position of left, or starting position of right (specify ome or the other)
# window - the function will look either for a left that ends at end +/- window, or a right that starts at start +/- window
# returns [{start, end, left_seq, gap_seq, right_seq, likelihood}] for each left, right combination that meets the threshold criterion
# If force has a value, always report a right starting at that position, with the best available left. This is used to force an L-PART2 to be discovered
# at the point that the BLAST alignment indicates the v-region should start
def find_compound_motif(assembly, conserved_motif_seqs, left_motif, right_motif, min_gap, max_gap, window, start=None, end=None, right_force=None):
    """
    Find compound motifs, such as heptamer plus nonamer or l-part1 and l-part2.

    :param assembly: The assembly sequence.
    :type assembly: str
    :param conserved_motif_seqs: Dictionary containing conserved motif sequences.
    :type conserved_motif_seqs: dict
    :param left_motif: The left motif object.
    :type left_motif: Motif
    :param right_motif: The right motif object.
    :type right_motif: Motif
    :param min_gap: The minimum gap length between the left and right motifs.
    :type min_gap: int
    :param max_gap: The maximum gap length between the left and right motifs.
    :type max_gap: int
    :param window: The window size used for motif search.
    :type window: int
    :param start: The starting position (1-based coordinate) for the left motif or None. (default: None)
    :type start: int, optional
    :param end: The ending position (1-based coordinate) for the right motif or None. (default: None)
    :type end: int, optional
    :param right_force: The forced starting position (1-based coordinate) for the right motif or None. (default: None)
    :type right_force: int, optional

    :return: List of compound motif results.
    :rtype: list of MotifResult

    .. note::
        The function searches for compound motifs within the specified parameters.
        It returns a list of `MotifResult` objects representing the discovered compound motifs.

    .. seealso::
        - `find_single_motif`: Function for finding a single motif within a range of start positions.
        - `MotifResult`: Class representing a motif result.
    """
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
            left_results = find_single_motif(assembly, conserved_motif_seqs, left_motif, start - window, start + window)
        else:
            left_results = find_single_motif(assembly, conserved_motif_seqs, left_motif, min_start, min_start + 2*window)

        right_results = []
        for left_result in left_results:
            right_results = find_single_motif(assembly, conserved_motif_seqs, right_motif, left_result.start + len(left_motif.consensus) + min_gap, left_result.start + len(left_result.seq) + max_gap)
            for right_result in right_results:
                res.append(MotifResult(assembly, left_result, right_result))

        # put in a dummy if we missed the lower likelihood motif
        if left_results and not right_results:
            best_left = sorted(left_results, key=attrgetter('likelihood'), reverse=True)[0]
            res.append(MotifResult(assembly, best_left, SingleMotifResult(assembly, conserved_motif_seqs, right_motif, best_left.start + len(best_left.seq) + min_gap, 0), notes=[f'{right_motif.name} not found']))
    else:
        if end is not None:
            right_results = find_single_motif(assembly, conserved_motif_seqs, right_motif, end - len(right_motif.consensus) - window, end + window)
        else:
            right_results = find_single_motif(assembly, conserved_motif_seqs, right_motif, min_end - len(right_motif.consensus), max_end - len(right_motif.consensus))

        # if right_force is set, if necessary, force discovery of a right motif at the forced position

        if right_force and right_force not in [r.start for r in right_results]:
            right_results.append(SingleMotifResult(assembly, conserved_motif_seqs, right_motif, right_force, 0))

        left_results = []
        for right_result in right_results:
            left_results = find_single_motif(assembly, conserved_motif_seqs, left_motif, right_result.start - max_gap - len(left_motif.consensus), right_result.start - min_gap - len(left_motif.consensus))
            for left_result in left_results:
                res.append(MotifResult(assembly, left_result, right_result))

        # put in a dummy if we missed the lower likelihood motif
        if right_results and not left_results:
            best_right = sorted(right_results, key=attrgetter('likelihood'), reverse=True)[0]
            res.append(MotifResult(assembly, SingleMotifResult(assembly, conserved_motif_seqs, left_motif, best_right.start - min_gap - len(left_motif.consensus), 0), best_right, notes=[f'{left_motif.name} not found']))

    return res


