# Compare the alleles identified in one annotation - say IMGT's - with another 'reference'
# - which were found in the second annotation, which missed or identified differently, with respect to the first
import argparse
import csv
import os
from os.path import basename

from matplotlib_venn import venn2_unweighted
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf


def get_parser():
    parser = argparse.ArgumentParser(description='Compare digger results to an IMGT annotation')
    parser.add_argument('digger_results', help='Digger results file')
    parser.add_argument('annotation_file', help='IMGT annotation produced by parse_imgt_assembly_x.py')
    parser.add_argument('sense', help='Sense of annotation compared to digger results (forward or reverse)')
    parser.add_argument('outfile', help='Output file name (will create .csv, .jpg, .txt')
    parser.add_argument('-nc', help='include sequences for leader and rss', action='store_true')
    parser.add_argument('--filter_annot', help='filter IMGT annotations by sense (forward or reverse)')
    parser.add_argument('--comp_name', help='name to use for comparison (default IMGT)')
    parser.add_argument('--target_locus', help='Only consider IMGT matches to the target locus (used for TRA/TRD)')

    return parser


class DictIterator(dict):
    def __iter__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ind = -1
        self.sorted_keys = sorted(list(self.keys()))
        return self

    def __next__(self):
        if self.ind < len(self.sorted_keys)-1:
            self.ind += 1
            return self[self.sorted_keys[self.ind]]
        else:
            return None

def near(i, j, diff):
    return abs(i - j) < diff


def main():
    args = get_parser().parse_args()
    candidate = args.digger_results
    reference = args.annotation_file
    output = os.path.splitext(args.outfile)[0]

    if args.sense not in ['forward', 'reverse']:
        print("Error - sense must be 'forward' or 'reverse'")
        exit(0)

    if args.filter_annot and args.filter_annot not in ['forward', 'reverse']:
        print("Error - sense must be 'forward' or 'reverse'")
        exit(0)

    if args.sense == 'forward':
        c_start = 'start'
        c_end = 'end'
    else:
        c_start = 'start_rev'
        c_end = 'end_rev'

    results = create_results(c_end, c_start, candidate, output + '.csv', reference, args.nc, args.filter_annot, args.target_locus)
    plot_results(results, output + '.jpg', args.comp_name if args.comp_name else 'IMGT', args.target_locus)
    report_results(results, output + '.txt', args.comp_name if args.comp_name else 'IMGT', args.target_locus)


def plot_results(results, plotfile, comp_name, target_locus):
    fig = plt.figure(figsize=(10, 10))

    i = 0
    for chain in 'V', 'D', 'J':
        # eliminate false positives where digger identifies a match that IMGT assigns to another locus
        imgt_starts = set([r['imgt_start'] for r in results if r['type'] == chain and r['imgt_start'] and r['imgt_func'] in ['functional'] and (not target_locus or r['imgt_allele'][:3] == target_locus)])
        digger_starts = set([r['digger_start'] for r in results if r['type'] == chain and r['digger_start'] and r['digger_func'].lower() in ['functional'] and (not target_locus or r['imgt_allele'][:3] == target_locus)])

        # if there is a start in digger with no exact match in imgt, but there is a match within 10nt, adjust the two to be
        # identical - this way we count overlapping functional sequences with minor length differences as matches

        for digger_start in digger_starts:
            if digger_start not in imgt_starts:
                for imgt_start in list(imgt_starts):
                    if near(imgt_start, digger_start, 10):
                        imgt_starts.add(digger_start)
                        imgt_starts.remove(imgt_start)
                        break

        if len(imgt_starts) > 0 or len(digger_starts) > 0:
            fig.add_subplot(2, 2, i + 1)
            venn2_unweighted([imgt_starts, digger_starts], (comp_name, 'digger'))
            plt.title(f'{chain} - Functional genes identified')
            i += 1

        print(imgt_starts-digger_starts)
        print(digger_starts-imgt_starts)

    plt.savefig(plotfile)
    plt.close()
    return

def report_result_set(result_set, chain, fo):
    for result in result_set:
        start = result['imgt_start'] if result['imgt_start'] else result['digger_start']
        end = result['imgt_end'] if result['imgt_end'] else result['digger_end']

        fo.write(f"{start},{end}  imgt_match: {result['imgt_allele']}  imgt_func: {result['imgt_func']} digger_match: {result['digger_allele']}  digger_func: {result['digger_func'].lower()}  notes: {result['digger_notes']}\n")


def report_results(results, reportfile, comp_name, target_locus):
    with open(reportfile, 'w') as fo:
        for chain in 'V', 'D', 'J':
            fo.write(f"\n\ntype: {chain}\n")
            fo.write(f'Functional sequences reported by digger but not by {comp_name}:\n')
            # eliminate false positives where digger identifies a match that IMGT assigns to another locus
            report_result_set([r for r in results if r['digger_func'].lower() == 'functional' and r['imgt_func'].lower() != 'functional' and r['type'] == chain and (not target_locus or r['imgt_allele'][:3] == target_locus)], chain, fo)

            fo.write(f'\nFunctional sequences reported by {comp_name} but not by digger:\n')
            report_result_set([r for r in results \
                               if r['digger_func'].lower() != 'functional' 
                                and r['imgt_func'].lower() == 'functional' 
                                and r['type'] == chain
                                and (not target_locus or r['imgt_allele'][:3] == target_locus)
                              ], 
                              chain, fo
                            )

            fo.write(f'\nFunctional sequences reported by both {comp_name} and digger but with different sequences:\n')
            report_result_set([r for r in results if r['digger_func'].lower() == 'functional' and r['imgt_func'].lower() == 'functional' and r['seq_matches'] == 'N' and r['type'] == chain], chain, fo)


def create_results(c_end, c_start, candidate, output, reference, show_non_coding, filter_annot, target_locus):
    with open(reference, 'r') as ref_i, open(candidate, 'r') as cand_i, open(output, 'w', newline='') as fo:
        results = []

        if show_non_coding:
            headers = ['imgt_start', 'imgt_end', 'type', 'digger_start', 'digger_end', 'imgt_allele', 'digger_allele', 'imgt_func', 'digger_func', 'call_matches',
                       'func_matches', 'seq_matches', 'start_matches', 'digger_notes', 'imgt_l_part1', 'digger_l_part1', 'imgt_l_part2', 'digger_l_part2',
                       'imgt_heptamer', 'digger_heptamer', 'imgt_nonamer', 'digger_nonamer']
        else:
            headers = ['imgt_start', 'imgt_end', 'type', 'digger_start', 'digger_end', 'imgt_allele', 'digger_allele', 'imgt_func', 'digger_func', 'call_matches',
                       'func_matches', 'seq_matches', 'start_matches', 'digger_notes']

        result_writer = csv.DictWriter(fo, fieldnames=headers, extrasaction='ignore')
        result_writer.writeheader()

        ref_reader = csv.DictReader(ref_i)
        ref_alleles = {}
        for row in ref_reader:
            if filter_annot and row['sense'] != filter_annot:
                continue

            #if target_locus and target_locus not in row['allele']:
            #    continue

            try:
                row['start'] = int(row['start'])
            except:
                row['start'] = 0

            try:
                row['end'] = int(row['end'])
            except:
                row['end'] = 0

            ref_alleles[row['start']] = row
        ref_iter = iter(DictIterator(ref_alleles))

        cand_reader = csv.DictReader(cand_i)
        cand_alleles = {}
        for row in cand_reader:
            row[c_start] = int(row[c_start])
            row[c_end] = int(row[c_end])

            if row[c_start] > row[c_end]:
                (row[c_start], row[c_end]) = (row[c_end], row[c_start])
                notes = row['notes'].split(', ') if row['notes'] else []
                notes = ['Gene reversed'] + notes
                row['notes'] = ', '.join(notes)

            cand_alleles[int(row[c_start])] = row
        cand_iter = iter(DictIterator(cand_alleles))

        cand_row = next(cand_iter)
        ref_row = next(ref_iter)

        while True:
            call_matches = 'No'
            func_matches = 'No'
            seq_matches = 'No'
            start_matches = 'No'

            # if not ref_row:
            #    breakpoint()

            if cand_row is None and ref_row is None:
                break
            elif cand_row is not None and ref_row is not None and near(ref_row['start'], cand_row[c_start], 20):
                call_matches = 'Yes' if ref_row['allele'] == cand_row['imgt_match'] else 'No'
                func_matches = 'Yes' if ref_row['functional'].lower() == cand_row['functional'].lower() else 'No'

                chain = ''
                if cand_row['gene_type'][3] == 'V':
                    ref_region = 'v-region'
                    chain = 'V'
                elif cand_row['gene_type'][3] == 'D':
                    ref_region = 'd-region'
                    chain = 'D'
                elif cand_row['gene_type'][3] == 'J':
                    ref_region = 'j-region'
                    chain = 'J'

                ref_row[ref_region] = ref_row[ref_region].lower()
                cand_row['seq'] = cand_row['seq'].lower()
                if ref_row[ref_region] == cand_row['seq']:
                    seq_matches = 'Yes'
                elif ref_row[ref_region] in cand_row['seq'] or cand_row['seq'] in ref_row[ref_region]:
                    sense = 'shorter' if len(ref_row[ref_region]) > len(cand_row['seq']) else 'longer'
                    seq_matches = f"Digger is {sense} by {abs(len(ref_row[ref_region]) - len(cand_row['seq']))} nt"


                start_matches = 'Yes' if ref_row['start'] == cand_row[c_start] else 'No'

                res = {
                    'imgt_start': ref_row['start'],
                    'imgt_end': ref_row['end'],
                    'type': chain,
                    'digger_start': cand_row[c_start],
                    'digger_end': cand_row[c_end],
                    'imgt_allele': ref_row['allele'],
                    'digger_allele': cand_row['imgt_match'],
                    'imgt_func': ref_row['functional'],
                    'digger_func': cand_row['functional'],
                    'digger_notes': cand_row['notes'],
                    'call_matches': call_matches,
                    'func_matches': func_matches,
                    'seq_matches': seq_matches,
                    'start_matches': start_matches,
                }

                if show_non_coding:
                    res += {
                        'imgt_l_part1': ref_row['l-part1'],
                        'digger_l_part1': cand_row['l_part1'],
                        'imgt_l_part2': ref_row['l-part2'],
                        'digger_l_part2': cand_row['l_part2'],
                        'imgt_heptamer': ref_row['v-heptamer'],
                        'imgt_nonamer': ref_row['v-nonamer'],
                        'digger_heptamer': cand_row['v_heptamer'],
                        'digger_nonamer': cand_row['v_nonamer'],
                    }

                results.append(res)
                cand_row = next(cand_iter)
                ref_row = next(ref_iter)

            elif ref_row is None or (cand_row is not None and cand_row[c_start] < ref_row['start']):
                chain = ''
                if cand_row['gene_type'][3] == 'V':
                    chain = 'V'
                elif cand_row['gene_type'][3] == 'D':
                    chain = 'D'
                elif cand_row['gene_type'][3] == 'J':
                    chain = 'J'

                res = {
                    'imgt_start': '',
                    'imgt_end': '',
                    'type': chain,
                    'digger_start': cand_row[c_start],
                    'digger_end': cand_row[c_end],
                    'imgt_allele': '',
                    'digger_allele': cand_row['imgt_match'],
                    'imgt_func': '',
                    'digger_func': cand_row['functional'],
                    'digger_notes': cand_row['notes'],
                    'call_matches': call_matches,
                    'func_matches': func_matches,
                    'seq_matches': seq_matches,
                    'start_matches': start_matches,
                }

                results.append(res)
                cand_row = next(cand_iter)
            elif ref_row is not None:
                chain = ''
                if 'V' in ref_row['allele']:
                    chain = 'V'
                elif 'D' in ref_row['allele']:
                    chain = 'D'
                elif 'J' in ref_row['allele']:
                    chain = 'J'

                res = {
                    'imgt_start': ref_row['start'],
                    'imgt_end': ref_row['end'],
                    'type': chain,
                    'digger_start': '',
                    'digger_end': '',
                    'imgt_allele': ref_row['allele'],
                    'digger_allele': '',
                    'imgt_func': ref_row['functional'],
                    'digger_func': '',
                    'digger_notes': '',
                    'call_matches': call_matches,
                    'func_matches': func_matches,
                    'seq_matches': seq_matches,
                    'start_matches': start_matches,
                }
                results.append(res)
                ref_row = next(ref_iter)

        result_writer.writerows(results)
    return results

if __name__ == "__main__":
    main()

