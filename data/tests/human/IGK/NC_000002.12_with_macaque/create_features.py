# Build features file from manually-created comparison of light chains annotated by Watson et al. (2015i)

from receptor_utils import simple_bio_seq as simple
import argparse

parser = argparse.ArgumentParser(description='Compare digger results to an IMGT annotation')
parser.add_argument('digger_results', help='Digger results file')
parser.add_argument('watson_annotation', help='Annotation from Watson et al.')
parser.add_argument('coordinate_bias', help='Amount to subtract from Watson et al. co-ordinates to align with digger co-ordinates')
parser.add_argument('germline_ref', help='V(D)J germline reference set (Vs ungapped)')
parser.add_argument('--strand_filter', help='only process records with matching strand')
parser.add_argument('output_file', help='output_file')


args = parser.parse_args()

digger_results = simple.read_csv(args.digger_results)
ref_results = simple.read_csv(args.watson_annotation)
coord_bias = int(args.coordinate_bias)
germ_ref = simple.read_fasta(args.germline_ref)
strand_filter = args.strand_filter

results = []


for ref_result in ref_results:
    res = {'start': 0, 'end': 0, 'digger_coord_match': '', 'sense': '', 'allele': '', 'functional': '', 'v-region': '', 'd-region': '', 'j-region': ''}
    if not ref_result['Position']:
        continue

    if strand_filter and ref_result['strand'] != strand_filter:
        continue

    ref_start, ref_end = ref_result['Position'].split('-')
    ref_end = int(ref_end) - coord_bias
    ref_start = int(ref_start.split(':')[1]) - coord_bias

    germ_type = ''
    if 'V' in  ref_result['CH17 allele']:
        germ_type = 'v'
    elif 'D' in ref_result['CH17 allele']:
        germ_type = 'd'
    elif 'J' in  ref_result['CH17 allele']:
        germ_type = 'j'
    else:
        continue

    for dig_result in digger_results:
        dig_result['start'] = int(dig_result['start'])
        dig_result['end'] = int(dig_result['end'])

        if dig_result['end'] < dig_result['start']:
            (dig_result['start'], dig_result['end']) = (dig_result['end'], dig_result['start'])

        if abs(dig_result['start'] - ref_start) < 5000 and ref_result['CH17 allele'].replace('D', '') == dig_result['imgt_match'].replace('D', ''):
            res['start'] = dig_result['start']
            res['end'] = dig_result['end']
            res['digger_coord_match'] = 'T'
            break

    if not res['start']:
        res['start'] = ref_start
        res['end'] = ref_end
        res['digger_coord_match'] = 'F'

    res['sense'] = 'forward' if ref_result['strand'] == '+' else 'reverse'
    res['allele'] = ref_result['CH17 allele']
    res['functional'] = 'functional' if ref_result['CH17 Status'] == 'F' else ref_result['CH17 Status']
    res[f'{germ_type}-region'] = germ_ref[ref_result['CH17 allele']] if ref_result['CH17 allele'] in germ_ref else 'NNNNNNN'

    results.append(res)

simple.write_csv(args.output_file, results)





