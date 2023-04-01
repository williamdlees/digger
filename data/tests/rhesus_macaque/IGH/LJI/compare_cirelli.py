# Compare annotation published by Cirelli et al. (2019) with that produced by digger

import argparse
from receptor_utils import simple_bio_seq as simple


parser = argparse.ArgumentParser(description='Compare annotation published by Cirelli et al. (2019) with that produced by digger')
parser.add_argument('cirelli_file', help='supplemental table S1 from Cirelli et al. in csv format')
parser.add_argument('digger_file', help='digger annotation of a single locus')
parser.add_argument('locus', help='locus (IGH, IGK, IGL)')

args = parser.parse_args()

cirelli_recs = simple.read_csv(args.cirelli_file)
cirelli_recs = [r for r in cirelli_recs if args.locus in r['Sequence.ID'] and r['Annotation.Type'] == 'Primary']
cirelli_by_id = {r['Sequence.ID']: r for r in cirelli_recs}
cirelli_by_seq = {r['V.Region.Sequence']: r for r in cirelli_recs}

digger_recs = simple.read_csv(args.digger_file)
digger_by_seq = {}

for rec in digger_recs:
    if rec['seq'] not in digger_by_seq:
        digger_by_seq[rec['seq']] = []
    digger_by_seq[rec['seq']].append(rec)


# Cirelli annotated with digger status

res = []

for row in cirelli_recs:
    r = dict(row)
    dn = []

    if row['V.Region.Sequence'] in digger_by_seq:
        for hit in digger_by_seq[row['V.Region.Sequence']]:
            dn.append(f"{hit['contig']} {hit['start']} ({hit['functional']})")
    else:
        for seq in digger_by_seq.keys():
            if (row['V.Region.Sequence'] in seq or seq in row['V.Region.Sequence']) and abs(len(row['V.Region.Sequence']) - len(seq)) <= 25:
                for hit in digger_by_seq[seq]:
                    note = f"{hit['contig']} {hit['start']} ({hit['functional']}, length mismatch)"
                    if note not in dn:
                        dn.append(f"{hit['contig']} {hit['start']} ({hit['functional']}, length mismatch)")

    r['Digger'] = ', '.join(dn)
    res.append(r)

simple.write_csv('cirelli_table_1_with_digger_hits.csv', res)

# Digger F/ORF annotated with Cirelli status

res = []

for row in digger_recs:
    if row['functional'] in ('Functional', 'ORF'):
        r = dict(row)
        cn = []

        if row['seq'] in cirelli_by_seq:
            hit = cirelli_by_seq[row['seq']]
            cn.append(f"{hit['Sequence.ID']} {hit['Rh.Contig.Name']}")

        r['Cirelli'] = ', '.join(cn)
        res.append(r)

simple.write_csv('digger_F_ORF_with_cirelli_annots.csv', res)
