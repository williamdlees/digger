import csv

locus = 'IGH'
input_file = 'reference_sets/cirelli.csv'
output_file = 'reference_sets/cirelli_IGH.fasta'

with open(input_file, 'r') as fi, open(output_file, 'w') as fo:
    reader = csv.reader(fi, dialect='excel')
    for row in reader:
        if locus in row[0]:
            fo.write('>%s\n%s\n' % (row[0], row[1]))

