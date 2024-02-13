from receptor_utils import simple_bio_seq as simple

recs = simple.read_csv('IMGT000024_genes_lpart1.csv')

for rec in recs:
	rec['AA'] = simple.translate(rec['l-part1'])
	
simple.write_csv('translated.csv', recs)



