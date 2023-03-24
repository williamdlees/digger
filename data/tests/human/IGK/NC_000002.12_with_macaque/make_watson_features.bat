# these results need to be massaged into one file: watson_results.csv

cd work
python ../create_features.py --strand_filter "-" NC_000002.12_reverse.csv ../watson_supp_table_2.csv 89156874 Homo_sapiens_IGKVJ.fasta ../watson_features_rev.csv
# take the first two sequences from the following (these are before the interdup region)
python ../create_features.py --strand_filter "+" NC_000002.12.csv ../watson_supp_table_2.csv 89156874 Homo_sapiens_IGKVJ.fasta ../watson_features_fw1.csv
# take the other sequences from this run (these are after the interdup region, co-ord bias needs adjusting)
python ../create_features.py --strand_filter "+" NC_000002.12.csv ../watson_supp_table_2.csv 88896166 Homo_sapiens_IGKVJ.fasta ../watson_features_fw2.csv
