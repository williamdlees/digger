
makeblastdb -in Macaca_mulatta_IGKV_fixed.fasta -dbtype nucl

blastn -db Macaca_mulatta_IGKV_fixed.fasta -query IMGT000063.fasta -out rhemac10_IGKV.out -outfmt 7 -gapopen 5 -gapextend 5 -penalty -1 -word_size 11

python ../../../../src/blastresults_to_csv.py rhemac10_IGKV.out rhemac10_IGKV_

python ../../../../src/find_alignments.py Macaca_mulatta_IGKV_fixed.fasta IMGT000063.fasta "rhemac10_IGKV_igk.csv" ../../../rhesus_motifs/IGK  -ref imgt,Macaca_mulatta_IGKV_fixed.fasta -locus IGK -align Macaca_mulatta_IGKV_gapped_fixed.fasta results.csv

python ../../../../src/find_alignments.py Macaca_mulatta_IGKV_fixed.fasta IMGT000063.fasta "rhemac10_IGKV_igk.csv" ../../../rhesus_motifs/IGK  -ref imgt,Macaca_mulatta_IGKV_fixed.fasta -locus IGK -align Macaca_mulatta_IGKV_gapped_fixed.fasta results_rev.csv -sense reverse

python ../../../../src/compare_annotations.py results.csv "IMGT genes.csv" forward rhemac10_IGK_digger_annotation_compared_to_imgt_fw.csv
python ../../../../src/compare_annotations.py results_rev.csv "IMGT genes.csv" forward rhemac10_IGK_digger_annotation_compared_to_imgt_rev.csv
