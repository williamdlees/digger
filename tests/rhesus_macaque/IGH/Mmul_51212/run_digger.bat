mkdir work
cd work

extract_refs -L IGH "Macaca mulatta"
fix_macaque_gaps Macaca_mulatta_IGHV_gapped.fasta Macaca_mulatta_IGHV_gapped_fixed.fasta IGH 
cat Macaca_mulatta_IGHV.fasta Macaca_mulatta_IGHD.fasta Macaca_mulatta_IGHJ.fasta >Macaca_mulatta_IGHVDJ.fasta


parse_imgt_annotations --save_sequence NW_001157919.fasta "https://www.imgt.org/ligmdb/view.action?format=IMGT&id=NW_001157919" NW_001157919_genes.csv IGH
parse_imgt_annotations --save_sequence NW_001122023.fasta "https://www.imgt.org/ligmdb/view.action?format=IMGT&id=NW_001122023" NW_001122023_genes.csv IGH
parse_imgt_annotations --save_sequence NW_001122024.fasta "https://www.imgt.org/ligmdb/view.action?format=IMGT&id=NW_001122024" NW_001122024_genes.csv IGH
parse_imgt_annotations --save_sequence NW_001121239.fasta "https://www.imgt.org/ligmdb/view.action?format=IMGT&id=NW_001121239" NW_001121239_genes.csv IGH
parse_imgt_annotations --save_sequence NW_001121240.fasta "https://www.imgt.org/ligmdb/view.action?format=IMGT&id=NW_001121240" NW_001121240_genes.csv IGH

cat NW_001157919.fasta NW_001122023.fasta NW_001122024.fasta NW_001121239.fasta NW_001121240.fasta >Mmul_051212.fasta

mkdir motifs
cd motifs
parse_imgt_annotations "http://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000064" IMGT000064_genes.csv IGH
calc_motifs IMGT000064_genes.csv
cp ../../conserved_motifs.fasta .
cd ..

makeblastdb -in Macaca_mulatta_IGHV.fasta -dbtype nucl
makeblastdb -in Macaca_mulatta_IGHD.fasta -dbtype nucl
makeblastdb -in Macaca_mulatta_IGHJ.fasta -dbtype nucl

blastn -db Macaca_mulatta_IGHV.fasta -query Mmul_051212.fasta -out mmul_IGHV.out -outfmt 7 -gapopen 5 -gapextend 5 -penalty -1 -word_size 11
blastn -db Macaca_mulatta_IGHD.fasta -query Mmul_051212.fasta -out mmul_IGHD.out -outfmt 7 -gapopen 5 -gapextend 5 -penalty -1 -word_size 7 -evalue 100
blastn -db Macaca_mulatta_IGHJ.fasta -query Mmul_051212.fasta -out mmul_IGHJ.out -outfmt 7 -gapopen 5 -gapextend 5 -penalty -1 -word_size 7

blastresults_to_csv mmul_IGHV.out mmul_ighvdj_   
blastresults_to_csv mmul_IGHD.out mmul_ighvdj_ -a
blastresults_to_csv mmul_IGHJ.out mmul_ighvdj_ -a

find_alignments Macaca_mulatta_IGHVDJ.fasta Mmul_051212.fasta "mmul_ighvdj_nw_*.csv" -ref imgt,Macaca_mulatta_IGHVDJ.fasta -align Macaca_mulatta_IGHV_gapped_fixed.fasta -motif_dir motifs Mmul_051212.csv

head -n 1 Mmul_051212.csv > mmul_header.csv

cp mmul_header.csv NW_001157919_digger.csv
grep NW_001157919 Mmul_051212.csv >> NW_001157919_digger.csv

cp mmul_header.csv NW_001122023_digger.csv
grep NW_001122023 Mmul_051212.csv >> NW_001122023_digger.csv

cp mmul_header.csv NW_001122024_digger.csv
grep NW_001122024 Mmul_051212.csv >> NW_001122024_digger.csv

cp mmul_header.csv NW_001121239_digger.csv
grep NW_001121239 Mmul_051212.csv >> NW_001121239_digger.csv

cp mmul_header.csv NW_001121240_digger.csv
grep NW_001121240 Mmul_051212.csv >> NW_001121240_digger.csv

compare_annotations NW_001157919_digger.csv NW_001157919_genes.csv forward ../NW_001157919_comp
compare_annotations NW_001122023_digger.csv NW_001122023_genes.csv forward ../NW_001122023_comp
compare_annotations NW_001122024_digger.csv NW_001122024_genes.csv forward ../NW_001122024_comp
compare_annotations NW_001121239_digger.csv NW_001121239_genes.csv forward ../NW_001121239_comp
compare_annotations NW_001121240_digger.csv NW_001121240_genes.csv forward ../NW_001121240_comp

cd ..