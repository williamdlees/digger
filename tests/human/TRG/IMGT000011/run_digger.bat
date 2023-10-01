mkdir work
cd work
extract_refs -L TRG "Homo sapiens"
cat Homo_sapiens_TRGV.fasta Homo_sapiens_TRGJ.fasta > Homo_sapiens_TRGVJ.fasta

parse_imgt_annotations --save_sequence IMGT000011.fasta "https://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000011" IMGT000011_genes.csv TRG


digger IMGT000011.fasta -v_ref Homo_sapiens_TRGV.fasta -locus TRG -j_ref Homo_sapiens_TRGJ.fasta -v_ref_gapped Homo_sapiens_TRGV_gapped.fasta -ref imgt,Homo_sapiens_TRGVJ.fasta -species Human IMGT000011.csv
compare_annotations IMGT000011.csv IMGT000011_genes.csv forward ../comparison_results
cd ..
