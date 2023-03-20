mkdir work
cd work
extract_refs -L IGK "Homo sapiens"
cat Homo_sapiens_IGKV.fasta Homo_sapiens_IGKJ.fasta > Homo_sapiens_IGKVJ.fasta

parse_imgt_annotations --save_sequence IMGT000100.fasta "http://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000100" IMGT000100_genes.csv IGK
digger IMGT000100.fasta -v_ref Homo_sapiens_IGKV.fasta -j_ref Homo_sapiens_IGKJ.fasta -v_ref_gapped Homo_sapiens_IGKV_gapped.fasta -ref imgt,Homo_sapiens_IGKVJ.fasta -searchd ../../../../../../SEARCH-D ../../../../human_motifs/IGK IMGT000100.csv
compare_annotations IMGT000100.csv IMGT000100_genes.csv forward ../rhemac10_IGK_digger_annotation_compared_to_imgt_searchd
