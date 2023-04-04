mkdir work
cd work
extract_refs -L IGL "Macaca mulatta"
fix_macaque_gaps Macaca_mulatta_IGLV_gapped.fasta Macaca_mulatta_IGLV_gapped_fixed.fasta IGL
cat Macaca_mulatta_IGLV.fasta Macaca_mulatta_IGLJ.fasta > Macaca_mulatta_IGLVJ.fasta

parse_imgt_annotations --save_sequence IMGT000062.fasta "http://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000062" IMGT000062_genes.csv IGK

digger IMGT000062.fasta -locus IGL -v_ref Macaca_mulatta_IGLV.fasta -v_ref_gapped Macaca_mulatta_IGLV_gapped_fixed.fasta -j_ref Macaca_mulatta_IGLJ.fasta -ref imgt,Macaca_mulatta_IGLVJ.fasta -sense forward -species rhesus_macaque IMGT000062.csv


compare_annotations IMGT000062.csv IMGT000062_genes.csv forward ../rhemac10_IGL_digger_annotation_compared_to_imgt.csv
