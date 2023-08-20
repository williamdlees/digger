mkdir work
cd work
extract_refs -L IGH "Macaca mulatta"
fix_macaque_gaps Macaca_mulatta_IGHV_gapped.fasta Macaca_mulatta_IGHV_gapped_fixed.fasta IGH
cat Macaca_mulatta_IGHV.fasta Macaca_mulatta_IGHD.fasta Macaca_mulatta_IGHJ.fasta > Macaca_mulatta_IGHVDJ.fasta

parse_imgt_annotations --save_sequence IMGT000064.fasta "https://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000064" IMGT000064_genes.csv IGH
digger IMGT000064.fasta -v_ref Macaca_mulatta_IGHV.fasta -d_ref Macaca_mulatta_IGHD.fasta -j_ref Macaca_mulatta_IGHJ.fasta -v_ref_gapped Macaca_mulatta_IGHV_gapped_fixed.fasta -ref imgt,Macaca_mulatta_IGHVDJ.fasta -species rhesus_macaque -locus IGH IMGT000064.csv
compare_annotations IMGT000064.csv IMGT000064_genes.csv forward ../rhemac10_IGH_digger_annotation_compared_to_imgt

cd ..