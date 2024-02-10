mkdir work
cd work
extract_refs -L IGH "Macaca mulatta"
fix_macaque_gaps Macaca_mulatta_IGHV_gapped.fasta Macaca_mulatta_IGHV_gapped_fixed.fasta IGH
rm Macaca_mulatta_IGHV_gapped.fasta
cat Macaca_mulatta_IGHV.fasta Macaca_mulatta_IGHD.fasta Macaca_mulatta_IGHJ.fasta > Macaca_mulatta_IGHVDJ.fasta

parse_imgt_annotations --save_sequence IMGT000035.fasta "https://www.imgt.org/ligmdb/view?format=IMGT&id=BK063799" IMGT000035_genes.csv IGH
digger -keepwf IMGT000035.fasta -v_ref Macaca_mulatta_IGHV.fasta -d_ref Macaca_mulatta_IGHD.fasta -j_ref Macaca_mulatta_IGHJ.fasta -v_ref_gapped Macaca_mulatta_IGHV_gapped_fixed.fasta -ref imgt,Macaca_mulatta_IGHVDJ.fasta -species human IMGT000035.csv
compare_annotations IMGT000035.csv IMGT000035_genes.csv forward ../comparison_results
cp IMGT000035.csv ..
cd ..
