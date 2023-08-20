mkdir work
cd work
extract_refs -L IGH "Homo sapiens"
cat Homo_sapiens_IGHV.fasta Homo_sapiens_IGHD.fasta Homo_sapiens_IGHJ.fasta > Homo_sapiens_IGHVDJ.fasta

parse_imgt_annotations --save_sequence IMGT000035.fasta "https://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000035" IMGT000035_genes.csv IGH
digger IMGT000035.fasta -v_ref Homo_sapiens_IGHV.fasta -d_ref Homo_sapiens_IGHD.fasta -j_ref Homo_sapiens_IGHJ.fasta -v_ref_gapped Homo_sapiens_IGHV_gapped.fasta -ref imgt,Homo_sapiens_IGHVDJ.fasta -species human IMGT000035.csv
compare_annotations IMGT000035.csv IMGT000035_genes.csv forward ../comparison_results
cd ..
