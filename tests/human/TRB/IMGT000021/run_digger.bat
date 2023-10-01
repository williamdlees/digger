mkdir work
cd work
extract_refs -L TRB "Homo sapiens"
cat Homo_sapiens_TRBV.fasta Homo_sapiens_TRBD.fasta Homo_sapiens_TRBJ.fasta > Homo_sapiens_TRBVDJ.fasta

parse_imgt_annotations --save_sequence IMGT000021.fasta "https://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000021" IMGT000021_genes.csv TRB


digger IMGT000021.fasta -v_ref Homo_sapiens_TRBV.fasta -locus TRB -d_ref Homo_sapiens_TRBD.fasta -j_ref Homo_sapiens_TRBJ.fasta -v_ref_gapped Homo_sapiens_TRBV_gapped.fasta -ref imgt,Homo_sapiens_TRBVDJ.fasta -species Human IMGT000021.csv
compare_annotations IMGT000021.csv IMGT000021_genes.csv forward ../comparison_results
cd ..
