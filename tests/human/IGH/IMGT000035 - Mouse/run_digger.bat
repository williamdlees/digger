mkdir work
cd work
extract_refs -L IGH "Mus musculus"
cat Mus_musculus_IGHV.fasta Mus_musculus_IGHD.fasta Mus_musculus_IGHJ.fasta > Mus_musculus_IGHVDJ.fasta

parse_imgt_annotations --save_sequence IMGT000035.fasta "https://www.imgt.org/ligmdb/view?format=IMGT&id=BK063799" IMGT000035_genes.csv IGH

digger -keepwf IMGT000035.fasta -v_ref Mus_musculus_IGHV.fasta -d_ref Mus_musculus_IGHD.fasta -j_ref Mus_musculus_IGHJ.fasta -v_ref_gapped Mus_musculus_IGHV_gapped.fasta -ref imgt,Mus_musculus_IGHVDJ.fasta -species human IMGT000035.csv

compare_annotations IMGT000035.csv IMGT000035_genes.csv forward ../comparison_results
cp IMGT000035.csv ..
cd ..
