mkdir work
cd work
extract_refs -L IGH "Homo sapiens"
cat Homo_sapiens_IGHV.fasta Homo_sapiens_IGHD.fasta Homo_sapiens_IGHJ.fasta > Homo_sapiens_IGHVDJ.fasta

parse_imgt_annotations --save_sequence IMGT000035.fasta "https://www.imgt.org/ligmdb/view?format=IMGT&id=BK063799" IMGT000035_genes.csv IGH

head -n 1 IMGT000035.fasta >IMGT000035_trunc.fasta
at_coords IMGT000035.fasta 1 15000 >>IMGT000035_trunc.fasta

digger IMGT000035_trunc.fasta -keepwf -v_ref Homo_sapiens_IGHV.fasta -d_ref Homo_sapiens_IGHD.fasta -j_ref Homo_sapiens_IGHJ.fasta -v_ref_gapped Homo_sapiens_IGHV_gapped.fasta -ref imgt,Homo_sapiens_IGHVDJ.fasta -species human IMGT000035.csv

cp IMGT000035_trunc.csv ..

cd ..
