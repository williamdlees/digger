mkdir work
cd work
extract_refs -L IGH "Homo sapiens"
cat Homo_sapiens_IGHV.fasta Homo_sapiens_IGHD.fasta Homo_sapiens_IGHJ.fasta > Homo_sapiens_IGHVDJ.fasta

digger ../IMGT000035_248000_270000.fasta -keepwf -v_ref Homo_sapiens_IGHV.fasta -d_ref Homo_sapiens_IGHD.fasta -j_ref Homo_sapiens_IGHJ.fasta -v_ref_gapped Homo_sapiens_IGHV_gapped.fasta -ref imgt,Homo_sapiens_IGHVDJ.fasta -species human IMGT000035.csv

cp IMGT000035.csv ..

cd ..
