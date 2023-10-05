mkdir work
cd work
extract_refs -L TRD "Homo sapiens"
cat Homo_sapiens_TRDV.fasta Homo_sapiens_TRDJ.fasta > Homo_sapiens_TRDVJ.fasta

parse_imgt_annotations --save_sequence IMGT000024.fasta "https://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000024" IMGT000024_genes.csv TRD


digger IMGT000024.fasta -v_ref Homo_sapiens_TRDV.fasta -locus TRD -j_ref Homo_sapiens_TRDJ.fasta -d_ref Homo_sapiens_TRDD.fasta -v_ref_gapped Homo_sapiens_TRDV_gapped.fasta -ref imgt,Homo_sapiens_TRDVJ.fasta -species Human IMGT000024.csv
compare_annotations IMGT000024.csv IMGT000024_genes.csv forward --target_locus TRD ../comparison_results
cp IMGT000024.csv ..

cd ..
