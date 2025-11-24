mkdir work
cd work
extract_refs -L TRA "Homo sapiens"
cat Homo_sapiens_TRAV.fasta Homo_sapiens_TRAJ.fasta > Homo_sapiens_TRAVJ.fasta
# IMGT now block 'robots' from downloading annotations. See notes in parse_imgt_annotations documentation about downloading manually.
# parse_imgt_annotations --save_sequence IMGT000024.fasta "https://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000024" IMGT000024_genes.csv TRA


digger IMGT000024.fasta -v_ref Homo_sapiens_TRAV.fasta -locus TRA -j_ref Homo_sapiens_TRAJ.fasta -v_ref_gapped Homo_sapiens_TRAV_gapped.fasta -ref imgt,Homo_sapiens_TRAVJ.fasta -species Human IMGT000024.csv
compare_annotations -nc IMGT000024.csv IMGT000024_genes.csv forward ../comparison_results --target_locus TRA
cp IMGT000024.csv ..
cd ..
