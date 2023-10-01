mkdir work
cd work
extract_refs -L IGL "Homo sapiens"
cat Homo_sapiens_IGLV.fasta Homo_sapiens_IGLJ.fasta > Homo_sapiens_IGLVJ.fasta

digger ../NC_000002.11.fasta -keepwf -locus IGL -sense forward -v_ref Homo_sapiens_IGLV.fasta -j_ref Homo_sapiens_IGLJ.fasta -v_ref_gapped Homo_sapiens_IGLV_gapped.fasta -ref imgt,Homo_sapiens_IGLVJ.fasta -species human NC_000002.11.csv
compare_annotations NC_000002.11.csv ../watson_features.csv forward ../comparison_results --filter forward --comp_name "Watson et al."

cd ..