mkdir work
cd work
extract_refs -L IGK "Homo sapiens"
cat Homo_sapiens_IGKV.fasta Homo_sapiens_IGKJ.fasta > Homo_sapiens_IGKVJ.fasta

digger ../NC_000002.12.fasta -locus IGK -sense forward -v_ref Homo_sapiens_IGKV.fasta -j_ref Homo_sapiens_IGKJ.fasta -v_ref_gapped Homo_sapiens_IGKV_gapped.fasta -ref imgt,Homo_sapiens_IGKVJ.fasta -species human NC_000002.12.csv
digger ../NC_000002.12.fasta -locus IGK -sense reverse -v_ref Homo_sapiens_IGKV.fasta -j_ref Homo_sapiens_IGKJ.fasta -v_ref_gapped Homo_sapiens_IGKV_gapped.fasta -ref imgt,Homo_sapiens_IGKVJ.fasta -species human NC_000002.12_reverse.csv

compare_annotations NC_000002.12.csv ../watson_features.csv forward ../comparison_results_rev --filter forward --comp_name "Watson et al."
compare_annotations NC_000002.12_reverse.csv ../watson_features.csv forward ../comparison_results_fw --filter reverse --comp_name "Watson et al."
cp NC_000002.12.csv  ..
cp NC_000002.12_reverse.csv  ..
cd ..