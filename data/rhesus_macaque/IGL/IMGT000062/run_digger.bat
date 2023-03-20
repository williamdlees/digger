python ../../../../src/digger.py IMGT000062.fasta -locus IGL -v_ref Macaca_mulatta_IGLV_fixed.fasta -v_ref_gapped Macaca_mulatta_IGLV_gapped_fixed.fasta -j_ref Macaca_mulatta_IGLJ.fasta -ref imgt,Macaca_mulatta_IGLVJ.fasta -sense forward ../../../rhesus_motifs/IGL IMGT000062.csv
python ../../../../src/digger.py IMGT000062.fasta -locus IGL -v_ref Macaca_mulatta_IGLV_fixed.fasta -v_ref_gapped Macaca_mulatta_IGLV_gapped_fixed.fasta -j_ref Macaca_mulatta_IGLJ.fasta -ref imgt,Macaca_mulatta_IGLVJ.fasta -sense reverse ../../../rhesus_motifs/IGL IMGT000062_rev.csv


python ../../../../src/compare_annotations.py IMGT000062.csv "IMGT genes.csv" forward rhemac10_IGL_digger_annotation_compared_to_imgt_fw.csv
python ../../../../src/compare_annotations.py IMGT000062_rev.csv "IMGT genes.csv" reverse rhemac10_IGL_digger_annotation_compared_to_imgt_rev.csv
