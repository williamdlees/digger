mkdir work
cd work
extract_refs -L IGK "Macaca mulatta"
fix_macaque_gaps Macaca_mulatta_IGKV_gapped.fasta Macaca_mulatta_IGKV_gapped_fixed.fasta IGK
cat Macaca_mulatta_IGKV.fasta Macaca_mulatta_IGKJ.fasta > Macaca_mulatta_IGKVJ.fasta

parse_imgt_annotations --save_sequence IMGT000063.fasta "http://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000063" IMGT000063_genes.csv IGK

digger IMGT000063.fasta -locus IGK -v_ref Macaca_mulatta_IGKV.fasta -v_ref_gapped Macaca_mulatta_IGKV_gapped_fixed.fasta -j_ref Macaca_mulatta_IGKJ.fasta -ref imgt,Macaca_mulatta_IGKVJ.fasta -sense forward -species rhesus_macaque -locus IGK IMGT000063_fw.csv
digger IMGT000063.fasta -locus IGK -v_ref Macaca_mulatta_IGKV.fasta -v_ref_gapped Macaca_mulatta_IGKV_gapped_fixed.fasta -j_ref Macaca_mulatta_IGKJ.fasta -ref imgt,Macaca_mulatta_IGKVJ.fasta -sense reverse -species rhesus_macaque -locus IGK IMGT000063_rev.csv


compare_annotations --filter_annot forward IMGT000063_fw.csv IMGT000063_genes.csv forward ../IMGT000063_digger_annotation_compared_to_imgt_fw.csv
compare_annotations --filter_annot reverse IMGT000063_rev.csv IMGT000063_genes.csv forward ../IMGT000063_digger_annotation_compared_to_imgt_rev.csv
