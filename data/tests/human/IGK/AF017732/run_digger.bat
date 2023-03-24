mkdir work
cd work
extract_refs -L IGK "Homo sapiens"
cat Homo_sapiens_IGKV.fasta Homo_sapiens_IGKJ.fasta > Homo_sapiens_IGKVJ.fasta

parse_imgt_annotations --save_sequence AF017732.fasta "http://www.imgt.org/ligmdb/view?format=IMGT&id=AF017732" AF017732_genes.csv IGK
digger AF017732.fasta -locus IGK -sense forward -v_ref Homo_sapiens_IGKV.fasta -j_ref Homo_sapiens_IGKJ.fasta -v_ref_gapped Homo_sapiens_IGKV_gapped.fasta -ref imgt,Homo_sapiens_IGKVJ.fasta ../../../../human_motifs/IGK AF017732.csv
compare_annotations --filter_annot forward AF017732.csv AF017732_genes.csv forward ../rhemac10_IGK_digger_annotation_compared_to_imgt_forward
digger AF017732.fasta -locus IGK -sense reverse -v_ref Homo_sapiens_IGKV.fasta -j_ref Homo_sapiens_IGKJ.fasta -v_ref_gapped Homo_sapiens_IGKV_gapped.fasta -ref imgt,Homo_sapiens_IGKVJ.fasta ../../../../human_motifs/IGK AF017732_reverse.csv
compare_annotations --filter_annot reverse AF017732_reverse.csv AF017732_genes.csv forward ../rhemac10_IGK_digger_annotation_compared_to_imgt_reverse

cd ..