# Download IMGT reference sets for required species
# This script requires the python receptor_utils package: https://williamdlees.github.io/receptor_utils/_build/html/introduction.html
# It can be installed with pip: pip install receptor_utils. Biopython (https://biopython.org/) is required and must be installed first.

rm *.fasta
wget --no-check-certificate -O imgt_refs.fasta http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP

extract_refs -L IGH imgt_refs.fasta "Macaca mulatta"
cat Macaca_mulatta_IGHV.fasta Macaca_mulatta_IGHD.fasta Macaca_mulatta_IGHJ.fasta > Macaca_mulatta_IGHVDJ.fasta
python ../src/fix_macaque_gaps.py Macaca_mulatta_IGHV_gapped.fasta Macaca_mulatta_IGHV_gapped_fixed.fasta IGH

extract_refs -L IGK imgt_refs.fasta "Macaca mulatta"
cat Macaca_mulatta_IGKV.fasta Macaca_mulatta_IGKJ.fasta > Macaca_mulatta_IGKVJ.fasta
python ../src/fix_macaque_gaps.py Macaca_mulatta_IGKV_gapped.fasta Macaca_mulatta_IGKV_gapped_fixed.fasta IGK

extract_refs -L IGL imgt_refs.fasta "Macaca mulatta"
cat Macaca_mulatta_IGLV.fasta Macaca_mulatta_IGLJ.fasta > Macaca_mulatta_IGLVJ.fasta
python ../src/fix_macaque_gaps.py Macaca_mulatta_IGLV_gapped.fasta Macaca_mulatta_IGLV_gapped_fixed.fasta IGL

extract_refs -L IGH imgt_refs.fasta "Homo sapiens"
cat Homo_sapiens_IGHV.fasta Homo_sapiens_IGHD.fasta Homo_sapiens_IGHJ.fasta > Homo_sapiens_IGHVDJ.fasta

extract_refs -L IGK imgt_refs.fasta "Homo sapiens"
cat Homo_sapiens_IGKV.fasta Homo_sapiens_IGKJ.fasta > Homo_sapiens_IGKVJ.fasta

extract_refs -L IGL imgt_refs.fasta "Homo sapiens"
cat Homo_sapiens_IGLV.fasta Homo_sapiens_IGLJ.fasta > Homo_sapiens_IGLVJ.fasta

