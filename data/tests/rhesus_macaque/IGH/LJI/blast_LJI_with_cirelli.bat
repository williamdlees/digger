REM D genes determined by SEARCH-D
REM Because the cirelli assembly consists of multiple contigs, run_digger can't be used

mkdir work
cd work
cp ../cirelli_*.fasta .
cp ../../../../Macaca_mulatta_IGHVDJ.fasta .

cat cirelli_IGHV.fasta cirelli_IGHD.fasta cirelli_IGHJ.fasta >cirelli_IGH.fasta

makeblastdb -in cirelli_IGHV.fasta -dbtype nucl
makeblastdb -in cirelli_IGHD.fasta -dbtype nucl
makeblastdb -in cirelli_IGHJ.fasta -dbtype nucl

blastn -db cirelli_IGHV.fasta -query ../LJI_Rh.PacBio_All.IGH.Contigs.fasta -out LJI_IGHV.out -outfmt 7 -gapopen 5 -gapextend 5 -penalty -1 -word_size 11
blastn -db cirelli_IGHD.fasta -query ../LJI_Rh.PacBio_All.IGHDJ.Contigs.fasta -out LJI_IGHD.out -outfmt 7 -gapopen 5 -gapextend 5 -penalty -1 -word_size 7
blastn -db cirelli_IGHJ.fasta -query ../LJI_Rh.PacBio_All.IGHDJ.Contigs.fasta -out LJI_IGHJ.out -outfmt 7 -gapopen 5 -gapextend 5 -penalty -1 -word_size 7
python ../../../../../src/blastresults_to_csv.py LJI_IGHV.out LJI_IGHVJ_
python ../../../../../src/blastresults_to_csv.py LJI_IGHD.out LJI_IGHVJ_ -a
python ../../../../../src/blastresults_to_csv.py LJI_IGHJ.out LJI_IGHVJ_ -a

cd ../../../../../../SEARCH-D
python search_d.py "../digger/data/rhesus_macaque/IGH/LJI/LJI_Rh.PacBio_All.IGHDJ.Contigs.fasta" "../digger/data/rhesus_macaque/IGH/LJI/work/cirelli_matches_LJI_search-d.fasta"
cd ../digger/data/rhesus_macaque/IGH/LJI/work

python ../../../../../src/find_alignments.py cirelli_IGH.fasta ../LJI_Rh.PacBio_All.IGH.Contigs.fasta LJI_IGHVJ_*quiver.csv ../../../../rhesus_motifs/IGH -ref searchd,cirelli_matches_LJI_search-d.fasta -ref imgt,Macaca_mulatta_IGHVDJ.fasta -ref cirelli,cirelli_IGH.fasta cirelli_matches_LJI_IGHVDJ_blast.csv

python ../../../../../src/merge_search-d.py cirelli_matches_LJI_IGHVDJ_blast.csv cirelli_matches_LJI_search-d.fasta "+" 1 1051363 1051363 ../LJI_Rh.PacBio_All.IGHDJ.Contigs.fasta -ref imgt,Macaca_mulatta_IGHVDJ.fasta -ref cirelli,cirelli_IGH.fasta cirelli_matches_LJI_IGHVDJ.csv

cp cirelli_matches_LJI_IGHVDJ.csv ..