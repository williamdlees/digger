.. human_igh:

Annotating the rhesus macaque IGH locus
=======================================

In this example we will see how to:
- How to create PWMs from existing IMGT annotations
- Use the underlying commands that digger calls, and how they might be useful when annotating multiple contigs or scaffolds

IMGT has identified scaffolds in the 2006 rhesus macaque reference assembly, Mmul_051212, which lie within the IGH locus. Here we will bring them together in a single file and annotate them with PWMs derived from the current reference assembly, rhemac10 (Mmul_10).
While this example is somewhat artificial, in that the scaffolds could equally well be handled individually using the digger command, the approach is useful where the number of sequences to be processed
is large. It also serves to show how the individual commands in the package can be used. This provides some additional flexibility, for example in tuning the blast searches, and also illsutates how they work together.

Data
****

As in the previous example, the rhesus IGH reference set can be downloaded from IMGT with the `receptor_utils <https://williamdlees.github.io/receptor_utils/_build/html/introduction.html>`__ command 
``extract_refs`` (receptor_utils is installed as part of digger's installation). However, the rhesus IG gapped V-genes provided by IMGT contain additional inserted codons relative to
the conventional IMGT alignment. As digger (alongside other tools) expects the conventional alignment, a further step is needed to realign the gapped sequences, using the receptor_utils 
tool ``fix_macaque_gaps``. The following commands will prepare the reference data::

   > extract_refs -L IGH "Macaca mulatta"
   > fix_macaque_gaps Macaca_mulatta_IGHV_gapped.fasta \
       Macaca_mulatta_IGHV_gapped_fixed.fasta IGH
   > cat Macaca_mulatta_IGHV.fasta Macaca_mulatta_IGHD.fasta Macaca_mulatta_IGHJ.fasta \
       > Macaca_mulatta_IGHVDJ.fasta
	   
The Mmul_51212 can be downloaded from IMGT as follows::

   > parse_imgt_annotations --save_sequence NW_001157919.fasta \
      "https://www.imgt.org/ligmdb/view.action?format=IMGT&id=NW_001157919" \
	  NW_001157919_genes.csv IGH
   > parse_imgt_annotations --save_sequence NW_001122023.fasta \
      "https://www.imgt.org/ligmdb/view.action?format=IMGT&id=NW_001122023" \
      NW_001122023_genes.csv IGH
   > parse_imgt_annotations --save_sequence NW_001122024.fasta \
      "https://www.imgt.org/ligmdb/view.action?format=IMGT&id=NW_001122024" \
	  NW_001122024_genes.csv IGH
   > parse_imgt_annotations --save_sequence NW_001121239.fasta \
      "https://www.imgt.org/ligmdb/view.action?format=IMGT&id=NW_001121239" \
	  NW_001121239_genes.csv IGH
   > parse_imgt_annotations --save_sequence NW_001121240.fasta \
      "https://www.imgt.org/ligmdb/view.action?format=IMGT&id=NW_001121240" \
	  NW_001121240_genes.csv IGH
	 
   >cat NW_001157919.fasta NW_001122023.fasta NW_001122024.fasta \
	  NW_001121239.fasta NW_001121240.fasta > Mmul_051212.fasta



Preparing position-weighted matrices
************************************

Digger already has PWMs for rhesus IGH, but for the purpose of this example, we will create a set using the features listed in IMGT's annotation of the rhemac10 IGH locus, which 
has the IMGT accession number IMGT000064. This is how digger's built-in rhesus PWMs were created. The following commands download the annotation, determine the features, and calculate the PWMs from 
features of functional annotations::

   > mkdir motifs
   > cd motifs
   > parse_imgt_annotations \
	   "http://www.imgt.org/ligmdb/view?format=IMGT&id=IMGT000064" \
	   IMGT000064_genes.csv IGH
   > calc_motifs IMGT000064_genes.csv
   
``calc_motifs`` will create 10 motif files in the directory.

The motifs directory may optionally contain a FASTA file ``conserved_motifs.fasta`` defining strongly-conserved nucleotides in the RSS and leader. Only those features 
with conserved residues need to be listed in the file. The names follow the filenames used for the PWMs.
The following sequences were derived from Figure 3 of Ngoune et al. (2022) and will be used in this example::

   >V-HEPTAMER
   CAC---G
   >V-NONAMER
   -----AACC
   >5'D-HEPTAMER
   ----GTG
   >5'D-NONAMER
   ---T-----
   >3'D-HEPTAMER
   C-C---G
   >3'D-NONAMER
   -C----A--
   >J-HEPTAMER
   C--TGTG
   >J-NONAMER
   -GTT--TG-
   
Again, this file is provided for download at the location provided near the top of this example.
   
While the presence or absence of conserved residues can be a useful guide to the likely functionality of a sequence, please bear in mind that it is a guide only:
exceptions can be expected, particularly where the definitions have been built on limited data.

Annotating the Assembly
***********************

The digger command is not able to handle a FASTA file containing multiple contigs, so we will call the underlying tools directly. We start by creating the blast databases and querying against the assembly, 
using the reference genes determined in the study::

   > makeblastdb -in Macaca_mulatta_IGHV.fasta -dbtype nucl
   > makeblastdb -in Macaca_mulatta_IGHD.fasta -dbtype nucl
   > makeblastdb -in Macaca_mulatta_IGHJ.fasta -dbtype nucl

   > blastn -db Macaca_mulatta_IGHV.fasta -query Mmul_051212.fasta -out mmul_IGHV.out \
      -outfmt 7 -gapopen 5 -gapextend 5 -penalty -1 -word_size 11
   > blastn -db Macaca_mulatta_IGHD.fasta -query Mmul_051212.fasta -out mmul_IGHD.out \
      -outfmt 7 -gapopen 5 -gapextend 5 -penalty -1 -word_size 7 -evalue 100
   > blastn -db Macaca_mulatta_IGHJ.fasta -query Mmul_051212.fasta -out mmul_IGHJ.out \
      -outfmt 7 -gapopen 5 -gapextend 5 -penalty -1 -word_size 7


Note that a higher evalue is used for the D genes, as they can be quite short.

Next we call ``blastresults_to_csv`` to convert to a more convenient format::

    > blastresults_to_csv mmul_IGHV.out mmul_ighvdj_   
    > blastresults_to_csv mmul_IGHD.out mmul_ighvdj_ -a
    > blastresults_to_csv mmul_IGHJ.out mmul_ighvdj_ -a

The commands instruct the tool to create merged files containing V,D and J hits. This is achieved by specifying the same prefix on each command ``(mmul_ighvdj_)`` and using the ``-a`` (append) option.
The records created by blastn contain the name of the contig in which a hit was found. ``blastresults_to_csv`` will create one file per contig. The names contain the ID of the contig in 
``Mmul_051212.fasta``, except that they are modified where necessary to ensure file system compatibility.

We now call find_alignments to process the annotations::

    > find_alignments Macaca_mulatta_IGHVDJ.fasta \
	   Mmul_051212.fasta \
	   "mmul_ighvdj_nw_*.csv" \
	   -ref imgt,Macaca_mulatta_IGHVDJ.fasta \
	   -align Macaca_mulatta_IGHV_gapped_fixed.fasta \
	   -motif_dir motifs \
	   Mmul_051212.csv

Note that the third argument, ``"mmul_ighvdj_nw_*.csv"``, contains a wildcard that will match all the files produced in the previous step. It is quoted to avoid expansion by the shell. 
V-genes in the annotation will be annotated and gapped using the IMGT set as a template (with fixed gaps).
``find_alignments`` will attempt to deduce the sense in which to annotate each segment. This is helpful in this case as the contigs vary in their orientation.  Note that we are
specifying the location of the motifs directory created in the previous step rather than the species and locus, which would cause digger to use the built-in tables.


Comparing the output to the study's annotation
**********************************************

``compare_annotations`` is not capable of handling the output from multiple sequences in the same file, so unfortunately we need to split the results up for the comparison:

    > head -n 1 Mmul_051212.csv > mmul_header.csv

    > cp mmul_header.csv NW_001157919_digger.csv
    > grep NW_001157919 Mmul_051212.csv >> NW_001157919_digger.csv

    > cp mmul_header.csv NW_001122023_digger.csv
    > grep NW_001122023 Mmul_051212.csv >> NW_001122023_digger.csv

    > cp mmul_header.csv NW_001122024_digger.csv
    > grep NW_001122024 Mmul_051212.csv >> NW_001122024_digger.csv

    > cp mmul_header.csv NW_001121239_digger.csv
    > grep NW_001121239 Mmul_051212.csv >> NW_001121239_digger.csv

    > cp mmul_header.csv NW_001121240_digger.csv
    > grep NW_001121240 Mmul_051212.csv >> NW_001121240_digger.csv

    > compare_annotations NW_001157919_digger.csv NW_001157919_genes.csv forward NW_001157919_comp
    > compare_annotations NW_001122023_digger.csv NW_001122023_genes.csv forward NW_001122023_comp
    > compare_annotations NW_001122024_digger.csv NW_001122024_genes.csv forward NW_001122024_comp
    > compare_annotations NW_001121239_digger.csv NW_001121239_genes.csv forward NW_001121239_comp
    > compare_annotations NW_001121240_digger.csv NW_001121240_genes.csv forward NW_001121240_comp



The output file ``cirelli_table_s1_with_digger_hits.csv`` recapitulates the relevant data from the study table S1, with an additional column showing whether and where where the sequence was found in the digger annotation.
It shows that all genes identified in the study's annotation of the IGH locus using the 'primary' method of annotation were listed also by digger, with the exception of one D gene, LJI.Rh_IGHD4.22.
This gene sequence was identified by BLAST in the digger run, but was not listed in the file as neither 3' nor 5' RSS passed the PWM threshold and were therefore both regarded as invalid. As D gene
sequences are short, a sequence match with invalid RSS at each end is a frequent false positive. 

Of the remaining 112 genes, which were all classified as F/ORF in the study, digger categorised 97 as functional, 13 as ORF and 2 as pseudogenes.  Among the functional and ORF genes, the length of the
coding sequence assigned by digger differed from that in Table S1 on 12 occasions, indicating differences in the identification of the RSS.
Of the two that digger classed as pseudogenes,
LJI.Rh_IGHV1.138 is noted as 'Leader missing initial ATG, Stop codon in leader' and LJI.Rh_IGHV3.107 as 'Stop codon in leader, First cysteine not found'. Interestingly, LJI.Rh_IGHV1.138
was observed in repertoires during the course of the study, suggesting either an error in digger's annotation, or a sequencing error. 

The output file ``digger_F_ORF_with_cirelli_annots.csv`` lists 3 functional V-genes and 7 functional D-genes identified by digger but not in Table S1. Three of the D-genes and two of the V-genes have exact
matches in the IMGT reference set.

Overall, the results are in good agreement, with, nevertheless, some interesting points of detail that merit further examination.

References
**********

Cirelli et al., 2019, Slow Delivery Immunization Enhances HIV Neutralizing Antibody and Germinal Center Responses via Modulation of Immunodominance. *Cell* `doi: 10.1016/j.cell.2019.04.012 <https://doi.org/10.1016/j.cell.2019.04.012>`__.

Ngoune et al., 2022, IMGT® Biocuration and Analysis of the Rhesus Monkey IG Loci. *Vaccines* `doi: 10.3390/vaccines10030394 <https://www.mdpi.com/2076-393X/10/3/394#>`__.