.. _annotation:

Anotation format
================

This page describes the annotation file produced by digger / find_alignments

Columns in the Annotation File
******************************

In addition to the columns in the first table, the file contains the columns in the second table, prefixed by the reference name, for each reference specified with a -ref argument.

.. csv-table::
   :file: annotation_cols.tsv
   :delim: tab
   :header-rows: 1
   :widths: 25, 75
      
Columns provided for each -ref:
   
.. csv-table::
   :file: annotation_ref_cols.tsv
   :delim: tab
   :header-rows: 1
   :widths: 25, 75
   

Functionality
*************

Functionality is assigned as follows:

Functional

* RSS and leader meet or exceed position-weighted matrix threshold
* Highly-conserved nucleotides agree with the definition for the locus, if a definition has been specified
* If a V-gene, leader starts with ATG, donor splice ends GT or CT, acceptor splice ends AG, and spliced leader has no stop codons
* If a V-gene, coding region has no stop codons before the cysteine at IMGT position 104
* If a V-gene, conserved nucleotides are at the expected locations
* If a J-gene, donor splice is as expected and coding region has no stop codons

ORF

* One or more of the above conditions are not met, but no stop codon has been detected
* If a V-gene, leader starts with ATG

Pseudo

* Coding region contains stop codon(s)
* Leader does not start with ATG

.. _leader_annotation:

V Leader Annotation
*******************

Exons 1 and 2 are annotated in accordance with the customary genomic annotation regarding splice sites. At first glance, it may be thought that the L-PART1 sequence
should be the same as the EXON1 sequence, and the L-PART2 sequence the same as the 5' end of the EXON2 sequence. However IMGT define L-PART1 and L-PART2 in a 
manner that requires each to occupy an entire number of codons. To achieve this, any extra nucleotides at the end of EXON1 that do not make up a complete codon 
are assigned to L-PART2. The approach is documented `here <https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/splicing>`_. In the table on that page, the first 
row, in which a single nucleotide is transferred to L-PART2, is the case normally encountered in V-genes.

It should be clear from the above that L-PART1 and L-PART2 are not usually separated by the V-INTRON, as is often depicted in the literature.

We provide coordinates for EXON1 and EXON2 but not for L-PART1 and L-PART2, as the latter do not have simple 'start' and 'end' coordinates in the genomic sequence.
They are best thought of as protein-based features.