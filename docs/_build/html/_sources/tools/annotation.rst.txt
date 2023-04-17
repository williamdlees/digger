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
* If a V-gene, leader starts with ATG, and spliced leader has no stop codons
* If a V-gene, coding region has no stop codons before the cysteine at IMGT position 104
* If a V-gene, conserved nucleotides are at the expected locations
* If a J-gene, donor splice is as expected and coding region has no stop codons

ORF

* One or more of the above conditions are not met, but no stop codon has been detected
* If a V-gene, leader starts with ATG

Pseudo

* Coding region contains stop codon(s)
* Leader does not start with ATG

