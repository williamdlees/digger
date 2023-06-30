.. _targeted_annotation:

Targeted Annotation
===================

In the examples covered so far, :ref:`digger` has been used to identify as many receptor genes as possible, 
using wide-ranging BLAST searches. A companion tool, :ref:`dig_sequence`, can be used to identify and annotate the closest match to a single specified
sequence. 

As an example of its use, an online `BLAST search <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`__ identifies a 100% sequence match to the human receptor gene IGHV1-18*04 at 
Genbank accession number `KC713938 <https://www.ncbi.nlm.nih.gov/nucleotide/KC713938.1?report%253Dgenbank>`__. This can be annotated by the ``dig_sequence single`` command:

.. code-block:: console

    >dig_sequence single -align Homo_sapiens_IGHV_gapped.fasta -species human IGHV1-18*01 Homo_sapiens_IGHV.fasta KC713938
    Using motif files from C:\Users\William\miniconda3\envs\digby_genomics310\lib\site-packages\digger\motifs\human\IGH
    target_allele: IGHV1-18*01
    genbank_acc: KC713938
    genbank_seq: length: 935nt
    gene_seq: ATGGACTGGACCTGGAGCATCCTTTTCTTGGTGGCAGCAGCAACAGGTAACGGACTCCCCAGTCCCAGGGCTGAGAGAGAAACCAGGCCAGTCATGTGAGACTTCACCCACTCCTGTGTCCTCTCCACAGGTGCCCACTCCCAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTACGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGACACAGTGTGAAAACCCACATCCTGAGGGTTTCAGAAACC
    seq: CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTACGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA
    alignment_score: 99.7
    nt_diff: 1
    snps: _t111c
    start: 392
    end: 687
    sense: +
    functional: Functional
    notes:
    l_part1: ATGGACTGGACCTGGAGCATCCTTTTCTTGGTGGCAGCAGCAACAG
    l_part2: GTGCCCACTCC
    v_heptamer: CACAGTG
    v_nonamer: TCAGAAACC

The command takes as arguments the id of the target sequence to search for (in this case IGHV1-18*01), the fasta file in
which the sequence can be found, and the Genbank ID to search. Optional arguments include the species (otherwise human 
is assumed) and a file of gapped reference sequences, which is used as a guide to gap V sequences. If an output file 
is specified 
with the ``-out_file`` argument, the output is written there in CSV format, with information matching that provided by 
``digger``. Otherwise a summary is provided to standard output. 

The command searches the specified Genbank accession for the closest match to the target sequenece (in this case
differing by a single nucleotide), and returns annotation details. These will include regulatory regions where they are 
available. Only a single closest-match sequence is annotated.

Variants of :ref:`dig_sequence` allow the sequence to be specified directly rather than by Genbank accession number, and
allow multiple searches to be specified in a csv file.