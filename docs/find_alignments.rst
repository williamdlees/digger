.. _find_alignments_label:

find_alignments
===============

``find_alignments`` takes the output of a blast search (as formatted by blastresults_to_csv.py), checks each identified location for the presence of a gene, and annotates if found.

.. argparse::
   :filename: ../src/digger/find_alignments.py
   :func: get_parser
   :prog: find_alignments

