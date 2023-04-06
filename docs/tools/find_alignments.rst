.. _find_alignments:

find_alignments
===============

``find_alignments`` takes the output of a blast search (as formatted by blastresults_to_csv.py), checks each identified location for the presence of a gene, and annotates if found.
Please refer to :ref:`rhesus_igh` for example usage of this and the other 'individual' commands.

.. argparse::
   :filename: ../src/digger/find_alignments.py
   :func: get_parser
   :prog: find_alignments

