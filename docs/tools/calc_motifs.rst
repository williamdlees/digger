.. _calc_motifs:

calc_motifs
===========

``calc_motifs`` creates motif files for RSS and leader fields, based on the features produced by parse_imgt_annotations. Only annotations
of sequences marked as 'functional' are considered.

.. argparse::
   :filename: ../src/digger/calc_motifs.py
   :func: get_parser
   :prog: calc_motifs
