.. _calc_motifs:

calc_motifs
===========

``calc_motifs`` creates motif files for RSS and leader fields, based on the features produced by :ref:`parse_imgt_annotations`. Only annotations
of sequences marked as 'functional' are considered. Please refer to :ref:`rhesus_igh` for example usage of this and the other 'individual' commands.

.. argparse::
   :filename: ../src/digger/calc_motifs.py
   :func: get_parser
   :prog: calc_motifs
