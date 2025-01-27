.. _calc_motifs:

calc_motifs
===========

``calc_motifs`` creates motif files for RSS and leader fields, based on the features produced by :ref:`parse_imgt_annotations`. Only annotations
of sequences marked as 'functional' are considered. 

``calc_motifs`` also creates a file containing definitions of other motif parameters for the locus. These are copied from the values used for
human locus, but may be modified if necessary. The following values are included:

.. csv-table::
   :file: motif_params.tsv
   :delim: tab
   :header-rows: 1
   :widths: 25, 75


Please refer to :ref:`rhesus_igh` for example usage of this and the other 'individual' commands.

When calc_motifs runs, it provides a summary of the sequences it has processed. Please check these for sanity: for example check that the consensus 
hpetamer is 7nt long, and so on. Errors in the motif table will impact annotation.

.. argparse::
   :filename: ../src/digger/calc_motifs.py
   :func: get_parser
   :prog: calc_motifs
