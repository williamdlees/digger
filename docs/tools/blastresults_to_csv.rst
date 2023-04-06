.. _blastresults_to_csv:

blastresults_to_csv
===================

``blastresults_to_csv`` converts the output of a blast search in 'format 7' into csv format. If there were multiple query sequences, results are split into separate files.
Please refer to :ref:`rhesus_igh` for example usage of this and the other 'individual' commands.

.. argparse::
   :filename: ../src/digger/blastresults_to_csv.py
   :func: get_parser
   :prog: blastresults_to_csv

