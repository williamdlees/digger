.. _blastresults_to_csv:

blastresults_to_csv
===================

``blastresults_to_csv`` converts the output of a blast search in 'format 7' into csv format. If there were multiple query sequences, results are split into separate files.

.. argparse::
   :filename: ../src/digger/blastresults_to_csv.py
   :func: get_parser
   :prog: blastresults_to_csv

