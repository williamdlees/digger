.. _parse_imgt_annotations_label:

parse_imgt_annotations
======================

``parse_imgt_annotations`` downloads an IMGT annotation file or or uses a file already downloaded. It parses the file to provide a list of annotated features.
Optionally it will also store the file downloaded, and create a FASTA file containing the annotated assembly.

.. argparse::
   :filename: ../src/digger/parse_imgt_annotations.py
   :func: get_parser
   :prog: parse_imgt_annotations

