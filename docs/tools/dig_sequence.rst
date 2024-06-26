.. _dig_sequence:

dig_sequence
============

``dig_sequence`` provides targeted annotation of genomic sequences. The sequences can be as small as a single coding region,
or as large as an entire locus. The tool will search for the sequence that best matches a specified target sequence, and
annotate just that best match.

Options are available to annotate a single sequence in a FASTA file, a single genbank ID (the sequence will be fetched from
Genbank), or a list of sequences or Genbank IDs specified in a CSV file. For GenBank requests, an email address must be
provided, as this is a requirement of the GenBank API.

Example usage is described in :ref:`targeted_annotation`.

.. argparse::
   :filename: ../src/digger/dig_sequence.py
   :func: get_parser
   :prog: dig_sequence
